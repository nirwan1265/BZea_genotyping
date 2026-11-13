#!/bin/bash
#SBATCH --job-name=filter_sites
#SBATCH --array=1-10
#SBATCH --output=logs/11_filter_sites.%A_%a.out
#SBATCH --error=logs/11_filter_sites.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=12:00:00
#SBATCH --mem=16G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# 11 — Site filtering + teosinte panel, one job per chromosome (bcftools + awk).
#
#   Keeps only the pools we trust: the 10 clean B73 controls and the donor pools
#   not in sample_lists/exclude_all.txt (see QC_NOTES.md).
#
#   A site is removed if:
#     - the clean B73 controls have > MAX_CTRL_ALT ALT reads (repeat / paralog / reference error)
#     - total depth > MAX_DP_MULT x chromosome median, or < MIN_DP_MULT x median
#     - more than MAX_OTH_FRAC of its reads show a third base
#     - the teosinte panel has the position but with different alleles
#   Remaining sites are split by whether the teosinte panel (227 plants) has the
#   same SNP; panel sites get the panel's ALT count (PANEL_AC of PANEL_AN chromosomes).
#
#   Output in ${D11}, per chromosome:
#     BZea_<chr>.PASS.bcf            all sites that passed, trusted pools only
#     BZea_<chr>.PASS.panel.bcf      ... that are also in the teosinte panel (+ INFO/PANEL_AC, PANEL_AN)
#     BZea_<chr>.PASS.nonpanel.bcf   ... that are not in the panel
#     BZea_<chr>.pool_reads.tsv      per pool: B73/ALT reads on panel and non-panel PASS sites
#     BZea_<chr>.ctrl_background.tsv B73-control ALT rate, measured fairly (see below)
#     BZea_<chr>.summary.tsv         sites removed by each filter
#   and once: pools.tsv (role, taxon group, accession of every pool)
#
#   Fair background: a control's ALT count on sites that passed a filter it helped
#   build is 0 by construction. So the 10 controls are split in two halves; each
#   half is counted on sites filtered using only the other half.
#
#   Submit: sbatch 11_filter_sites.sh
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"
need bcftools

IDX=${SLURM_ARRAY_TASK_ID:-${1:?usage: bash $0 <chrom index 1-10>}}
CHR=${CHROMS[$(( IDX - 1 ))]}
AD=${D09}/BZea_${CHR}.AD.bcf
PANEL=$(teo_panel "$CHR")
CTRL=${SAMPLE_LISTS}/b73_controls_clean.txt
EXCL=${SAMPLE_LISTS}/exclude_all.txt
for f in "$AD" "$PANEL" "${D07}/samples.txt" "$SAMPLE_META" "$CTRL" "$EXCL"; do
  [[ -s $f ]] || die "missing ${f}"
done
mkdir -p "$D11"
OUT=${D11}/BZea_${CHR}
TMP=${D11}/.tmp_${CHR}
mkdir -p "$TMP"
trap 'rm -rf "$TMP"' EXIT

# ---- pool table: SAMPLE ROLE GROUP ACCESSION (samples.txt order) -------------
tr -d '\r' < "$SAMPLE_META" \
  | awk -F'\t' -v OFS='\t' '
      FILENAME == ARGV[1] { ctrl[$1] = 1; next }
      FILENAME == ARGV[2] { ex[$1] = 1; next }
      FILENAME == ARGV[3] { if (FNR > 1) g[$1] = ($1 in g) ? g[$1] "|" $2 : $2; next }
      FNR == 1 { print "SAMPLE", "ROLE", "GROUP", "ACCESSION" }
      {
        s = $1; gt = g[s]
        if (s in ctrl)                           print s, "control", "B73", "B73"
        else if (s in ex)                        print s, "excluded", ".", "."
        else if (gt ~ /^Z[a-z]\./ && gt !~ /\|/) { acc = gt; sub(/_.*/, "", acc); print s, "donor", substr(gt, 1, 2), acc }
        else                                     print s, "other", (gt == "" ? "." : gt), "."
      }' "$CTRL" "$EXCL" /dev/stdin "${D07}/samples.txt" > "${TMP}/pools.tsv"
awk -F'\t' 'NR > 1 && ($2 == "control" || $2 == "donor") {print $1}' "${TMP}/pools.tsv" > "${TMP}/keep.txt"
if (( IDX == 1 )); then cp "${TMP}/pools.tsv" "${D11}/pools.tsv"; fi
echo "pools by role: $(awk -F'\t' 'NR > 1 {n[$2]++} END {for (r in n) printf "%s=%d ", r, n[r]}' "${TMP}/pools.tsv")"
NCTRL=$(wc -l < "$CTRL")

# ---- teosinte panel: biallelic SNPs only ----------------------------------------
bcftools view -m2 -M2 -v snps -Oz -o "${TMP}/panel.vcf.gz" "$PANEL"
bcftools index -f "${TMP}/panel.vcf.gz"

# ---- depth thresholds from the chromosome median --------------------------------
MED=$(bcftools query -f '%INFO/DP\n' "$AD" | sort -n | awk '{a[NR] = $1} END {print a[int((NR + 1) / 2)]}')
HI=$(awk -v m="$MED" -v x="$MAX_DP_MULT" 'BEGIN {printf "%d", m * x}')
LO=$(awk -v m="$MED" -v x="$MIN_DP_MULT" 'BEGIN {printf "%d", m * x}')
SITE_EXPR="INFO/DP>${HI} || INFO/DP<${LO} || INFO/OTH>${MAX_OTH_FRAC}*INFO/DP"
echo "median depth ${MED}: keep ${LO} <= DP <= ${HI}, third-base reads <= ${MAX_OTH_FRAC} x DP"

# ---- positions to remove ---------------------------------------------------------
ctrl_alt_sites() {   # $1 = control list -> CHROM POS of sites where those controls have ALT reads
  bcftools view -S "$1" -Ou "$AD" | bcftools query -i "SUM(FMT/AD[*:1])>${MAX_CTRL_ALT}" -f '%CHROM\t%POS\n'
}
ctrl_alt_sites "$CTRL" > "${TMP}/ctrl_alt.txt"
bcftools isec -c all  -n=2 -w1 "$AD" "${TMP}/panel.vcf.gz" | bcftools query -f '%CHROM\t%POS\n' | sort > "${TMP}/panel_any.txt"
bcftools isec -c none -n=2 -w1 "$AD" "${TMP}/panel.vcf.gz" | bcftools query -f '%CHROM\t%POS\n' | sort > "${TMP}/panel_same.txt"
comm -23 "${TMP}/panel_any.txt" "${TMP}/panel_same.txt" > "${TMP}/panel_mismatch.txt"

# keep sites: not in the exclude list, depth / third-base OK; then keep trusted pools
keep_sites() {   # $1 = position exclude list, $2 = sample list, $3 = output BCF
  local targets=()
  if [[ -s $1 ]]; then targets=(-T "^$1"); fi
  bcftools view "${targets[@]}" -e "$SITE_EXPR" -Ou "$AD" \
    | bcftools view -I -S "$2" -Ob -o "$3"
  bcftools index -f "$3"
}
sort -u -k1,1 -k2,2n "${TMP}/ctrl_alt.txt" "${TMP}/panel_mismatch.txt" > "${TMP}/exclude_pos.txt"
keep_sites "${TMP}/exclude_pos.txt" "${TMP}/keep.txt" "${OUT}.PASS.bcf"

# ---- split by teosinte panel ------------------------------------------------------
printf '%s\n' \
  '##INFO=<ID=PANEL_AC,Number=1,Type=Integer,Description="Teosinte panel (Schnable 2023, 227 plants): ALT allele count">' \
  '##INFO=<ID=PANEL_AN,Number=1,Type=Integer,Description="Teosinte panel: number of chromosomes called">' > "${TMP}/panel_hdr.txt"
bcftools isec -c none -n=2 -w1 -Ob -o "${TMP}/pass_panel.bcf" "${OUT}.PASS.bcf" "${TMP}/panel.vcf.gz"
bcftools index -f "${TMP}/pass_panel.bcf"          # annotate needs an indexed input, not a pipe
bcftools annotate -a "${TMP}/panel.vcf.gz" -h "${TMP}/panel_hdr.txt" \
  -c INFO/PANEL_AC:=INFO/AC,INFO/PANEL_AN:=INFO/AN -Ob -o "${OUT}.PASS.panel.bcf" "${TMP}/pass_panel.bcf"
bcftools index -f "${OUT}.PASS.panel.bcf"
bcftools isec -c none -C -w1 -Ob -o "${OUT}.PASS.nonpanel.bcf" "${OUT}.PASS.bcf" "${TMP}/panel.vcf.gz"
bcftools index -f "${OUT}.PASS.nonpanel.bcf"

# ---- reads per pool on panel / non-panel PASS sites -------------------------------
pool_sums() {   # $1 = BCF -> SAMPLE REF ALT
  paste <(bcftools query -l "$1") \
        <(bcftools query -f '[%AD{0}\t%AD{1}\t]\n' "$1" \
            | awk -F'\t' '{for (i = 1; i < NF; i++) s[i] += $i; n = NF - 1}
                          END {for (i = 1; i <= n; i += 2) print s[i] + 0 "\t" s[i + 1] + 0}')
}
join -t $'\t' <(pool_sums "${OUT}.PASS.panel.bcf" | sort -k1,1) <(pool_sums "${OUT}.PASS.nonpanel.bcf" | sort -k1,1) \
  | awk -F'\t' -v OFS='\t' -v c="$CHR" 'BEGIN {print "CHROM", "SAMPLE", "PANEL_REF", "PANEL_ALT", "NONPANEL_REF", "NONPANEL_ALT"} {print c, $0}' \
  > "${OUT}.pool_reads.tsv"

# ---- fair B73-control background (split-half) -------------------------------------
head -n $(( NCTRL / 2 )) "$CTRL" > "${TMP}/ctrlA.txt"
tail -n +$(( NCTRL / 2 + 1 )) "$CTRL" > "${TMP}/ctrlB.txt"
{
  printf "CHROM\tHALF\tPANEL_REF\tPANEL_ALT\tNONPANEL_REF\tNONPANEL_ALT\n"
  for h in A B; do
    o=$([[ $h == A ]] && echo B || echo A)
    ctrl_alt_sites "${TMP}/ctrl${o}.txt" | sort -u -k1,1 -k2,2n - "${TMP}/panel_mismatch.txt" > "${TMP}/excl_${h}.txt"
    keep_sites "${TMP}/excl_${h}.txt" "${TMP}/ctrl${h}.txt" "${TMP}/bg_${h}.bcf"
    bcftools isec -c none -n=2 -w1 -Ob -o "${TMP}/bg_${h}.panel.bcf" "${TMP}/bg_${h}.bcf" "${TMP}/panel.vcf.gz"; bcftools index -f "${TMP}/bg_${h}.panel.bcf"
    bcftools isec -c none -C   -w1 -Ob -o "${TMP}/bg_${h}.nonpanel.bcf" "${TMP}/bg_${h}.bcf" "${TMP}/panel.vcf.gz"; bcftools index -f "${TMP}/bg_${h}.nonpanel.bcf"
    p=$(pool_sums "${TMP}/bg_${h}.panel.bcf"    | awk '{r += $2; a += $3} END {print r + 0 "\t" a + 0}')
    q=$(pool_sums "${TMP}/bg_${h}.nonpanel.bcf" | awk '{r += $2; a += $3} END {print r + 0 "\t" a + 0}')
    printf "%s\t%s\t%s\t%s\n" "$CHR" "$h" "$p" "$q"
  done
} > "${OUT}.ctrl_background.tsv"

# ---- summary ----------------------------------------------------------------------
n() { bcftools index -n "$1"; }
NALL=$(n "$AD"); NPASS=$(n "${OUT}.PASS.bcf")
NCA=$(wc -l < "${TMP}/ctrl_alt.txt"); NMM=$(wc -l < "${TMP}/panel_mismatch.txt")
{
  printf "CHROM\tSITES\tCTRL_ALT\tPANEL_MISMATCH\tDEPTH_OR_THIRD_BASE\tPASS\tPASS_PANEL\tPASS_NONPANEL\tMEDIAN_DP\tDP_MIN\tDP_MAX\n"
  printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\n" "$CHR" "$NALL" "$NCA" "$NMM" \
    $(( NALL - NPASS - $(wc -l < "${TMP}/exclude_pos.txt") )) "$NPASS" \
    "$(n "${OUT}.PASS.panel.bcf")" "$(n "${OUT}.PASS.nonpanel.bcf")" "$MED" "$LO" "$HI"
} > "${OUT}.summary.tsv"
cat "${OUT}.summary.tsv"
cat "${OUT}.ctrl_background.tsv"
echo "WROTE: ${OUT}.PASS.bcf ${OUT}.PASS.panel.bcf ${OUT}.PASS.nonpanel.bcf ${OUT}.pool_reads.tsv"
