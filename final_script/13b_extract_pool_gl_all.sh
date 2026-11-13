#!/bin/bash
#SBATCH --job-name=pool_gl_all
#SBATCH --array=1-10
#SBATCH --output=logs/13b_pool_gl_all.%A_%a.out
#SBATCH --error=logs/13b_pool_gl_all.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mem=16G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# 13b — Per-pool GL files at ALL PASS sites (panel and non-panel), for the
#       family-informed HMM (14b). Same as step 13 otherwise.
#
#   Non-panel sites are only useful where a pool's BC1 family has a known donor allele
#   (step 13c); the HMM drops them elsewhere. PANEL_AC/PANEL_AN = -1 at non-panel sites.
#   Output: ${D13B}/<chr>/<SAMPLE>.tsv.gz
#     CHROM POS B73 ALT GL_0 GL_1 GL_2 PANEL_AC PANEL_AN
#
#   Submit after 11:  sbatch 13b_extract_pool_gl_all.sh
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"
need bcftools

IDX=${SLURM_ARRAY_TASK_ID:-${1:?usage: bash $0 <chrom index 1-10>}}
CHR=${CHROMS[$(( IDX - 1 ))]}
PASS=${D11}/BZea_${CHR}.PASS.bcf
PANELBCF=${D11}/BZea_${CHR}.PASS.panel.bcf
[[ -s $PANELBCF ]] || die "missing ${PANELBCF} (run 11)"
DISC=$(discovery_bcf "$CHR")
[[ -s $PASS ]] || die "missing ${PASS} (run 11)"
[[ -s $DISC ]] || die "missing discovery BCF ${DISC}"
[[ -s ${D11}/pools.tsv ]] || die "missing ${D11}/pools.tsv (run 11)"

OUTDIR=${D13B}/${CHR}
rm -rf "$OUTDIR"; mkdir -p "$OUTDIR"
TMP=${OUTDIR}/.tmp
mkdir -p "$TMP"
trap 'rm -rf "$TMP"' EXIT
ulimit -n "$(ulimit -Hn)" 2> /dev/null || true      # one open output file per pool

# sites to use: all PASS sites, with panel AC/AN where the site is a panel SNP (else -1)
bcftools query -f '%POS\t%INFO/PANEL_AC\t%INFO/PANEL_AN\n' "$PANELBCF" > "${TMP}/panel.tsv"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$PASS" \
  | awk -F'\t' -v OFS='\t' 'NR == FNR {ac[$1] = $2; an[$1] = $3; next}
                           {print $1, $2, $3, $4, (($2 in ac) ? ac[$2] : -1), (($2 in an) ? an[$2] : -1)}' \
        "${TMP}/panel.tsv" - > "${TMP}/sites.tsv"
cut -f1,2 "${TMP}/sites.tsv" > "${TMP}/targets.tsv"
awk -F'\t' 'NR > 1 && ($2 == "control" || $2 == "donor") {print $1}' "${D11}/pools.tsv" > "${TMP}/keep.txt"
bcftools view -h "$DISC" | tail -n 1 | tr '\t' '\n' | tail -n +10 | grep -Fxf "${TMP}/keep.txt" > "${TMP}/names.txt" || true
echo "${CHR}: $(wc -l < "${TMP}/sites.tsv") PASS sites (panel + non-panel) | pools in discovery BCF and kept: $(wc -l < "${TMP}/names.txt") of $(wc -l < "${TMP}/keep.txt")"

bcftools view -T "${TMP}/targets.tsv" -S "${TMP}/names.txt" -Ou "$DISC" \
  | bcftools query -f '%POS\t%REF\t%ALT[\t%GL]\n' \
  | awk -F'\t' -v OFS='\t' -v chr="$CHR" -v dir="$OUTDIR" -v sites="${TMP}/sites.tsv" -v names="${TMP}/names.txt" '
      BEGIN {
        while ((getline l < sites) > 0) { split(l, s, "\t"); b[s[2]] = s[3]; al[s[2]] = s[4]; ac[s[2]] = s[5]; an[s[2]] = s[6] }
        while ((getline l < names) > 0) { n++; f[n] = dir "/" l ".tsv" }
      }
      {
        p = $1
        if (!(p in b)) { nsk++; next }
        if ($2 == b[p] && $3 == al[p]) flip = 0
        else if ($3 == b[p] && $2 == al[p]) { flip = 1; nflip++ }
        else { nmis++; next }
        nsite++
        for (i = 1; i <= n; i++) {
          if (split($(3 + i), g, ",") != 3) continue
          if (g[1] == 0 && g[2] == 0 && g[3] == 0) continue          # no reads in this pool
          if (flip) print chr, p, b[p], al[p], g[3], g[2], g[1], ac[p], an[p] > f[i]
          else      print chr, p, b[p], al[p], g[1], g[2], g[3], ac[p], an[p] > f[i]
        }
      }
      END {
        for (i = 1; i <= n; i++) close(f[i])
        printf "sites used: %d (GL flipped to B73 orientation at %d) | allele mismatch: %d | not a PASS site: %d\n", nsite, nflip, nmis, nsk > "/dev/stderr"
      }'
gzip -f "${OUTDIR}"/*.tsv
echo "pool files: $(find "${OUTDIR}" -maxdepth 1 -name "*.tsv.gz" | wc -l)"
echo "WROTE: ${OUTDIR}/<SAMPLE>.tsv.gz"
