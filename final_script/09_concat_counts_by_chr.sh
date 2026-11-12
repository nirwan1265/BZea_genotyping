#!/bin/bash
#SBATCH --job-name=concat_counts
#SBATCH --array=1-10
#SBATCH --output=logs/09_concat_counts.%A_%a.out
#SBATCH --error=logs/09_concat_counts.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mem=8G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# 09 — Join the 20 count parts of each chromosome (one job per chromosome).
#
#   Output in ${D09}:
#     BZea_<chr>.AD.bcf (+ .csi)   FINAL read counts: REF = B73 base, FORMAT/AD = B73,ALT reads
#     BZea_<chr>.sitesum.tsv.gz    per-site totals (DP, OTH, NCOV, NALT, B73-control reads)
#     BZea_<chr>.sampleqc.tsv      per-pool totals for this chromosome
#     BZea_<chr>.siteqc.tsv        one-line site summary for step 10
#
#   Quick look:  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' BZea_chr1.AD.bcf | head
#   Submit after 08 finishes:  sbatch --dependency=afterok:<08 job id> 09_concat_counts_by_chr.sh
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"
need bcftools

IDX=${SLURM_ARRAY_TASK_ID:-${1:?usage: bash $0 <chrom index 1-10>}}
CHR=${CHROMS[$(( IDX - 1 ))]}
mkdir -p "$D09"
OUT=${D09}/BZea_${CHR}

prefixes=()
while read -r PART START END; do
  P=${D08}/BZea_${CHR}.part${PART}.${START}_${END}
  [[ -f ${P}.done ]] || die "part not finished: ${P} (rerun 08 for it)"
  prefixes+=("$P")
done < <(parts_of_chr "$CHR")
echo "==== ${CHR}: ${#prefixes[@]} parts"

# BCF: identical headers (written with --no-version) -> fast naive concat
printf '%s.AD.bcf\n' "${prefixes[@]}" > "${OUT}.parts.list"
bcftools concat --naive -f "${OUT}.parts.list" -Ob -o "${OUT}.AD.bcf"
bcftools index -f "${OUT}.AD.bcf"
rm -f "${OUT}.parts.list"

# per-site summary: one header, then all bodies in order
{
  gzip -dc "${prefixes[0]}.sitesum.tsv.gz" | head -n 1
  for P in "${prefixes[@]}"; do gzip -dc "${P}.sitesum.tsv.gz" | tail -n +2; done
} | gzip -c > "${OUT}.sitesum.tsv.gz"

# per-pool totals summed over parts
qcfiles=()
for P in "${prefixes[@]}"; do qcfiles+=("${P}.sampleqc.tsv"); done
awk -F'\t' -v OFS='\t' '
  FNR == 1 { next }
  { if (!($1 in c)) ord[++n] = $1; c[$1] += $2; r[$1] += $3; a[$1] += $4; o[$1] += $5 }
  END {
    print "SAMPLE", "SITES_COVERED", "REF_READS", "ALT_READS", "OTHER_READS"
    for (i = 1; i <= n; i++) { s = ord[i]; print s, c[s], r[s], a[s], o[s] }
  }' "${qcfiles[@]}" > "${OUT}.sampleqc.tsv"

# one-line site summary
gzip -dc "${OUT}.sitesum.tsv.gz" | awk -F'\t' -v OFS='\t' -v chr="$CHR" '
  NR == 1 { next }
  { n++; dp += $5; oth += $6; if ($12 > 0) ctrl_alt++; if ($10 == 0) no_alt++
    if ($5 > 0 && $6 / $5 > 0.10) oth10++ }
  END {
    print "CHROM", "N_SITES", "MEAN_DP_ALL_POOLS", "SITES_B73CTRL_ALT", "SITES_NO_ALT_READ", "SITES_OTH_GT10PCT", "OTH_FRAC"
    print chr, n + 0, (n ? dp / n : 0), ctrl_alt + 0, no_alt + 0, oth10 + 0, (dp ? oth / dp : 0)
  }' > "${OUT}.siteqc.tsv"

N_BCF=$(bcftools index -n "${OUT}.AD.bcf")
N_SUM=$(( $(gzip -dc "${OUT}.sitesum.tsv.gz" | wc -l) - 1 ))
(( N_BCF == N_SUM )) || die "BCF (${N_BCF}) and site summary (${N_SUM}) disagree"
cat "${OUT}.siteqc.tsv"
echo "WROTE: ${OUT}.AD.bcf"
