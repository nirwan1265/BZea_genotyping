#!/bin/bash
#SBATCH --job-name=filter_sites
#SBATCH --array=1-10
#SBATCH --output=logs/11_filter_sites.%A_%a.out
#SBATCH --error=logs/11_filter_sites.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=12:00:00
#SBATCH --mem=24G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# 11 — Site filtering + teosinte-panel annotation + read sums per donor group,
#      one job per chromosome.
#
#   Uses only the pools we trust: the 10 clean B73 controls and the donor pools
#   not in sample_lists/exclude_all.txt (see QC_NOTES.md).
#
#   A site is removed if:
#     - any clean B73 control has an ALT read        (repeat / paralog / reference error)
#     - total depth > MAX_DP_MULT x median, or < MIN_DP_MULT x median
#     - more than MAX_OTH_FRAC of its reads show a third base
#     - it is in the teosinte panel but with different alleles
#   Every site is also marked IN_PANEL (also a SNP in the 227-plant teosinte panel,
#   same alleles) with the panel's ALT count AC out of AN chromosomes.
#
#   Output in ${D11}, per chromosome:
#     BZea_<chr>.sites.tsv.gz     all sites: filter result, panel AC/AN, B73/ALT reads per taxon
#     BZea_<chr>.acc_counts.npz   PASS sites: B73/ALT reads per donor accession
#     BZea_<chr>.pool_class.tsv   per pool: reads on panel vs non-panel PASS sites (input to 12)
#     BZea_<chr>.summary.tsv      sites removed by each filter
#   and once: pools.tsv (role, taxon group, accession of every pool)
#
#   Submit: sbatch 11_filter_sites.sh
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"
need bcftools python3
python3 -c "import numpy, pandas" 2> /dev/null \
  || die "python3 needs numpy and pandas (conda install -c conda-forge numpy pandas)"

IDX=${SLURM_ARRAY_TASK_ID:-${1:?usage: bash $0 <chrom index 1-10>}}
CHR=${CHROMS[$(( IDX - 1 ))]}
AD=${D09}/BZea_${CHR}.AD.bcf
SITESUM=${D09}/BZea_${CHR}.sitesum.tsv.gz
PANEL=$(teo_panel "$CHR")
for f in "$AD" "$SITESUM" "$PANEL" "${D07}/samples.txt" "$SAMPLE_META" \
         "${SAMPLE_LISTS}/b73_controls_clean.txt" "${SAMPLE_LISTS}/exclude_all.txt"; do
  [[ -s $f ]] || die "missing ${f}"
done
mkdir -p "$D11"
OUT=${D11}/BZea_${CHR}
TMP=${D11}/.tmp_${CHR}
mkdir -p "$TMP"
trap 'rm -rf "$TMP"' EXIT

# pool table (identical for every chromosome; each task builds its own copy, task 1 publishes it)
python3 "${FINAL_DIR}/11_filter_sites.py" pools \
  --samples "${D07}/samples.txt" --metadata "$SAMPLE_META" \
  --controls "${SAMPLE_LISTS}/b73_controls_clean.txt" --exclude "${SAMPLE_LISTS}/exclude_all.txt" \
  --out "${TMP}/pools.tsv" --keep "${TMP}/keep.txt"
if (( IDX == 1 )); then cp "${TMP}/pools.tsv" "${D11}/pools.tsv"; fi

# teosinte panel: biallelic SNPs with AC/AN
bcftools view -m2 -M2 -v snps -Ou "$PANEL" \
  | bcftools query -f '%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\n' > "${TMP}/panel.tsv"

# AD counts of the kept pools -> filters, panel, group sums
bcftools query -s "$(paste -sd, "${TMP}/keep.txt")" \
    -f '%POS\t%REF\t%ALT\t%INFO/DP\t%INFO/OTH[\t%AD{0}\t%AD{1}]\n' "$AD" \
  | python3 "${FINAL_DIR}/11_filter_sites.py" filter \
      --pools "${TMP}/pools.tsv" --keep "${TMP}/keep.txt" --sitesum "$SITESUM" \
      --panel "${TMP}/panel.tsv" --chrom "$CHR" --out "$OUT" \
      --max-dp-mult "$MAX_DP_MULT" --min-dp-mult "$MIN_DP_MULT" \
      --max-oth "$MAX_OTH_FRAC" --max-ctrl-alt "$MAX_CTRL_ALT"

echo "WROTE: ${OUT}.sites.tsv.gz ${OUT}.acc_counts.npz ${OUT}.pool_class.tsv ${OUT}.summary.tsv"
