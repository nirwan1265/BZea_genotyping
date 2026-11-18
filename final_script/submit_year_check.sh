#!/bin/bash
# =============================================================================
# Check: the k-mer ancestry GWAS (25-27) and the joint QTL model (29.2) on each year's BLUEs separately,
# to see whether the QTL of the two-year analysis are found in 2023 and in 2025 on their own.
# Run on the login node, from final_script/, with the GWAS conda env active (jobs inherit it):
#   bash submit_year_check.sh                       # DTS, both years
#   TRAITS="DTS PH" bash submit_year_check.sh       # other traits (SL, BL, BW exist only in 2023)
# Per year, everything goes to its own folders (taxon as the fixed effect, 1,000 permutations):
#   ${D21}/assoc_y<YEAR>/                      inputs
#   ${D21}/assoc_gridlmm_spats_y<YEAR>/        scan, permutations, threshold, peaks, genes
#   ${D29}/noperm_y<YEAR>/                     joint QTL model
# The two-year results are not touched.
# =============================================================================
set -euo pipefail
export TRAITS=${TRAITS:-DTS}
source config.sh
mkdir -p logs
Rscript -e 'library(GridLMM); library(qtl2)' > /dev/null 2>&1 || { echo "GridLMM/qtl2 not found: activate the GWAS conda env first"; exit 1; }
NTR=$(wc -w <<< "$TRAITS")
for Y in 2023 2025; do
  P=${FINAL_DIR}/../data/phenotypes/BZea_spats_BLUEs_${Y}.csv
  [[ -s "$P" ]] || { echo "missing $P"; exit 1; }
  E="--export=ALL,ASSOC_TAG=y$Y,PHENO=$P,TRAITS=$TRAITS,ASSOC_FIXED_COL=TAXON,TAXON_COL=TAXON"
  J25=$(sbatch --parsable $E 25_assoc_inputs.sh)
  JS=$(sbatch --parsable $E --dependency=afterok:$J25 26_assoc_scan.sh)
  JP=$(sbatch --parsable $E --dependency=afterok:$J25 --array=1-100 26_assoc_perm.sh)
  JF=$(sbatch --parsable $E --dependency=afterok:$JS:$JP 26_assoc_finish.sh)
  JQ=$(sbatch --parsable $E --dependency=afterok:$JF --array=1-$NTR 29.2_joint_qtl_noperm.sh)
  echo "$Y ($TRAITS): inputs $J25 | scan $JS | permutations $JP | finish $JF | joint QTL $JQ"
done
echo "results: ${D21}/assoc_gridlmm_spats_y{2023,2025} and ${D29}/noperm_y{2023,2025}"
