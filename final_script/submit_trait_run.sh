#!/bin/bash
# =============================================================================
# Run the k-mer ancestry GWAS (25-27) and the joint QTL model (29.2) for traits in their own phenotype file,
# into their own folders (the main two-year results are not touched). Example: predicted leaf N, 2025 block 2.
# Run on the login node, from final_script/, with the GWAS conda env active (jobs inherit it):
#   conda activate /usr/local/usrapps/maize/ntanduk/bzea_gwas
#   bash submit_trait_run.sh                                                   # PredN (defaults below)
#   PHENO_FILE=../data/phenotypes/X.csv TRAITS="A B" TAG=x bash submit_trait_run.sh
# Outputs (taxon as the fixed effect, 1,000 permutations):
#   ${D21}/assoc_<TAG>/                      inputs
#   ${D21}/assoc_gridlmm_spats_<TAG>/        scan, permutations, threshold, peaks, genes
#   ${D29}/noperm_<TAG>/                     joint QTL model
# =============================================================================
set -euo pipefail
export TRAITS=${TRAITS:-PredN}
TAG=${TAG:-N2025}
P=$(readlink -f "${PHENO_FILE:-../data/phenotypes/BZea_spats_BLUEs_N_2025.csv}")
source config.sh
mkdir -p logs
[[ -s "$P" ]] || { echo "missing $P"; exit 1; }
Rscript -e 'library(GridLMM); library(qtl2)' > /dev/null 2>&1 || { echo "GridLMM/qtl2 not found: activate the GWAS conda env first"; exit 1; }
NTR=$(wc -w <<< "$TRAITS")
E="--export=ALL,ASSOC_TAG=$TAG,PHENO=$P,TRAITS=$TRAITS,ASSOC_FIXED_COL=TAXON,TAXON_COL=TAXON"
J25=$(sbatch --parsable $E 25_assoc_inputs.sh)
JS=$(sbatch --parsable $E --dependency=afterok:$J25 26_assoc_scan.sh)
JP=$(sbatch --parsable $E --dependency=afterok:$J25 --array=1-100 26_assoc_perm.sh)
JF=$(sbatch --parsable $E --dependency=afterok:$JS:$JP 26_assoc_finish.sh)
JQ=$(sbatch --parsable $E --dependency=afterok:$JF --array=1-$NTR 29.2_joint_qtl_noperm.sh)
echo "$TAG ($TRAITS): inputs $J25 | scan $JS | permutations $JP | finish $JF | joint QTL $JQ"
echo "results: ${D21}/assoc_gridlmm_spats_$TAG and ${D29}/noperm_$TAG"
