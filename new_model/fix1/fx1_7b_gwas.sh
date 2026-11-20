#!/bin/bash
# =============================================================================
# fix1 / 7b — k-mer GWAS + 1,000 permutations + joint QTL on the new calls (MODEL=fix1): taxon and mexicana races
#   (submit_assoc_all.sh), each year (submit_year_check.sh), predicted leaf N 2025 (submit_trait_run.sh).
#   Outputs next to the calls: ${KMER_DIR}/new_model/fix1/with_levels/introgression/. Jobs use the env active when
#   you submit (the GWAS env: R with GridLMM, qtl2). Run on the login node:  cd new_model/fix1 && bash fx1_7b_gwas.sh
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")/../../final_script"
export MODEL=fix1
source config.sh
mkdir -p logs
[[ -s ${D21}/tracts_all.tsv.gz ]] || die "no calls in ${D21} (fx1_submit_levels.sh first)"
bash submit_assoc_all.sh
bash submit_year_check.sh
bash submit_trait_run.sh
echo "outputs in ${D21}"
