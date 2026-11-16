#!/bin/bash
# =============================================================================
# Submit the whole association analysis on the cluster in one go (run on the login node, from final_script/,
# with the GWAS conda env active; jobs inherit it):
#   25  inputs (two-year BLUEs, SUBTAXON from data/metadata_BZea.csv)
#   for the fixed effect TAXON and SUBTAXON (mexicana split by race), each in its own folder:
#     26  scan (all traits) + 1,000 permutations (100 tasks x 10) -> 26/27 finish (threshold, plots, peaks, genes)
#     29.2 joint QTL model on those peaks (one task per trait; per-taxon or per-race effects)
#   bash submit_assoc_all.sh            # both
#   FIXEDS="TAXON" bash submit_assoc_all.sh
# =============================================================================
set -euo pipefail
source config.sh
mkdir -p logs
FIXEDS=${FIXEDS:-"TAXON SUBTAXON"}
NTR=$(wc -w <<< "$TRAITS")
J25=$(sbatch --parsable 25_assoc_inputs.sh); echo "25 inputs: $J25"
for FX in $FIXEDS; do
  E="--export=ALL,ASSOC_FIXED_COL=$FX,TAXON_COL=$FX"
  JS=$(sbatch --parsable $E --dependency=afterok:$J25 26_assoc_scan.sh)
  JP=$(sbatch --parsable $E --dependency=afterok:$J25 --array=1-100 26_assoc_perm.sh)
  JF=$(sbatch --parsable $E --dependency=afterok:$JS:$JP 26_assoc_finish.sh)
  JQ=$(sbatch --parsable $E --dependency=afterok:$JF --array=1-$NTR 29.2_joint_qtl_noperm.sh)
  echo "$FX: scan $JS | permutations $JP | finish $JF | joint QTL $JQ"
done
echo "squeue -u $USER to follow; results in ${D21}/assoc_gridlmm_spats{,_subtaxon} and ${D29}/noperm{,_subtaxon}"
