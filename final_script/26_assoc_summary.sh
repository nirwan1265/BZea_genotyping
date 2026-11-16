#!/bin/bash
# 26 — Summary: permutation threshold, top windows, Manhattan plots (additive and dominance).
#   Run on the login node after 26_assoc_scan.sh and 26_assoc_perm.sh:  bash 26_assoc_summary.sh
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"
need Rscript
for T in $TRAITS; do
  GROUP_COL=$ASSOC_GROUP_COL FIXED_COL=$ASSOC_FIXED_COL MIN_N=$ASSOC_MIN_N MIN_CARRIERS=$ASSOC_GROUP_MIN_CARRIERS \
    Rscript "${FINAL_DIR}/26_assoc_gridlmm.R" summary "$D25" "$T" "$D26"
done
