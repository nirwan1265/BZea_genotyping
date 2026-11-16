#!/bin/bash
# 26 — Summary: permutation threshold, top windows, Manhattan plots (additive and dominance).
#   Run on the login node after 26_assoc_scan.sh and 26_assoc_perm.sh:  bash 26_assoc_summary.sh
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"
need Rscript
for T in $TRAITS; do
  Rscript "${FINAL_DIR}/26_assoc_gridlmm.R" summary "$D25" "$T" "$D26"
done
