#!/bin/bash
# 27 — Peaks and the B73 v5 genes inside them (27_assoc_genes.py), for every trait and scope, at the
#   step-26 permutation threshold. Also lists CCT-domain genes (ZmCCT10 etc.) and the scan near them.
#   Run on the login node after 26_assoc_summary.sh:  bash 27_assoc_genes.sh
#   Output: ${D26}/genes/peaks_<TAG>.tsv, peak_genes_<TAG>.tsv, find_<TAG>.tsv
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"
need python3
for T in $TRAITS; do
  tail -n +2 "${D26}/summary_${T}.tsv" | while IFS=$'\t' read -r trait scope plots windows nperm thr rest; do
    [[ "$thr" == "NA" || -z "$thr" ]] && { echo "$T $scope: no permutation threshold yet, skipped"; continue; }
    if [[ "$scope" == "ALL" ]]; then f="${D26}/scan_${T}.tsv"; else f="${D26}/scan_${T}_${scope}.tsv"; fi
    echo "== $T $scope (threshold $thr)"
    python3 "${FINAL_DIR}/27_assoc_genes.py" --scan "$f" --threshold "$thr" --gff "$GFF_GLOB" \
      --out_dir "${D26}/genes" --tag "${T}_${scope}" --find 'CCT'
  done
done
