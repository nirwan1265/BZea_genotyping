#!/bin/bash
#BSUB -J counts_qc
#BSUB -n 1
#BSUB -W 2:00
#BSUB -q sara
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=8GB]"
#BSUB -o logs/10_counts_qc.%J.out
#BSUB -e logs/10_counts_qc.%J.err
# =============================================================================
# 10 — Per-pool QC of the read counts (wrapper for 10_counts_qc.R).
#   Output: ${D10}/sample_qc.tsv, taxon_summary.tsv, site_qc.tsv, qc_plots.png
#   Submit after 09:  bsub -w "done(concat_counts)" < 10_counts_qc.sh
#   or interactively: bash 10_counts_qc.sh
# =============================================================================
set -euo pipefail
source "${LS_SUBCWD:-$PWD}/config.sh"

for c in "${CHROMS[@]}"; do
  [[ -s ${D09}/BZea_${c}.sampleqc.tsv ]] || die "missing ${D09}/BZea_${c}.sampleqc.tsv (run 09)"
done
Rscript "${FINAL_DIR}/10_counts_qc.R" "$D09" "$SAMPLE_META" "$D10"
