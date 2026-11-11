#!/bin/bash
#SBATCH --job-name=counts_qc
#SBATCH --output=logs/10_counts_qc.%j.out
#SBATCH --error=logs/10_counts_qc.%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --mem=8G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# 10 — Per-pool QC of the read counts (wrapper for 10_counts_qc.R).
#   Output: ${D10}/sample_qc.tsv, taxon_summary.tsv, site_qc.tsv, qc_plots.png
#   Submit after 09:  sbatch --dependency=afterok:<09 job id> 10_counts_qc.sh
#   or interactively: bash 10_counts_qc.sh
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"

for c in "${CHROMS[@]}"; do
  [[ -s ${D09}/BZea_${c}.sampleqc.tsv ]] || die "missing ${D09}/BZea_${c}.sampleqc.tsv (run 09)"
done
Rscript "${FINAL_DIR}/10_counts_qc.R" "$D09" "$SAMPLE_META" "$D10"
