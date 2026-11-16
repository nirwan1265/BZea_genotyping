#!/bin/bash
#SBATCH --job-name=assoc_scan
#SBATCH --output=logs/26_assoc_scan.%j.out
#SBATCH --error=logs/26_assoc_scan.%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem=16G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# 26 — GridLMM scans (26_assoc_gridlmm.R): additive + dominance, every trait in TRAITS.
#   ~1 min per trait on 2 cores. Output: ${D26}/scan_<TRAIT>.tsv
#   Submit after 25:  sbatch 26_assoc_scan.sh
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"
need Rscript
for T in $TRAITS; do
  Rscript "${FINAL_DIR}/26_assoc_gridlmm.R" scan "$D25" "$T" "$D26" "${SLURM_CPUS_PER_TASK:-1}"
done
