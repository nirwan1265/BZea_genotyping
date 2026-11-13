#!/bin/bash
#SBATCH --job-name=collect_hmm
#SBATCH --output=logs/15_collect_hmm.%j.out
#SBATCH --error=logs/15_collect_hmm.%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --mem=16G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# 15 — Collect the per-pool HMM results (wrapper for 15_collect_hmm.R).
#   Output in ${D15}: pool_summary.tsv, tracts_all.tsv.gz
#   Submit after 14:  sbatch 15_collect_hmm.sh
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"
need Rscript
Rscript "${FINAL_DIR}/15_collect_hmm.R" "$D14" "${D11}/pools.tsv" "$D15"
