#!/bin/bash
#SBATCH --job-name=assoc_finish
#SBATCH --output=logs/26_assoc_finish.%j.out
#SBATCH --error=logs/26_assoc_finish.%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --mem=16G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# 26 + 27 as one job: after 26_assoc_scan.sh and the 26_assoc_perm.sh array, for every trait in $TRAITS:
#   summary (permutation threshold, top windows, Manhattan plots) and peaks + B73 v5 genes (step 27).
#   sbatch --dependency=afterok:<scan job>:<perm job> 26_assoc_finish.sh
# =============================================================================
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
bash 26_assoc_summary.sh
bash 27_assoc_genes.sh
