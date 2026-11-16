#!/bin/bash
#SBATCH --job-name=assoc_perm
#SBATCH --output=logs/26_assoc_perm.%A_%a.out
#SBATCH --error=logs/26_assoc_perm.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=04:00:00
#SBATCH --mem=8G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# 26 — Permutation threshold (26_assoc_gridlmm.R perm): phenotypes shuffled within taxon,
#   genome-wide max -log10 p per permutation, for all plots and every group. ~1 min per permutation
#   (all scopes) on 2 cores. One task = ASSOC_PERM_PER_TASK permutations of every trait (seed = task id).
#   Submit after 25:  sbatch --array=1-100 26_assoc_perm.sh       (100 x 10 = 1,000 permutations, ~20 min per task)
#   Then:             bash 26_assoc_summary.sh   (threshold, top windows, Manhattan plots)
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"
need Rscript
for T in $TRAITS; do
  GROUP_COL=$ASSOC_GROUP_COL MIN_N=$ASSOC_MIN_N MIN_CARRIERS=$ASSOC_GROUP_MIN_CARRIERS \
    Rscript "${FINAL_DIR}/26_assoc_gridlmm.R" perm "$D25" "$T" "$D26" "$ASSOC_PERM_PER_TASK" "$SLURM_ARRAY_TASK_ID" "${SLURM_CPUS_PER_TASK:-1}"
done
