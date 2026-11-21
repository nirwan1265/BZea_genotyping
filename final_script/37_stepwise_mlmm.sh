#!/bin/bash
#SBATCH --job-name=stepwise
#SBATCH --output=logs/37_stepwise.%A_%a.out
#SBATCH --error=logs/37_stepwise.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --mem=16G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# 37 — Stepwise multi-locus mixed model (37_stepwise_mlmm.R), one trait per array task.
#   Uses the step-25 inputs (${D25}) and the step-26 permutation threshold (${D26}/summary_<TRAIT>.tsv);
#   output: ${D21}/stepwise/stepwise_<TRAIT>_qtl.tsv and _steps.tsv. A few minutes per forward step on 4 cores.
#   Traits: SW_TRAITS (default below). Jobs use the env active at submission (GridLMM, data.table).
#   Submit (new model):
#     cd final_script && export MODEL=fix1
#     sbatch --array=1-12%100 37_stepwise_mlmm.sh
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"
SW_TRAITS=${SW_TRAITS:-"DTS DTA ASI NBR SPAD2 PH EH EN Prolif SL BL BW"}
T=$(echo $SW_TRAITS | cut -d' ' -f${SLURM_ARRAY_TASK_ID})
echo "trait $T | inputs $D25 | threshold from $D26/summary_${T}.tsv"
Rscript "${FINAL_DIR}/37_stepwise_mlmm.R" "$D25" "$T" "$D26/summary_${T}.tsv" "${D21}/stepwise" "${SLURM_CPUS_PER_TASK:-1}"
