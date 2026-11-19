#!/bin/bash
#SBATCH --job-name=acc_levels
#SBATCH --output=logs/20d_accession_levels.%A_%a.out
#SBATCH --error=logs/20d_accession_levels.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH --mem=24G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# new_model / 20d — DONOR level per window and ACCESSION (20d_accession_levels.py), one array task per taxon.
#   Needs 20c first (the taxon levels in ${DNEW}/levels are the prior). Writes ${DNEW}/levels_acc/<ACCESSION>.npy
#   Run from new_model/ (env with numpy, e.g. machine_learning):  sbatch --array=1-5 20d_accession_levels.sh
# =============================================================================
set -euo pipefail
HERE=$(cd "${SLURM_SUBMIT_DIR:-$PWD}" && pwd)                 # new_model/
source "${HERE}/../final_script/config.sh"
FINAL_DIR=$(cd "${HERE}/../final_script" && pwd); MAP_DIR=${FINAL_DIR}/genetic_maps
python3 -c "import numpy" 2> /dev/null || die "python3 here has no numpy"
DNEW=${DNEW:-${KMER_DIR}/new_model}
TAXA=(Zd Zh Zl Zv Zx)
T=${TAXA[$(( ${SLURM_ARRAY_TASK_ID:-1} - 1 ))]}
[[ -s ${DNEW}/levels/${T}.npy ]] || die "no taxon level ${DNEW}/levels/${T}.npy (run 20c first)"
mkdir -p "${DNEW}/levels_acc"
python3 "${HERE}/20d_accession_levels.py" --hmm_dir "$D20" --rho_file "${D21}/rho_loglik.tsv" \
  --taxon_levels "${DNEW}/levels" --meta "$SAMPLE_META" --out_dir "${DNEW}/levels_acc/${T}" --taxa "$T" --plants "$KMER_PLANTS"
cp "${DNEW}/levels_acc/${T}/"*.npy "${DNEW}/levels_acc/"
