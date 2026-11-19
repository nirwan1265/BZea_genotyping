#!/bin/bash
#SBATCH --job-name=pool_levels
#SBATCH --output=logs/20c_pool_levels.%A_%a.out
#SBATCH --error=logs/20c_pool_levels.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mem=24G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# new_model / 20c — DONOR level per window learned from all pools of a taxon (20c_pool_levels.py).
#   One array task per taxon. Reads the CURRENT run (config defaults: ROUND=0, SPREAD=1 -> hmm_spread,
#   introgression_spread); writes ${DNEW}/levels/<TAXON>.npy. Nothing in the current run is changed.
#   Run from new_model/ with any env that has numpy (e.g. machine_learning):
#     mkdir -p logs && sbatch --array=1-5 20c_pool_levels.sh
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/../final_script/config.sh"
python3 -c "import numpy" 2> /dev/null || die "python3 here has no numpy"
DNEW=${DNEW:-${KMER_DIR}/new_model}
TAXA=(Zd Zh Zl Zv Zx)
T=${TAXA[$(( ${SLURM_ARRAY_TASK_ID:-1} - 1 ))]}
ASMS=$(for a in "${DONOR_ASSEMBLIES[@]}"; do IFS=: read -r n t _ <<< "$a"; echo "${n}:${t}"; done | paste -sd,)
mkdir -p "${DNEW}/levels"
echo "taxon ${T}: current run ${D20} -> ${DNEW}/levels"
python3 "${SLURM_SUBMIT_DIR:-$PWD}/20c_pool_levels.py" --hmm_dir "$D20" --rho_file "${D21}/rho_loglik.tsv" \
  --ref_dir "$D19" --assemblies "$ASMS" --out_dir "${DNEW}/levels" --taxa "$T" --plants "$KMER_PLANTS"
