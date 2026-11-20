#!/bin/bash
#SBATCH --job-name=fx1_collect
#SBATCH --output=logs/fx1_5_collect.%j.out
#SBATCH --error=logs/fx1_5_collect.%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --mem=16G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# fix1 / 5 — pick rho and write the calls (final_script/21_kmer_collect.py, unchanged). Output: ${FX}/introgression/
# =============================================================================
set -euo pipefail
HERE=$(cd "${SLURM_SUBMIT_DIR:-$PWD}" && pwd)
source "${HERE}/../../final_script/config.sh"; source "${HERE}/fx1_config.sh"
FINAL_DIR=$(cd "${HERE}/../../final_script" && pwd)
python3 "${FINAL_DIR}/21_kmer_collect.py" --hmm_dir "${FXR}/hmm" --ref_dir "$D19" --out_dir "${FXR}/introgression" \
  --window "$KMER_WINDOW" --noisy_sigma "$KMER_NOISY_SIGMA"
