#!/bin/bash
#SBATCH --job-name=fx1_choose
#SBATCH --output=logs/fx1_3_choose.%A_%a.out
#SBATCH --error=logs/fx1_3_choose.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=06:00:00
#SBATCH --mem=24G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# fix1 / 3 — per taxon (array 1-5: Zd Zh Zl Zv Zx): which windows switch to the informative k-mers, and their
#   DONOR level and spread (fx1_3_choose.py). Reads the current run (hmm_spread, introgression_spread).
#   Output: ${FX}/choice/<T>.npz, choice_<T>.tsv.gz, summary_<T>.tsv
# =============================================================================
set -euo pipefail
HERE=$(cd "${SLURM_SUBMIT_DIR:-$PWD}" && pwd)
source "${HERE}/../../final_script/config.sh"; source "${HERE}/fx1_config.sh"
python -c "import numpy, scipy, pandas" 2> /dev/null || die "python here has no numpy/scipy/pandas: activate your python env before submitting"
T=${FX_TAXA[$(( ${SLURM_ARRAY_TASK_ID:-1} - 1 ))]}
mkdir -p "${FX}/choice"
python "${HERE}/fx1_3_choose.py" --taxon "$T" --fx_dir "$FX" --d16 "$D16" --d19 "$D19" --pools_tsv "${D11}/pools.tsv" \
  --hmm_dir "$D20" --rho_file "${D21}/rho_loglik.tsv" --tracts "${D21}/tracts_all.tsv.gz" --out_dir "${FX}/choice" \
  --min_new "$FX_MIN_NEW" --min_pools "$FX_MIN_POOLS" --gain "$FX_GAIN" --noisy_sigma "$KMER_NOISY_SIGMA"
