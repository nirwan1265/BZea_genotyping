#!/bin/bash
#SBATCH --job-name=kmer_collect_new
#SBATCH --output=logs/21n_kmer_collect.%j.out
#SBATCH --error=logs/21n_kmer_collect.%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --mem=16G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# new_model / 21n — pick rho and write the calls of the new model (runs final_script/21_kmer_collect.py, unchanged).
#   Output: ${DNEW}/introgression/ (rho_loglik.tsv, pool_summary.tsv, tracts_all.tsv.gz, blocks_<GROUP>.pdf)
#   Run from new_model/ after 20n:  sbatch 21n_kmer_collect.sh
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/../final_script/config.sh"
need python3
DNEW=${DNEW:-${KMER_DIR}/new_model}
python3 "${FINAL_DIR}/21_kmer_collect.py" --hmm_dir "${DNEW}/hmm" --ref_dir "$D19" --out_dir "${DNEW}/introgression" \
  --window "$KMER_WINDOW" --noisy_sigma "$KMER_NOISY_SIGMA"
