#!/bin/bash
#SBATCH --job-name=kmer_collect
#SBATCH --output=logs/21_kmer_collect.%j.out
#SBATCH --error=logs/21_kmer_collect.%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --mem=16G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# 21 — Pick rho (highest total log-likelihood), write calls and figures (21_kmer_collect.py).
#   Output in ${D21}: rho_loglik.tsv, pool_summary.tsv, tracts_all.tsv.gz, blocks_<GROUP>.pdf
#   ROUND (default 0) picks the round: sbatch --export=ALL,ROUND=1 21_kmer_collect.sh
#   Submit after 20:  sbatch 21_kmer_collect.sh
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"
need python3
python3 "${FINAL_DIR}/21_kmer_collect.py" --hmm_dir "$D20" --ref_dir "$D19" --out_dir "$D21" \
  --window "$KMER_WINDOW" --noisy_sigma "$KMER_NOISY_SIGMA"
