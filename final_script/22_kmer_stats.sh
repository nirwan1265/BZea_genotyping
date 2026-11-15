#!/bin/bash
#SBATCH --job-name=kmer_stats
#SBATCH --output=logs/22_kmer_stats.%j.out
#SBATCH --error=logs/22_kmer_stats.%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --mem=8G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# 22 — Observed vs expected BC2S3 proportions with a family bootstrap (22_kmer_stats.py).
#   Output in ${D21}: stats_overall.tsv, stats_by_chrom.tsv
#   Submit after 21:  sbatch 22_kmer_stats.sh   (or run on the login node: bash 22_kmer_stats.sh)
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"
need python3
python3 "${FINAL_DIR}/22_kmer_stats.py" --pool_summary "${D21}/pool_summary.tsv" --tracts "${D21}/tracts_all.tsv.gz" \
  --meta "$SAMPLE_META" --out_dir "$D21"
