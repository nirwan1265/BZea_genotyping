#!/bin/bash
#SBATCH --job-name=kmer_figs
#SBATCH --output=logs/24_kmer_figures.%j.out
#SBATCH --error=logs/24_kmer_figures.%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --mem=4G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# 24 — Result figures (24_kmer_figures.py): observed vs expected BC2S3 overall, by taxon, by family,
#   by chromosome, and the per-pool donor share distribution. Needs matplotlib.
#   Output: ${D21}/figures/*.png (+ family_values.tsv), also copied to figures/results_${RUN_TAG}/ in the repo.
#   Submit after 22:  sbatch 24_kmer_figures.sh      (SPREAD/ROUND as for the other steps)
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"
need python3
python3 -c "import matplotlib" || { echo "matplotlib missing in this env" >&2; exit 1; }
python3 "${FINAL_DIR}/24_kmer_figures.py" --stats_dir "$D21" --meta "$SAMPLE_META" --out_dir "${D21}/figures"
mkdir -p "${FINAL_DIR}/figures/results_${RUN_TAG}"
cp "${D21}"/figures/*.png "${D21}/figures/family_values.tsv" "${FINAL_DIR}/figures/results_${RUN_TAG}/"
echo "copied to ${FINAL_DIR}/figures/results_${RUN_TAG}/ (git add + commit to share)"
