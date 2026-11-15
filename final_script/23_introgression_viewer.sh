#!/bin/bash
#SBATCH --job-name=intro_viewer
#SBATCH --output=logs/23_introgression_viewer.%j.out
#SBATCH --error=logs/23_introgression_viewer.%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --mem=4G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# 23 — Interactive introgression viewer (23_introgression_viewer.py): one HTML file with drop-down
#   menus by taxon, family and pool; family view lines up all pools of a family; SVG/PNG download.
#   Output: ${D21}/introgression_viewer.html  (also copied to figures/ in the repo, to commit)
#   ROUND (default 0) picks the round: sbatch --export=ALL,ROUND=1 23_introgression_viewer.sh
#   Submit after 21:  sbatch 23_introgression_viewer.sh   (or on the login node: bash 23_introgression_viewer.sh)
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"
need python3
python3 "${FINAL_DIR}/23_introgression_viewer.py" --tracts "${D21}/tracts_all.tsv.gz" \
  --pool_summary "${D21}/pool_summary.tsv" --meta "$SAMPLE_META" --out "${D21}/introgression_viewer.html" \
  --title "BZea BC2S3 introgressions (k-mer HMM, ${RUN_TAG})"
cp "${D21}/introgression_viewer.html" "${FINAL_DIR}/figures/introgression_viewer_${RUN_TAG}.html"
echo "copied to ${FINAL_DIR}/figures/introgression_viewer_${RUN_TAG}.html (git add + commit to share)"
