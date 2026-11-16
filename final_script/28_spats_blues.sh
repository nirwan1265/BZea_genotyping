#!/bin/bash
#SBATCH --job-name=spats_blues
#SBATCH --output=logs/28_spats_blues.%j.out
#SBATCH --error=logs/28_spats_blues.%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --mem=8G
# =============================================================================
# 28 — SpATS BLUEs for 2023 (CLY23_D4) and 2025 (CLY25_B5), all traits, + combined across years.
#   ~35 min. Output in data/phenotypes/ (BZea_spats_BLUEs_{2023,2025,combined}.csv, summary, PDFs).
#   Needs R with SpATS, dplyr, tidyr. Locally: cd <repo> && Rscript final_script/28_spats_blues.R
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"
cd "${FINAL_DIR}/.."
Rscript final_script/28_spats_blues.R
