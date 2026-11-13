#!/bin/bash
#SBATCH --job-name=call_introgression
#SBATCH --output=logs/14_call_introgression.%j.out
#SBATCH --error=logs/14_call_introgression.%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mem=32G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# 14 — Call introgressions for every pool (wrapper for 14_call_introgression.R).
#
#   Each window of each pool gets B73 / SEG (segregating) / DONOR, with posterior
#   probabilities and an expected donor dosage (0..1). See the R script header.
#
#   Output in ${D14}: calls.tsv.gz, dosage_matrix.tsv.gz, segments.tsv, pool_summary.tsv
#   Submit after 13:  sbatch 14_call_introgression.sh
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"
need Rscript

[[ -s ${D11}/pools.tsv ]] || die "missing ${D11}/pools.tsv (run 11)"
for c in "${CHROMS[@]}"; do
  [[ -s ${D13}/BZea_${c}.win_counts.tsv.gz ]] || die "missing ${D13}/BZea_${c}.win_counts.tsv.gz (run 13)"
done
Rscript "${FINAL_DIR}/14_call_introgression.R" "$D13" "${D11}/pools.tsv" "${REF}.fai" "$D14" \
  "$WIN_SIZE" "$SWITCH_PROB" "$(IFS=,; echo "${CHROMS[*]}")"
