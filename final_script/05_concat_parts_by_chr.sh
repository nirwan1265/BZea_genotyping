#!/bin/bash
#SBATCH --job-name=concat_chr
#SBATCH --array=1-10
#SBATCH --output=logs/05_concat_chr.%A_%a.out
#SBATCH --error=logs/05_concat_chr.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=10-00:00:00
#SBATCH --mem=20G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# 05 — Concatenate the 20 discovery parts of each chromosome, in coordinate order.
#      (replaces 0.4_combine_parts_chrom.sh; one array job per chromosome)
#
#   Part file names are computed from the .fai, not globbed, so a missing or
#   duplicated part stops the job instead of leaving a silent gap.
#
#   Output: ${D05}/BZea_<chr>.full.bcf (+ .csi)
#   Submit: sbatch 05_concat_parts_by_chr.sh
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"
need bcftools

IDX=${SLURM_ARRAY_TASK_ID:-${1:?usage: bash $0 <chrom index 1-10>}}
CHR=${CHROMS[$(( IDX - 1 ))]}
mkdir -p "$D05"
OUT=${D05}/BZea_${CHR}.full.bcf
LISTF=${D05}/BZea_${CHR}.parts.list
: > "$LISTF"

while read -r PART START END; do
  P=${D01}/BZea_${CHR}.part${PART}.${START}_${END}
  if [[ -f ${P}.done ]] || { [[ -s ${P}.bcf ]] && angsd_finished "$P"; }; then
    echo "${P}.bcf" >> "$LISTF"
  else
    die "part not finished: ${P} (run 02-04)"
  fi
done < <(parts_of_chr "$CHR")

echo "==== ${CHR}: $(wc -l < "$LISTF") parts"
bcftools concat -f "$LISTF" -Ob -o "$OUT"
bcftools index -f "$OUT"
echo "SNPs in ${CHR}: $(bcftools index -n "$OUT")"
echo "WROTE: ${OUT}"
