#!/bin/bash
#BSUB -J "concat_chr[1-10]%10"
#BSUB -n 1
#BSUB -W 14400
#BSUB -q sara
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=20GB]"
#BSUB -o logs/05_concat_chr.%J.%I.out
#BSUB -e logs/05_concat_chr.%J.%I.err
# =============================================================================
# 05 — Concatenate the 20 discovery parts of each chromosome, in coordinate order.
#      (replaces 0.4_combine_parts_chrom.sh; one array job per chromosome)
#
#   Part file names are computed from the .fai, not globbed, so a missing or
#   duplicated part stops the job instead of leaving a silent gap.
#
#   Output: ${D05}/BZea_<chr>.full.bcf (+ .csi)
#   Submit: bsub < 05_concat_parts_by_chr.sh
# =============================================================================
set -euo pipefail
source "${LS_SUBCWD:-$PWD}/config.sh"

IDX=${LSB_JOBINDEX:-${1:?usage: bash $0 <chrom index 1-10>}}
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
