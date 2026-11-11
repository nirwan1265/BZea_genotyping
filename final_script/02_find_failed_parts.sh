#!/bin/bash
# =============================================================================
# 02 — Find discovery parts that did not finish and split each into sub-intervals.
#      (replaces 0.1_find_failed_runs.sh; checks outputs instead of grepping job logs)
#
#   A part counts as finished if <prefix>.done exists, or the .bcf exists and
#   ANGSD wrote "[ALL done]" to <prefix>.arg (works for the original run too).
#
#   Output: ${D03}/failed_subintervals.tsv   (CHR  START  END  PART  SUBTAG)
#
#   Run on a login/interactive node (fast):  bash 02_find_failed_parts.sh
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"

mkdir -p "$D03"
LIST=${D03}/failed_subintervals.tsv
: > "$LIST"

ntot=$(( ${#CHROMS[@]} * NPARTS ))
nfail=0
for (( i = 1; i <= ntot; i++ )); do
  read -r CHR PART <<< "$(index_to_chr_part "$i")"
  read -r START END <<< "$(part_interval "$CHR" "$PART")"
  OUT=${D01}/BZea_${CHR}.part${PART}.${START}_${END}

  if [[ -f ${OUT}.done ]] || { [[ -s ${OUT}.bcf ]] && angsd_finished "$OUT"; }; then
    continue
  fi

  nfail=$(( nfail + 1 ))
  echo "unfinished: ${CHR} part${PART} ${START}-${END}"
  len=$(( END - START + 1 ))
  sub=$(( (len + RERUN_NSUB - 1) / RERUN_NSUB ))
  for (( s = 1; s <= RERUN_NSUB; s++ )); do
    ss=$(( START + (s - 1) * sub ))
    se=$(( START + s * sub - 1 ))
    if (( se > END )); then se=$END; fi
    printf "%s\t%d\t%d\t%d\tsub%d\n" "$CHR" "$ss" "$se" "$PART" "$s" >> "$LIST"
  done
done

nsub=$(wc -l < "$LIST")
echo "Unfinished parts: ${nfail} / ${ntot}"
echo "Sub-intervals:    ${nsub}  -> ${LIST}"
if (( nsub > 0 )); then
  echo "Next: sbatch --array=1-${nsub}%60 03_angsd_rerun_failed_split.sh"
fi
