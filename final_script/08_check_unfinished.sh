#!/bin/bash
# =============================================================================
# 08_check_unfinished.sh — list step-08 array indices without a .done marker and
#                          print the bsub command that resubmits only those.
#   Run on a login node:  bash 08_check_unfinished.sh
# =============================================================================
set -euo pipefail
source "${LS_SUBCWD:-$PWD}/config.sh"

todo=()
ntot=$(( ${#CHROMS[@]} * NPARTS ))
for (( i = 1; i <= ntot; i++ )); do
  read -r CHR PART <<< "$(index_to_chr_part "$i")"
  read -r START END <<< "$(part_interval "$CHR" "$PART")"
  [[ -f ${D08}/BZea_${CHR}.part${PART}.${START}_${END}.done ]] || todo+=("$i")
done

echo "finished: $(( ntot - ${#todo[@]} )) / ${ntot}"
if (( ${#todo[@]} )); then
  list=$(IFS=,; echo "${todo[*]}")
  echo "unfinished indices: ${list}"
  echo "resubmit: bsub -J \"angsd_counts[${list}]\" < 08_angsd_counts_by_part.sh"
fi
