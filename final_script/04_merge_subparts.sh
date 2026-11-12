#!/bin/bash
#SBATCH --job-name=merge_subparts
#SBATCH --output=logs/04_merge_subparts.%j.out
#SBATCH --error=logs/04_merge_subparts.%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=10-00:00:00
#SBATCH --mem=20G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# 04 — Concatenate the sub-interval BCFs from 03 back into whole-part BCFs and
#      put them next to the step-01 parts, so step 05 sees a complete set.
#      (replaces 0.3_combine_sub_parts.sh; handles any number of sub-intervals)
#
#   Output: ${D01}/BZea_<chr>.part<k>.<start>_<end>.bcf (+ .csi, .done)
#   Submit: sbatch 04_merge_subparts.sh
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"
need bcftools

LIST=${D03}/failed_subintervals.tsv
[[ -s $LIST ]] || { echo "No ${LIST}; nothing to merge."; exit 0; }

# unique CHR/PART pairs, in file order
cut -f1,4 "$LIST" | awk '!seen[$0]++' | while read -r CHR PART; do
  read -r PSTART PEND <<< "$(part_interval "$CHR" "$PART")"
  OUT=${D01}/BZea_${CHR}.part${PART}.${PSTART}_${PEND}

  # sub-interval files for this part, sorted by start
  mapfile -t subs < <(
    awk -v c="$CHR" -v p="$PART" -F'\t' '$1==c && $4==p {print $2"\t"$3"\t"$5}' "$LIST" \
      | sort -k1,1n \
      | while IFS=$'\t' read -r s e tag; do echo "${D03}/BZea_${CHR}.part${PART}.${tag}.${s}_${e}"; done
  )

  for s in "${subs[@]}"; do
    [[ -f ${s}.done ]] || die "sub-interval not finished: ${s} (rerun 03)"
  done

  echo "==== ${CHR} part${PART}: ${#subs[@]} sub-intervals -> ${OUT}.bcf"
  printf '%s.bcf\n' "${subs[@]}" > "${OUT}.sublist"
  bcftools concat -f "${OUT}.sublist" -Ob -o "${OUT}.bcf"
  bcftools index -f "${OUT}.bcf"
  rm -f "${OUT}.sublist"
  touch "${OUT}.done"
done

echo "DONE."
