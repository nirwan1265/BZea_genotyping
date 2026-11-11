#!/bin/bash
#BSUB -J "angsd_rerun[1-60]%60"
#BSUB -n 4
#BSUB -W 5000
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=48GB]"
#BSUB -o logs/03_angsd_rerun.%J.%I.out
#BSUB -e logs/03_angsd_rerun.%J.%I.err
# =============================================================================
# 03 — Rerun unfinished discovery parts as smaller sub-intervals with more memory.
#      (replaces 0.2_run_failed_parts.sh)
#
#   Reads ${D03}/failed_subintervals.tsv from step 02. Set the array size to the
#   number of lines in that file (02 prints the exact bsub command):
#       bsub -J "angsd_rerun[1-N]%60" < 03_angsd_rerun_failed_split.sh
#
#   Uses RERUN_* settings from config.sh. The original rerun used SNP_pval 1e-10,
#   minInd 20, setMinDepthInd 4 (stricter than step 01), so the killed parts have
#   a stricter SNP set than the rest of the genome. See README.
#
#   Output: ${D03}/BZea_<chr>.part<k>.sub<s>.<start>_<end>.bcf (+ .csi, .done)
# =============================================================================
set -euo pipefail
source "${LS_SUBCWD:-$PWD}/config.sh"

LIST=${D03}/failed_subintervals.tsv
[[ -s $LIST ]] || die "missing or empty ${LIST} (run 02 first)"

IDX=${LSB_JOBINDEX:-${1:?usage: bash $0 <line number in failed_subintervals.tsv>}}
line=$(sed -n "${IDX}p" "$LIST")
[[ -n $line ]] || die "no line ${IDX} in ${LIST}"
read -r CHR START END PART SUB <<< "$line"

OUT=${D03}/BZea_${CHR}.part${PART}.${SUB}.${START}_${END}
echo "CHR=${CHR} PART=${PART} ${SUB} INTERVAL=${CHR}:${START}-${END}"
echo "OUT=${OUT}"

[[ -f ${OUT}.done ]] && { echo "Already done; skipping."; exit 0; }

angsd \
  -bam "$BAMLIST" \
  -ref "$REF" \
  -r "${CHR}:${START}-${END}" \
  -GL 1 \
  -doMajorMinor 1 \
  -doMaf 1 \
  -doBcf 1 \
  -SNP_pval "$RERUN_SNP_PVAL" \
  "${READ_FILTERS[@]}" \
  -minInd "$RERUN_MININD" \
  -skipTriallelic 1 \
  -setMinDepthInd "$RERUN_SETMINDEPTHIND" \
  -P "$THREADS" \
  -out "$OUT"

angsd_finished "$OUT" || die "ANGSD did not finish (no [ALL done] in ${OUT}.arg)"
bcftools index -f "${OUT}.bcf"
touch "${OUT}.done"
echo "DONE: ${OUT}.bcf"
