#!/bin/bash
#SBATCH --job-name=angsd_counts
#SBATCH --array=1-200%40
#SBATCH --output=logs/08_counts.%A_%a.out
#SBATCH --error=logs/08_counts.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=3-11:20:00
#SBATCH --mem=32G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# 08 — Per-pool B73/ALT read counts at the step-07 sites, one job per part.
#
#   1) ANGSD -doCounts 1 -dumpCounts 4 at the discovered sites only, with the
#      SAME read filters as discovery. No -minInd / -setMinDepthInd: every read
#      that passed the filters is counted, and zero-read pools stay as 0,0.
#   2) A/C/G/T per pool -> AD (B73,ALT) in a BCF, plus per-site and per-pool
#      summaries for QC and step-2 filtering.
#
#   Array index 1..200 = 10 chromosomes x 20 parts (same split as step 01).
#   Output in ${D08}, prefix BZea_<chr>.part<k>.<start>_<end>:
#     .AD.bcf (+ .csi)     B73/ALT read counts per pool (REF = B73 base)
#     .sitesum.tsv.gz      per-site totals incl. B73-control reads
#     .sampleqc.tsv        per-pool totals for this part
#     .pos.gz .arg .done   ANGSD bookkeeping (.counts.gz deleted unless KEEP_RAW_COUNTS=1)
#
#   Submit after 07 finishes:  sbatch --dependency=afterok:<07 job id> 08_angsd_counts_by_part.sh
#   Test one part:             bash 08_angsd_counts_by_part.sh 1
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"
need angsd bcftools python3

IDX=${SLURM_ARRAY_TASK_ID:-${1:?usage: bash $0 <index 1-200>}}
read -r CHR PART <<< "$(index_to_chr_part "$IDX")"
read -r START END <<< "$(part_interval "$CHR" "$PART")"

SITES=${D07}/BZea_${CHR}.sites
ALLELES=${D07}/BZea_${CHR}.alleles.tsv.gz
SAMPLES=${D07}/samples.txt
HEADER=${D07}/ad_header.vcf
CTRL=${D07}/b73_controls.idx
AWK=${FINAL_DIR}/08_counts_to_ad.awk
for f in "$SITES" "${SITES}.idx" "${SITES}.bin" "$ALLELES" "$SAMPLES" "$HEADER" "$CTRL" "$AWK"; do
  [[ -e $f ]] || die "missing ${f} (run 07 first, from inside final_script/)"
done

mkdir -p "$D08"
OUT=${D08}/BZea_${CHR}.part${PART}.${START}_${END}
NIND=$(wc -l < "$SAMPLES")
echo "CHR=${CHR} PART=${PART}/${NPARTS} INTERVAL=${CHR}:${START}-${END} POOLS=${NIND}"
echo "OUT=${OUT}"

[[ -f ${OUT}.done ]] && { echo "Already done; skipping."; exit 0; }

NSITES=$(awk -v s="$START" -v e="$END" '$2 >= s && $2 <= e' "$SITES" | wc -l)
echo "sites in interval: ${NSITES}"

if (( NSITES == 0 )); then
  # nothing to count: write header-only outputs so step 09 sees a complete part
  bcftools view --no-version -Ob -o "${OUT}.AD.bcf" "$HEADER"
  bcftools index -f "${OUT}.AD.bcf"
  printf "CHROM\tPOS\tB73\tALT\tDP\tOTH\tREF_READS\tALT_READS\tNCOV\tNALT\tCTRL_REF\tCTRL_ALT\tCTRL_NCOV\n" | gzip -c > "${OUT}.sitesum.tsv.gz"
  { printf "SAMPLE\tSITES_COVERED\tREF_READS\tALT_READS\tOTHER_READS\n"; awk '{print $0"\t0\t0\t0\t0"}' "$SAMPLES"; } > "${OUT}.sampleqc.tsv"
  touch "${OUT}.done"
  echo "DONE (empty part)."
  exit 0
fi

# ---- 1) ANGSD read counts ----------------------------------------------------
# ANGSD keeps every BAM open at once. Raise the open-file limit to the node's
# hard limit; if that is still too low for all BAMs, run ANGSD on consecutive
# slices of the BAM list and merge (08_merge_batches.py), keeping column order.
run_angsd_counts() {   # $1 = BAM list, $2 = output prefix
  angsd \
    -bam "$1" \
    -ref "$REF" \
    -r "${CHR}:${START}-${END}" \
    -sites "$SITES" \
    -doCounts 1 \
    -dumpCounts 4 \
    "${READ_FILTERS[@]}" \
    -P "$THREADS" \
    -out "$2"
  angsd_finished "$2" || die "ANGSD did not finish (no [ALL done] in $2.arg)"
}

if [[ -f ${OUT}.counts_ok && -s ${OUT}.counts.gz ]]; then
  echo "ANGSD counts already finished; converting only."
else
  BAMS_CLEAN=${OUT}.bams.txt
  tr -d '\r' < "$BAMLIST" | grep -v '^[[:space:]]*$' > "$BAMS_CLEAN"
  (( $(wc -l < "$BAMS_CLEAN") == NIND )) || die "BAM list and samples.txt have different lengths"

  ulimit -n "$(ulimit -Hn)" 2> /dev/null || true
  FDLIM=$(ulimit -n)
  if [[ $FDLIM == unlimited ]]; then PER=$NIND; else PER=$(( (FDLIM - 100) / 2 )); fi
  if [[ -n ${BAMS_PER_RUN:-} ]] && (( BAMS_PER_RUN < PER )); then PER=$BAMS_PER_RUN; fi
  (( PER >= 1 )) || die "open-file limit ${FDLIM} is too low to open any BAMs"
  if (( PER > NIND )); then PER=$NIND; fi
  NB=$(( (NIND + PER - 1) / PER ))
  echo "open-file limit: ${FDLIM} -> ${NB} ANGSD run(s) of <= ${PER} BAMs"

  if (( NB == 1 )); then
    run_angsd_counts "$BAMS_CLEAN" "$OUT"
  else
    split -d -a 3 -l "$PER" "$BAMS_CLEAN" "${OUT}.batch"
    merge_args=()
    for (( k = 0; k < NB; k++ )); do
      B=$(printf "%s.batch%03d" "$OUT" "$k")
      echo "-- batch $((k + 1))/${NB}: $(wc -l < "$B") BAMs"
      if angsd_finished "$B" && [[ -s ${B}.counts.gz ]]; then
        echo "   already finished"
      else
        run_angsd_counts "$B" "$B"
      fi
      merge_args+=(--batch "$B" "$(wc -l < "$B")")
    done
    python3 "${FINAL_DIR}/08_merge_batches.py" --out "$OUT" "${merge_args[@]}"
    rm -f "${OUT}".batch???.{pos.gz,counts.gz,arg} "${OUT}".batch???
  fi
  rm -f "$BAMS_CLEAN"
  touch "${OUT}.counts_ok"
fi

# ---- 2) A/C/G/T -> AD (B73,ALT) ---------------------------------------------
{
  cat "$HEADER"
  paste <(gzip -dc "${OUT}.pos.gz" | tail -n +2) <(gzip -dc "${OUT}.counts.gz" | tail -n +2) \
    | awk -v alleles="$ALLELES" -v start="$START" -v end="$END" -v nind="$NIND" \
          -v samples="$SAMPLES" -v ctrl="$CTRL" \
          -v sampleqc="${OUT}.sampleqc.tsv" -v sitesum="${OUT}.sitesum.tsv.gz" \
          -f "$AWK"
} | bcftools view --no-version -Ob -o "${OUT}.AD.bcf"
bcftools index -f "${OUT}.AD.bcf"

# ---- checks ------------------------------------------------------------------
N_ANGSD=$(( $(gzip -dc "${OUT}.pos.gz" | wc -l) - 1 ))
N_BCF=$(bcftools index -n "${OUT}.AD.bcf")
N_SUM=$(( $(gzip -dc "${OUT}.sitesum.tsv.gz" | wc -l) - 1 ))
echo "sites requested: ${NSITES} ; with reads (ANGSD): ${N_ANGSD} ; written: ${N_BCF}"
(( N_BCF == N_SUM )) || die "BCF (${N_BCF}) and site summary (${N_SUM}) disagree"
(( N_BCF == N_ANGSD )) || echo "WARNING: $(( N_ANGSD - N_BCF )) ANGSD sites not in the allele table"
[[ $(( $(wc -l < "${OUT}.sampleqc.tsv") - 1 )) -eq $NIND ]] || die "sample QC has wrong number of rows"

if [[ ${KEEP_RAW_COUNTS} != 1 ]]; then rm -f "${OUT}.counts.gz" "${OUT}.counts_ok"; fi
touch "${OUT}.done"
echo "DONE: ${OUT}.AD.bcf"
