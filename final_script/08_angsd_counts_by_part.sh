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
if angsd_finished "$OUT" && [[ -s ${OUT}.counts.gz ]]; then
  echo "ANGSD counts already finished; converting only."
else
  angsd \
    -bam "$BAMLIST" \
    -ref "$REF" \
    -r "${CHR}:${START}-${END}" \
    -sites "$SITES" \
    -doCounts 1 \
    -dumpCounts 4 \
    "${READ_FILTERS[@]}" \
    -P "$THREADS" \
    -out "$OUT"
  angsd_finished "$OUT" || die "ANGSD did not finish (no [ALL done] in ${OUT}.arg)"
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

if [[ ${KEEP_RAW_COUNTS} != 1 ]]; then rm -f "${OUT}.counts.gz"; fi
touch "${OUT}.done"
echo "DONE: ${OUT}.AD.bcf"
