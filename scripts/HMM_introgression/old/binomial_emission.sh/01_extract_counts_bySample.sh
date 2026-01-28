#!/usr/bin/env bash
#BSUB -J "counts_bySample[1-1600]%200"
#BSUB -n 1
#BSUB -W 08:00
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=6GB]"
#BSUB -o /rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf/logs/counts_bySample.%J.%I.out
#BSUB -e /rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf/logs/counts_bySample.%J.%I.err

set -euo pipefail

BCFDIR="/rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf"
OUTDIR="${BCFDIR}/HMM_inputs_counts"
LOGDIR="${BCFDIR}/logs"

mkdir -p "${OUTDIR}" "${LOGDIR}"

DPmin=2

# Sample list from chr1 header (THIS avoids your earlier “sample not found” issue)
SAMPLELIST="${BCFDIR}/samples.chr1.txt"
if [[ ! -f "${SAMPLELIST}" ]]; then
  bcftools query -l "${BCFDIR}/chr1.teosnp.clean.rehead.bcf" > "${SAMPLELIST}"
fi

SAMPLE="$(sed -n "${LSB_JOBINDEX}p" "${SAMPLELIST}")"
if [[ -z "${SAMPLE}" ]]; then
  echo "ERROR: No sample for JOBINDEX=${LSB_JOBINDEX} in ${SAMPLELIST}" >&2
  exit 2
fi

OUT="${OUTDIR}/${SAMPLE}.chr1-10.refalt.DP${DPmin}plus.tsv.gz"

echo "JOBINDEX=${LSB_JOBINDEX}"
echo "SAMPLE=${SAMPLE}"
echo "OUT=${OUT}"
echo "DPmin=${DPmin}"

tmp="$(mktemp)"

for c in {1..10}; do
  CHR="chr${c}"
  BCF="${BCFDIR}/${CHR}.teosnp.clean.rehead.bcf"
  [[ -f "${BCF}" ]] || { echo "Missing BCF: ${BCF}" >&2; exit 2; }

  echo "Extracting ${CHR} from ${BCF}"

  # DP4 = refF,refR,altF,altR
  bcftools query -s "${SAMPLE}" -f '%CHROM\t%POS\t%REF\t%ALT\t[%DP4]\n' "${BCF}" \
  | awk -F'\t' -v DPmin="${DPmin}" 'BEGIN{OFS="\t"}
      {
        split($5,a,",");
        refc=(a[1]+0)+(a[2]+0);
        altc=(a[3]+0)+(a[4]+0);
        dp=refc+altc;
        if(dp>=DPmin) print $1,$2,$3,$4,refc,altc,dp;
      }' >> "${tmp}"
done

# Global sort, compress, index
sort -k1,1 -k2,2n "${tmp}" | bgzip -c > "${OUT}"
tabix -f -s1 -b2 -e2 "${OUT}"
rm -f "${tmp}"

echo "DONE: ${OUT}"
