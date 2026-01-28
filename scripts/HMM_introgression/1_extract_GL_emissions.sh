#!/bin/bash
#BSUB -J "GL_bySample[1-1600]%200"
#BSUB -n 1
#BSUB -q sara
#BSUB -W 5000
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=6GB]"
#BSUB -o /rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf/logs/GL_bySample.%J.%I.out
#BSUB -e /rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf/logs/GL_bySample.%J.%I.err

set -euo pipefail

BCFDIR="/rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf"
SAMPLELIST="${BCFDIR}/samples.chr1.txt"

# IMPORTANT: if you want multiple DP versions, make DP-specific OUTDIRs instead
DPmin=8
OUTDIR="/rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf/HMM_inputs_GL"   # or HMM_inputs_GL_DP6
LOGDIR="${BCFDIR}/logs"

mkdir -p "${OUTDIR}" "${LOGDIR}"

SAMPLE="$(sed -n "${LSB_JOBINDEX}p" "${SAMPLELIST}")"
if [[ -z "${SAMPLE}" ]]; then
  echo "ERROR: No sample for index ${LSB_JOBINDEX} in ${SAMPLELIST}" >&2
  exit 2
fi

# NO DP in name anymore
OUT="${OUTDIR}/${SAMPLE}.chr1-10.GL.tsv.gz"
tmp="$(mktemp)"

echo "JOBINDEX=${LSB_JOBINDEX}"
echo "SAMPLE=${SAMPLE}"
echo "OUT=${OUT}"
echo "DPmin=${DPmin}"

for c in {1..10}; do
  CHR="chr${c}"
  BCF="${BCFDIR}/${CHR}.teosnp.clean.rehead.bcf"

  [[ -f "${BCF}" ]] || { echo "ERROR missing: ${BCF}" >&2; exit 2; }

  bcftools view -m2 -M2 -v snps -r "${CHR}" "${BCF}" -Ou \
  | bcftools query -s "${SAMPLE}" -f '%CHROM\t%POS\t%REF\t%ALT\t[%DP]\t[%GL]\n' \
  | awk -v DPmin="${DPmin}" 'BEGIN{OFS="\t"}
      {
        split($6,a,","); dp=$5+0;
        if(dp>=DPmin && length($6)>0){
          print $1,$2,$3,$4,dp,a[1],a[2],a[3]
        }
      }' >> "${tmp}"
done

sort -k1,1 -k2,2n "${tmp}" | bgzip -c > "${OUT}"
rm -f "${tmp}"

echo "DONE: ${OUT}"
