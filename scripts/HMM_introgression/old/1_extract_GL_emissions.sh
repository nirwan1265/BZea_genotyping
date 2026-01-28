#!/bin/bash
#BSUB -J "GL_bySample[1-1600]%200"
#BSUB -n 1
#BSUB -W 08:00
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=6GB]"
#BSUB -o /rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf/logs/GL_bySample.%J.%I.out
#BSUB -e /rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf/logs/GL_bySample.%J.%I.err


# DO this once:
#cd /rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf

#bcftools query -l chr1.teosnp.clean.rehead.bcf > samples.chr1.txt
#wc -l samples.chr1.txt
#head samples.chr1.txt



set -euo pipefail

BCFDIR="/rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf"
SAMPLELIST="${BCFDIR}/samples.chr1.txt"

OUTDIR="/rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf/HMM_inputs_GL"
LOGDIR="${BCFDIR}/logs"

DPmin=6   # change if you want

mkdir -p "${OUTDIR}" "${LOGDIR}"

SAMPLE="$(sed -n "${LSB_JOBINDEX}p" "${SAMPLELIST}")"
if [[ -z "${SAMPLE}" ]]; then
  echo "ERROR: No sample for index ${LSB_JOBINDEX} in ${SAMPLELIST}" >&2
  exit 2
fi

OUT="${OUTDIR}/${SAMPLE}.chr1-10.GL.DP${DPmin}plus.tsv.gz"
tmp="$(mktemp)"

echo "JOBINDEX=${LSB_JOBINDEX}"
echo "SAMPLE=${SAMPLE}"
echo "OUT=${OUT}"
echo "DPmin=${DPmin}"

for c in {1..10}; do
  CHR="chr${c}"
  BCF="${BCFDIR}/${CHR}.teosnp.clean.rehead.bcf"

  [[ -f "${BCF}" ]] || { echo "ERROR missing: ${BCF}" >&2; exit 2; }

  echo "Extracting ${CHR} from ${BCF}"

  # Output: CHROM POS REF ALT DP GL00 GL01 GL11
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

# global sort + compress
sort -k1,1 -k2,2n "${tmp}" | bgzip -c > "${OUT}"
rm -f "${tmp}"

echo "DONE: ${OUT}"









