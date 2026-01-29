#!/bin/bash
#BSUB -J "GL_bySample[1428-1460]%50"
#BSUB -n 1
##BSUB -q sara
#BSUB -W 5000
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=6GB]"
#BSUB -o /rsstu/users/r/rrellan/BZea/angsd_genotyping/angsd_full_snps/logs/GL_bySample.%J.%I.out
#BSUB -e /rsstu/users/r/rrellan/BZea/angsd_genotyping/angsd_full_snps/logs/GL_bySample.%J.%I.err

set -euo pipefail

BCFDIR="/rsstu/users/r/rrellan/BZea/angsd_genotyping/angsd_full_snps/4_per_chr_renamed"
SAMPLELIST="${BCFDIR}/samples.chr1.txt"

OUTDIR="/rsstu/users/r/rrellan/BZea/angsd_genotyping/angsd_full_snps/5_HMM_inputs_GL"
LOGDIR="/rsstu/users/r/rrellan/BZea/angsd_genotyping/angsd_full_snps/logs"

mkdir -p "${OUTDIR}" "${LOGDIR}"

SAMPLE="$(sed -n "${LSB_JOBINDEX}p" "${SAMPLELIST}")"
if [[ -z "${SAMPLE}" ]]; then
  echo "ERROR: No sample for index ${LSB_JOBINDEX} in ${SAMPLELIST}" >&2
  exit 2
fi

OUT="${OUTDIR}/${SAMPLE}.chr1-10.GL.tsv.gz"
tmp="$(mktemp)"

echo "JOBINDEX=${LSB_JOBINDEX}"
echo "SAMPLE=${SAMPLE}"
echo "OUT=${OUT}"

for c in {1..10}; do
  CHR="chr${c}"
  BCF="${BCFDIR}/BZea_${CHR}.full.renamed.bcf"

  [[ -f "${BCF}" ]] || { echo "ERROR missing: ${BCF}" >&2; exit 2; }

  bcftools view -m2 -M2 -v snps -r "${CHR}" "${BCF}" -Ou \
  | bcftools query -s "${SAMPLE}" -f '%CHROM\t%POS\t%REF\t%ALT\t[%GL]\n' \
  | awk 'BEGIN{OFS="\t"}
      {
        gl=$5
        if(gl=="." || gl=="") next
        split(gl,a,",")
        if(a[1]=="" || a[2]=="" || a[3]=="") next
        print $1,$2,$3,$4,a[1],a[2],a[3]
      }' >> "${tmp}"
done

sort -k1,1 -k2,2n "${tmp}" | bgzip -c > "${OUT}"
rm -f "${tmp}"

echo "DONE: ${OUT}"
