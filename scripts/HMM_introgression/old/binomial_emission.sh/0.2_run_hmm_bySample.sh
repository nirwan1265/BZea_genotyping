#!/usr/bin/env bash
#BSUB -J "hmm_bySample[1-1600]%200"
#BSUB -n 1
#BSUB -W 08:00
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=6GB]"
#BSUB -o /rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf/logs/hmm_bySample.%J.%I.out
#BSUB -e /rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf/logs/hmm_bySample.%J.%I.err

set -euo pipefail

BCFDIR="/rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf"
INPDIR="${BCFDIR}/HMM_inputs_counts"
OUTDIR="${BCFDIR}/HMM_outputs_statepaths"
LOGDIR="${BCFDIR}/logs"

MAPDIR="/rsstu/users/r/rrellan/sara/ref/NAM_genetic_map"

mkdir -p "${OUTDIR}" "${LOGDIR}"

SAMPLELIST="${BCFDIR}/samples.chr1.txt"
[[ -f "${SAMPLELIST}" ]] || { echo "Missing ${SAMPLELIST}. Run counts job first." >&2; exit 2; }

SAMPLE="$(sed -n "${LSB_JOBINDEX}p" "${SAMPLELIST}")"
if [[ -z "${SAMPLE}" ]]; then
  echo "ERROR: No sample for JOBINDEX=${LSB_JOBINDEX} in ${SAMPLELIST}" >&2
  exit 2
fi

IN="${INPDIR}/${SAMPLE}.chr1-10.refalt.DP2plus.tsv.gz"
OUT="${OUTDIR}/${SAMPLE}.statepath.tsv.gz"

if [[ ! -f "${IN}" ]]; then
  echo "ERROR: Missing input counts file: ${IN}" >&2
  exit 2
fi

echo "JOBINDEX=${LSB_JOBINDEX}"
echo "SAMPLE=${SAMPLE}"
echo "IN=${IN}"
echo "OUT=${OUT}"

python3 hmm_binom_introgress.py "${IN}" "${MAPDIR}" "${OUT}" "BC2S3"

echo "DONE: ${OUT}"
