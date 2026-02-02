#!/bin/bash
#BSUB -J "HMM_BC2S3[501-1600]%32"
#BSUB -n 1
#BSUB -W 5000
#BSUB -q sara
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=50GB]"
#BSUB -o /rsstu/users/r/rrellan/BZea/angsd_genotyping/angsd_full_snps/8_HMM_outputs/logs/hmm.%J.%I.out
#BSUB -e /rsstu/users/r/rrellan/BZea/angsd_genotyping/angsd_full_snps/8_HMM_outputs/logs/hmm.%J.%I.err

set -euo pipefail

GLDIR="/rsstu/users/r/rrellan/BZea/angsd_genotyping/angsd_full_snps/7_HMM_input_GL_filtered"
MAPDIR="/rsstu/users/r/rrellan/sara/ref/NAM_genetic_map"

OUTBASE="/rsstu/users/r/rrellan/BZea/angsd_genotyping/angsd_full_snps/8_HMM_outputs"
OUT_STATE="${OUTBASE}/statepaths"
OUT_TRACTS="${OUTBASE}/tracts"
LOGDIR="${OUTBASE}/logs"
CODEDIR="/share/maize/ntanduk/angsd_genotyping/HMM_introgression"

mkdir -p "${OUT_STATE}" "${OUT_TRACTS}" "${LOGDIR}"

# Do this once
#OUTBASE="/rsstu/users/r/rrellan/BZea/angsd_genotyping/angsd_full_snps/7_HMM_input_GL_filtered"
#GLDIR="${OUTBASE}"

#LIST="/rsstu/users/r/rrellan/BZea/angsd_genotyping/angsd_full_snps/8_HMM_outputs/gl_files.list"
#if [[ ! -f "${LIST}" ]]; then
#  ls -1 "${GLDIR}"/*.chr1-10.GL.filtered.tsv.gz | sort > "${LIST}"
#fi


LIST="${OUTBASE}/gl_files.list"


#if [[ ! -f "${LIST}" ]]; then
#  ls -1 "${GLDIR}"/*.chr1-10.GL.tsv.gz | sort > "${LIST}"
#fi

INFILE="$(sed -n "${LSB_JOBINDEX}p" "${LIST}")"
if [[ -z "${INFILE}" ]]; then
  echo "ERROR: No input file for JOBINDEX=${LSB_JOBINDEX} in ${LIST}" >&2
  exit 2
fi

BASE="$(basename "${INFILE}")"
SAMPLE="${BASE%%.chr1-10.GL.filtered.tsv.gz}"

OUT1="${OUT_STATE}/${SAMPLE}.statepath.tsv.gz"
OUT2="${OUT_TRACTS}/${SAMPLE}.tracts.bed.gz"

echo "JOBINDEX=${LSB_JOBINDEX}"
echo "INFILE=${INFILE}"
echo "SAMPLE=${SAMPLE}"
echo "OUT_STATE=${OUT1}"
echo "OUT_TRACTS=${OUT2}"

"${CODEDIR}/hmm_viterbi_bc2s3.py" \
  --in_gl_tsv_gz "${INFILE}" \
  --map_dir "${MAPDIR}" \
  --out_statepath_gz "${OUT1}" \
  --out_tracts_bed_gz "${OUT2}" \
  --prior_rr 0.85 --prior_rh 0.05 --prior_hh 0.10 \
  --min_morgan 1e-12 \
  --eta_hh_from_rh 0.25 \
  --rh_penalty 0 \
  --decode posterior_hysteresis \
  --rho 0.2 \
  --min_run_hh 3 \
  --min_run_rh 5

echo "DONE: ${SAMPLE}"
