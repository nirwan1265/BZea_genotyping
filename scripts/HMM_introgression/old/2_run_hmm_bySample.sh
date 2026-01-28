#!/usr/bin/env bash
#BSUB -J "HMM_BC2S3[1-1600]%200"
#BSUB -n 1
#BSUB -W 08:00
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=4GB]"
#BSUB -o /rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf/HMM_outputs/logs/hmm.%J.%I.out
#BSUB -e /rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf/HMM_outputs/logs/hmm.%J.%I.err

set -euo pipefail

GLDIR="/rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf/HMM_inputs_GL"
MAPDIR="/rsstu/users/r/rrellan/sara/ref/NAM_genetic_map"

OUTBASE="/rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf/HMM_outputs"
OUT_STATE="${OUTBASE}/statepaths"
OUT_TRACTS="${OUTBASE}/tracts"
LOGDIR="${OUTBASE}/logs"
CODEDIR="/share/maize/ntanduk/angsd_genotyping/HMM_introgression"

mkdir -p "${OUT_STATE}" "${OUT_TRACTS}" "${LOGDIR}"

# Build stable file list
LIST="${OUTBASE}/gl_files.list"
if [[ ! -f "${LIST}" ]]; then
  ls -1 "${GLDIR}"/*.chr1-10.GL.DP2plus.tsv.gz | sort > "${LIST}"
fi

INFILE="$(sed -n "${LSB_JOBINDEX}p" "${LIST}")"
if [[ -z "${INFILE}" ]]; then
  echo "ERROR: No input file for JOBINDEX=${LSB_JOBINDEX} in ${LIST}" >&2
  exit 2
fi

BASE="$(basename "${INFILE}")"
SAMPLE="${BASE%%.chr1-10.GL.DP2plus.tsv.gz}"

OUT1="${OUT_STATE}/${SAMPLE}.statepath.tsv.gz"
OUT2="${OUT_TRACTS}/${SAMPLE}.tracts.bed.gz"

echo "JOBINDEX=${LSB_JOBINDEX}"
echo "INFILE=${INFILE}"
echo "SAMPLE=${SAMPLE}"
echo "OUT_STATE=${OUT1}"
echo "OUT_TRACTS=${OUT2}"

# Old hmm
# "${CODEDIR}/hmm_viterbi_bc2s3.py" \
#  --in_gl_tsv_gz "${INFILE}" \
#  --map_dir "${MAPDIR}" \
#  --out_statepath_gz "${OUT1}" \
#  --out_tracts_bed_gz "${OUT2}" \
#  --prior_rr 0.92 --prior_rh 0.02 --prior_hh 0.06 \
#  --min_morgan 1e-8


"${CODEDIR}/hmm_viterbi_bc2s3.py" \
  --in_gl_tsv_gz "${INFILE}" \
  --map_dir "${MAPDIR}" \
  --out_statepath_gz "${OUT1}" \
  --out_tracts_bed_gz "${OUT2}" \
  --prior_rr 0.9 --prior_rh 0.03 --prior_hh 0.07 \
  --min_morgan 1e-8 \
  --eta_hh_from_rh 0.35 \
  --rh_penalty 0 \
  --rho 5 \
  --min_run_hh 1 \
  --min_run_rh 1


echo "DONE: ${SAMPLE}"
