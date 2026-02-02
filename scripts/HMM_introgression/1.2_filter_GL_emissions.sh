#!/bin/bash
#==============================================================================
# Filter GL TSVs: remove (0,0,0) and "nearly uninformative" sites
#
# Input:  /rsstu/users/r/rrellan/BZea/angsd_genotyping/angsd_full_snps/6_HMM_inputs_GL/*.GL.tsv.gz
# Output: /rsstu/users/r/rrellan/BZea/angsd_genotyping/angsd_full_snps/7_HMM_input_GL_filtered/*.GL.filtered.tsv.gz
#
# "Nearly uninformative" definition:
#   max(GL) - min(GL) < EPS  (default EPS=0.05 log10 units)
#
# LSF job array: 1 job per sample
#==============================================================================

#BSUB -n 1
#BSUB -W 06:00
#BSUB -q sara
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=4GB]"
#BSUB -J "GLfilter[201-1600]%32"
#BSUB -o logs/stdout.%J.%I
#BSUB -e logs/stderr.%J.%I

set -euo pipefail

IN_DIR="/rsstu/users/r/rrellan/BZea/angsd_genotyping/angsd_full_snps/6_HMM_inputs_GL"
OUT_DIR="/rsstu/users/r/rrellan/BZea/angsd_genotyping/angsd_full_snps/7_HMM_input_GL_filtered"
EPS="${EPS:-0.05}"   # log10 GL range threshold; override by exporting EPS=0.1 etc.

mkdir -p "${OUT_DIR}" "${OUT_DIR}/logs"

# Build file list deterministically
FILELIST="${OUT_DIR}/gl_files.list"
if [[ ! -s "${FILELIST}" ]]; then
  ls -1 "${IN_DIR}"/*.GL.tsv.gz | sort > "${FILELIST}"
fi

# Pick file for this array index
F=$(sed -n "${LSB_JOBINDEX}p" "${FILELIST}" || true)
if [[ -z "${F}" ]]; then
  echo "No input file at index ${LSB_JOBINDEX}. Check array range vs file count." >&2
  echo "N files = $(wc -l < "${FILELIST}")" >&2
  exit 1
fi

B=$(basename "${F}")
SAMPLE="${B%.chr1-10.GL.tsv.gz}"
OUT="${OUT_DIR}/${SAMPLE}.chr1-10.GL.filtered.tsv.gz"

echo "== GL filter =="
echo "Input : ${F}"
echo "Output: ${OUT}"
echo "EPS   : ${EPS}"

# Filter logic:
#   - drop if GL5=GL6=GL7=0 (exact)
#   - drop if (max-min) < EPS (nearly uninformative)
# NOTE: fields 5,6,7 are GL for RR, RH, HH in your TSV.
zcat "${F}" \
| awk -F'\t' -v OFS='\t' -v eps="${EPS}" '
  NF>=7 {
    a=$5; b=$6; c=$7;

    # drop exact uninformative 0,0,0
    if (a==0 && b==0 && c==0) next;

    # compute range
    max=a; if(b>max)max=b; if(c>max)max=c;
    min=a; if(b<min)min=b; if(c<min)min=c;

    # drop nearly uninformative (flat likelihoods)
    if ((max-min) < eps) next;

    print $0
  }
' \
| gzip -c > "${OUT}"

# Optional: quick stats
IN_N=$(zcat "${F}" | wc -l)
OUT_N=$(zcat "${OUT}" | wc -l)
REM=$((IN_N - OUT_N))

echo "Lines in : ${IN_N}"
echo "Lines out: ${OUT_N}"
echo "Removed  : ${REM}"
echo "Done."

