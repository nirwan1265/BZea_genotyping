#!/bin/bash
#BSUB -J "angsd_chr_parts[1-20]%20"
#BSUB -n 4
#BSUB -W 5000
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=16GB]"
#BSUB -o logs/angsd.%J.%I.out
#BSUB -e logs/angsd.%J.%I.err

set -euo pipefail

# ---------------- USER EDITS ----------------
BAMLIST=/rsstu/users/r/rrellan/DOE_CAREER/BZea/dedup_filtered/bam_list.txt
REF=/rsstu/users/r/rrellan/sara/ref/Zm-B73-REFERENCE-NAM-5.0.fa

# CHANGE THIS MANUALLY (no loop):
CHR="chr10"     # or "chr2", ... depending on your reference naming

# Output prefix directory
OUTDIR=/rsstu/users/r/rrellan/BZea/angsd_genotyping/angsd_full_snps
OUTPREFIX="BZea_${CHR}"

# SNP calling stringency + filters
SNP_PVAL="1e-8"
MINMAPQ=30
MINQ=20
MININD=5          # IMPORTANT: raise this if you want less sparse/noisy SNP set
P=4               # threads
# -------------------------------------------

mkdir -p "${OUTDIR}" logs

# Check .fai exists
if [[ ! -f "${REF}.fai" ]]; then
  echo "ERROR: Missing reference index ${REF}.fai (run: samtools faidx ${REF})" >&2
  exit 2
fi

# Get chromosome length from .fai
CHR_LEN=$(awk -v c="${CHR}" '$1==c{print $2}' "${REF}.fai")
if [[ -z "${CHR_LEN}" ]]; then
  echo "ERROR: CHR=${CHR} not found in ${REF}.fai" >&2
  exit 2
fi

# Split into 20 parts
NPARTS=20
PART=${LSB_JOBINDEX}

# ceil(CHR_LEN/NPARTS)
CHUNK=$(( (CHR_LEN + NPARTS - 1) / NPARTS ))

START=$(( (PART - 1) * CHUNK + 1 ))
END=$(( PART * CHUNK ))
if (( END > CHR_LEN )); then END=$CHR_LEN; fi

# In case last chunks are empty (shouldn't happen, but safe)
if (( START > CHR_LEN )); then
  echo "SKIP: PART=${PART} gives START=${START} > CHR_LEN=${CHR_LEN}"
  exit 0
fi

INTERVAL="${CHR}:${START}-${END}"
OUT="${OUTDIR}/${OUTPREFIX}.part${PART}.${START}_${END}"

echo "CHR=${CHR}  CHR_LEN=${CHR_LEN}"
echo "PART=${PART}/${NPARTS}  INTERVAL=${INTERVAL}"
echo "OUT=${OUT}"

# SNP calling + GL (BEAGLE)
angsd \
  -bam "${BAMLIST}" \
  -ref "${REF}" \
  -r "${INTERVAL}" \
  -GL 1 \
  -doMajorMinor 1 \
  -doMaf 1 \
  -doBcf 1 \
  -SNP_pval "${SNP_PVAL}" \
  -minMapQ "${MINMAPQ}" \
  -minQ "${MINQ}" \
  -uniqueOnly 1 \
  -remove_bads 1 \
  -only_proper_pairs 1 \
  -baq 1 \
  -C 50 \
  -minInd "${MININD}" \
  -skipTriallelic 1 \
  -setMinDepthInd 2 \
  -P "${P}" \
  -out "${OUT}"

echo "DONE: ${OUT}"
