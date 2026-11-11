#!/bin/bash
#BSUB -J "angsd_rerun_split2[1-60]%60"
#BSUB -n 4
#BSUB -W 5000
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=48GB]"
#BSUB -o /rsstu/users/r/rrellan/BZea/angsd_genotyping/angsd_full_snps/remaining/logs/angsd_rerun.%J.%I.out
#BSUB -e /rsstu/users/r/rrellan/BZea/angsd_genotyping/angsd_full_snps/remaining/logs/angsd_rerun.%J.%I.err

set -euo pipefail

# ---------------- USER EDITS ----------------
BAMLIST=/rsstu/users/r/rrellan/DOE_CAREER/BZea/dedup_filtered/bam_list.txt
REF=/rsstu/users/r/rrellan/sara/ref/Zm-B73-REFERENCE-NAM-5.0.fa

# this is the 60-line list we generated above
LIST=/share/maize/ntanduk/angsd_genotyping/HMM_introgression/logs/failed_jobs_subintervals.tsv

OUTDIR=/rsstu/users/r/rrellan/BZea/angsd_genotyping/angsd_full_snps/remaining

SNP_PVAL="1e-10"
MINMAPQ=30
MINQ=20
MININD=20
P=4
# -------------------------------------------

mkdir -p "${OUTDIR}/logs"

line="$(sed -n "${LSB_JOBINDEX}p" "${LIST}")"
if [[ -z "${line}" ]]; then
  echo "ERROR: no line for JOBINDEX=${LSB_JOBINDEX} in ${LIST}" >&2
  exit 2
fi

CHR="$(echo "$line" | awk '{print $1}')"
START="$(echo "$line" | awk '{print $2}')"
END="$(echo "$line" | awk '{print $3}')"
PARTNUM="$(echo "$line" | awk '{print $4}')"
SUBTAG="$(echo "$line" | awk '{print $5}')"

INTERVAL="${CHR}:${START}-${END}"
OUT="${OUTDIR}/BZea_${CHR}.part${PARTNUM}.${SUBTAG}.${START}_${END}"

echo "JOBINDEX=${LSB_JOBINDEX}"
echo "CHR=${CHR}"
echo "PARTNUM=${PARTNUM}"
echo "SUBTAG=${SUBTAG}"
echo "INTERVAL=${INTERVAL}"
echo "OUT=${OUT}"

# skip if already done + indexed
if [[ -f "${OUT}.bcf" ]]; then
  echo "FOUND existing ${OUT}.bcf"
  if [[ -f "${OUT}.bcf.csi" || -f "${OUT}.bcf.tbi" ]]; then
    echo "Index exists; skipping."
    exit 0
  fi
fi

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
  -setMinDepthInd 4 \
  -P "${P}" \
  -out "${OUT}"

# index for downstream bcftools query
bcftools index -f "${OUT}.bcf"

echo "DONE: ${OUT}.bcf"
