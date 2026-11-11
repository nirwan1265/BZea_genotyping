#!/bin/bash
#BSUB -J "angsd_disc[1-200]%40"
#BSUB -n 4
#BSUB -W 5000
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=16GB]"
#BSUB -o logs/01_angsd_disc.%J.%I.out
#BSUB -e logs/01_angsd_disc.%J.%I.err
# =============================================================================
# 01 — ANGSD SNP discovery + genotype likelihoods, one job per chromosome part.
#      (replaces scripts/HMM_introgression/0.0_angsd_genotyping.sh; no manual CHR edits)
#
#   Array index 1..200 = 10 chromosomes x 20 parts.  chr1 = 1-20, chr2 = 21-40, ...
#   Output: ${D01}/BZea_<chr>.part<k>.<start>_<end>.bcf (+ .csi, .arg, .done)
#
#   STATUS: already run (outputs in .../angsd_full_snps). Kept for reproducibility.
#   NOTE:   ANGSD's BCF alleles are MAJOR,MINOR (not B73 REF,ALT). Step 07 re-polarizes.
#
#   Submit:  bsub < 01_angsd_discovery_by_part.sh
#   Test one part interactively:  bash 01_angsd_discovery_by_part.sh 1
# =============================================================================
set -euo pipefail
source "${LS_SUBCWD:-$PWD}/config.sh"

IDX=${LSB_JOBINDEX:-${1:?usage: bash $0 <index 1-200>}}
read -r CHR PART <<< "$(index_to_chr_part "$IDX")"
read -r START END <<< "$(part_interval "$CHR" "$PART")"

mkdir -p "$D01"
OUT=${D01}/BZea_${CHR}.part${PART}.${START}_${END}
echo "CHR=${CHR} PART=${PART}/${NPARTS} INTERVAL=${CHR}:${START}-${END}"
echo "OUT=${OUT}"

[[ -f ${OUT}.done ]] && { echo "Already done; skipping."; exit 0; }
[[ -f ${REF}.fai ]] || die "missing ${REF}.fai (samtools faidx ${REF})"

angsd \
  -bam "$BAMLIST" \
  -ref "$REF" \
  -r "${CHR}:${START}-${END}" \
  -GL 1 \
  -doMajorMinor 1 \
  -doMaf 1 \
  -doBcf 1 \
  -SNP_pval "$SNP_PVAL" \
  "${READ_FILTERS[@]}" \
  -minInd "$MININD" \
  -skipTriallelic 1 \
  -setMinDepthInd "$SETMINDEPTHIND" \
  -P "$THREADS" \
  -out "$OUT"

angsd_finished "$OUT" || die "ANGSD did not finish (no [ALL done] in ${OUT}.arg)"
bcftools index -f "${OUT}.bcf"
touch "${OUT}.done"
echo "DONE: ${OUT}.bcf"
