#!/bin/bash
# filter_angsd_bcf_noDP.sh
# Filter ANGSD BCFs when DP field is completely missing

BCFDIR="/rsstu/users/r/rrellan/BZea/angsd_genotyping/angsd_full_snps/4_per_chr_renamed"
FILTERED_DIR="${BCFDIR}/filtered_noDP"
mkdir -p "${FILTERED_DIR}"

REFERENCE="/rsstu/users/r/rrellan/sara/ref/Zm-B73-REFERENCE-NAM-5.0.fa"

for c in {1..10}; do
  CHR="chr${c}"
  BCF="${BCFDIR}/BZea_${CHR}.full.renamed.bcf"
  OUT="${FILTERED_DIR}/BZea_${CHR}.filtered.bcf"
  
  echo "Filtering ${CHR}..."
  
  # Step-by-step filtering for ANGSD BCFs
  bcftools view "${BCF}" -Ou \
    | bcftools filter -e 'QUAL<20' -Ou 2>/dev/null \
    | bcftools filter -e 'F_MISSING>0.2' -Ou \
    | bcftools view -m2 -M2 -v snps -Ou \
    | bcftools +fill-tags -Ou -- -t AF,NS,ExcHet \
    | bcftools filter -e 'AF<0.05 || AF>0.95' -Ou \
    | bcftools filter -e 'INFO/NS<100' -Ou \
    | bcftools norm --check-ref w -f "${REFERENCE}" -Ou \
    | bcftools view -Oz -o "${OUT}"
  
  bcftools index "${OUT}"
  
  # Count before/after
  BEFORE=$(bcftools view -H "${BCF}" | wc -l)
  AFTER=$(bcftools view -H "${OUT}" | wc -l)
  echo "  ${CHR}: ${BEFORE} → ${AFTER} SNPs (kept $((AFTER*100/BEFORE))%)"
done
