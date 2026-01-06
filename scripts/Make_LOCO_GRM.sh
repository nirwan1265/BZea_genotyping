BASE="/Users/nirwantandukar/Documents/Research/data/BZea/genotype/BZea.beagle.imputed.allchr.renamed"
OUT="/Users/nirwantandukar/Documents/Research/data/BZea/genotype/LOCO_K"
mkdir -p "${OUT}/tmp"

# 1) LD prune ONCE (recommended)
# Tune these if you want (window/step/r2). This is a common default.
plink2 \
  --bfile "${BASE}" \
  --maf 0.01 \
  --geno 0.2 \
  --indep-pairwise 200 50 0.2 \
  --out "${OUT}/tmp/prune"

# 2) LOCO GRM per chromosome: exclude SNPs on that chromosome, compute rel matrix
for c in {1..10}; do
  echo "LOCO K for chr${c}"

  # list SNP IDs on chr c (works if BIM has "1".."10" OR "chr1".."chr10")
  awk -v c="${c}" '($1==c || $1=="chr"c){print $2}' "${BASE}.bim" > "${OUT}/tmp/chr${c}.snps"

  # compute GRM using all OTHER chromosomes
  plink2 \
    --bfile "${BASE}" \
    --extract "${OUT}/tmp/prune.prune.in" \
    --exclude "${OUT}/tmp/chr${c}.snps" \
    --make-rel square \
    --out "${OUT}/tmp/LOCO_chr${c}"
done

TMP="/Users/nirwantandukar/Documents/Research/data/BZea/genotype/LOCO_K/tmp"
OUT="/Users/nirwantandukar/Documents/Research/data/BZea/genotype/LOCO_K"

mkdir -p "${OUT}"

for c in {1..10}; do
  awk '{print $2}' "${TMP}/LOCO_chr${c}.rel.id" > "${TMP}/LOCO_chr${c}.iid"
  paste "${TMP}/LOCO_chr${c}.iid" "${TMP}/LOCO_chr${c}.rel" > "${OUT}/K_matrix_chr${c}.txt"
  echo "Wrote ${OUT}/K_matrix_chr${c}.txt"
done




Rscript GridLMM_BZea_RTIGER_GWAS.R \
  rtiger_stateprob 1 8 \
  /Users/nirwantandukar/Documents/Github/BZea_genotyping/data/phenotypes/BZea_FloweringTime_spats_BLUE_weighted.csv \
  DTS 1 \
  /Users/nirwantandukar/Documents/Research/data/BZea/genotype/LOCO_K \
  /Users/nirwantandukar/Documents/Research/data/BZea/genotype/GridLMM_out \
  /Users/nirwantandukar/Documents/Research/data/BZea/rtiger_results/Samples_1_100_RTIGER_results.rds \
  /Users/nirwantandukar/Documents/Research/data/BZea/rtiger_inputs_renamed \
  100000 \
  /Users/nirwantandukar/Documents/Research/data/BZea/genotype/chr_len.rds
