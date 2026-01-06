# Getting the allele table for the GridLMM
# Option 1 - Using VCF file# 2) Per-chromosome DS matrix (TSV.gz)

for c in chr{1..10}; do
  bcftools view -r ${c} -m2 -M2 -v snps BZea.beagle.imputed.allchr.renamed.vcf.gz \
    | bcftools query -f'%CHROM\t%POS[\t%DS]\n' \
    | bgzip -c > ${c}.DS.tsv.gz
  tabix -s1 -b2 -e2 ${c}.DS.tsv.gz
done

# Option 2 - Using PLINK files
for c in {1..10}; do
  plink2 --bfile BZea.beagle.imputed.allchr.renamed \
    --chr $c \
    --export A \
    --out chr${c}.A
done



# Get sample IDs form vcf file
bcftools query -l BZea.beagle.imputed.allchr.renamed.vcf.gz > samples.txt