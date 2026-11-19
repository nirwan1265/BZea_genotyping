#!/bin/bash
#SBATCH --job-name=targseq_call
#SBATCH --output=logs/34b_targseq_call.%A_%a.out
#SBATCH --error=logs/34b_targseq_call.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=08:00:00
#SBATCH --mem=8G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# 34b — Variable sites across the 96 targeted-sequencing samples, one array task per chromosome.
#   bcftools mpileup (MAPQ >= 20, base quality >= 20, allele depths AD) + multi-sample call (-mv), biallelic SNPs
#   with QUAL >= 30. Only regions with reads are visited, so this is fast for a targeted panel.
#   No donor genotype is needed: step 34c uses the allele depths per sample.
#   Output: ${D34}/vcf/targseq_<chr>.vcf.gz (+ .csi)
#   Needs bcftools. Submit after 34:  sbatch --array=1-10 34b_targseq_call.sh
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"
need bcftools
D34=${D34:-${D21}/targseq}
mkdir -p "$D34/vcf"
CHR=${CHROMS[$(( ${SLURM_ARRAY_TASK_ID:-1} - 1 ))]}
ls "$D34"/bam/*.bam > "$D34/vcf/bams.txt"
echo "$CHR | $(wc -l < "$D34/vcf/bams.txt") BAMs"
bcftools mpileup -f "$REF" -r "$CHR" -b "$D34/vcf/bams.txt" -q 20 -Q 20 -d 5000 -a FORMAT/AD,FORMAT/DP -Ou |
  bcftools call -mv -Ou |
  bcftools view -m2 -M2 -v snps -i 'QUAL>=30' -Oz -o "$D34/vcf/targseq_${CHR}.vcf.gz"
bcftools index "$D34/vcf/targseq_${CHR}.vcf.gz"
echo "$CHR: $(bcftools view -H "$D34/vcf/targseq_${CHR}.vcf.gz" | wc -l) SNPs"
