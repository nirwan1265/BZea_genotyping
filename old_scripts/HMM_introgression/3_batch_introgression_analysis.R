#!/bin/bash
#BSUB -n 1
#BSUB -W 14400
#BSUB -q sara 
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=20GB]"
#BSUB -J introgression_stats
#BSUB -o stdout.%J
#BSUB -e stderr.%J


# make the output directory /rsstu/users/r/rrellan/BZea/angsd_genotyping/angsd_full_snps/8_HMM_outputs/stats

Rscript batch_introgression_analysis.R /rsstu/users/r/rrellan/BZea/angsd_genotyping/angsd_full_snps/8_HMM_outputs/statepaths 10000 5000 /rsstu/users/r/rrellan/BZea/angsd_genotyping/angsd_full_snps/8_HMM_outputs/stats
