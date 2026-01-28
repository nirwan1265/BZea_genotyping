#!/bin/bash
#BSUB -n 1
#BSUB -W 14400
#BSUB -q sara 
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=20GB]"
#BSUB -J combine
#BSUB -o stdout.%J
#BSUB -e stderr.%J





cd /rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf/HMM_outputs/state_stats
head -n 1 *.stats.tsv | head -n 1 > all_samples.state_stats.tsv
tail -n +2 -q *.stats.tsv >> all_samples.state_stats.tsv
