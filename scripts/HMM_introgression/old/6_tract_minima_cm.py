#!/bin/bash
#BSUB -n 1
#BSUB -W 14400
#BSUB -q sara 
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=20GB]"
#BSUB -J minima_cm
#BSUB -o stdout.%J
#BSUB -e stderr.%J




python3 tract_minima_cm.py \
  /rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf/HMM_outputs/statepaths \
  /rsstu/users/r/rrellan/sara/ref/NAM_genetic_map \
  /rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf/HMM_outputs/state_stats/tract_minima_cm.tsv
