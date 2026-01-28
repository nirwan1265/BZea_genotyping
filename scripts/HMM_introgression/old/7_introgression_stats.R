#!/bin/bash
#BSUB -n 1
#BSUB -W 14400
#BSUB -q sara 
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=20GB]"
#BSUB -J Intro_stats
#BSUB -o stdout.%J
#BSUB -e stderr.%J




Rscript introgression_stats.R
