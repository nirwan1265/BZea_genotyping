#!/bin/bash
#SBATCH --job-name=assoc_inputs
#SBATCH --output=logs/25_assoc_inputs.%j.out
#SBATCH --error=logs/25_assoc_inputs.%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --mem=8G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# 25 — Association inputs from the k-mer calls (25_assoc_inputs.py): window dosage (0/1/2),
#   SEG indicator, LOCO kinship per chromosome, phenotype joined by Sequencing_ID.
#   Low-depth (NOISY) pools are dropped. Output: ${D25}
#   Submit after 21:  sbatch 25_assoc_inputs.sh
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"
need python3
python3 "${FINAL_DIR}/25_assoc_inputs.py" --tracts "${D21}/tracts_all.tsv.gz" --pool_summary "${D21}/pool_summary.tsv" \
  --pheno "$PHENO" --out_dir "$D25" --min_carriers "$ASSOC_MIN_CARRIERS" \
  $( [[ "${MARK_GENES:-0}" == 1 ]] && echo --genes "$GENES_TSV" --gff "$GFF_GLOB" )
[[ "${MARK_GENES:-0}" == 1 ]] || rm -f "$D25/gene_positions.tsv"     # no gene labels on the plots
