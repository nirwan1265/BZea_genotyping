#!/bin/bash
#SBATCH --job-name=kmer_taxon
#SBATCH --output=logs/32b_kmer_taxon_profile.%j.out
#SBATCH --error=logs/32b_kmer_taxon_profile.%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --mem=32G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# 32b + 32c — teosinte k-mer profile per pool (sourmash sketches from 32a) and its PCA.
#   Submit after 32a:  sbatch 32b_kmer_taxon_profile.sh
#   Output: ${D32}/kmer_taxon_profile.tsv, ${D32}/pca/{kmer_pca.png, kmer_taxon_shares.png, kmer_taxon_agreement.tsv}
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"
need python3 Rscript
ASMS=(); for x in "${DONOR_ASSEMBLIES[@]}"; do IFS=: read -r n t _ <<< "$x"; ASMS+=("$n:$t"); done
python3 "${FINAL_DIR}/32b_kmer_taxon_profile.py" --refs "${D32}/refs" --pools "${D32}/pools" \
  --assemblies "${ASMS[@]}" --pools_tsv "${D11}/pools.tsv" --out "${D32}/kmer_taxon_profile.tsv"
Rscript "${FINAL_DIR}/32c_kmer_pca.R" "${D32}/kmer_taxon_profile.tsv" "${D32}/pca"
