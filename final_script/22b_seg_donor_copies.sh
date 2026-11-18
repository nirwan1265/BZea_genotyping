#!/bin/bash
#SBATCH --job-name=seg_copies
#SBATCH --output=logs/22b_seg_donor_copies.%j.out
#SBATCH --error=logs/22b_seg_donor_copies.%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --mem=16G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# 22b — Donor allele share inside SEG calls, against the Mendelian 1/2 (22b_seg_donor_copies.py).
#   Output in ${D21}: seg_donor_copies.tsv (per pool), seg_donor_copies.tsv.taxon.tsv (per taxon)
#   Needs python3 with numpy, pandas, scipy (e.g. conda activate /usr/local/usrapps/maize/ntanduk/bzea_gwas).
#   Submit after 21:  sbatch 22b_seg_donor_copies.sh
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"
need python3
python3 -c "import numpy, pandas, scipy" 2> /dev/null || die "python3 here has no numpy/pandas/scipy"
ASMS=$(for a in "${DONOR_ASSEMBLIES[@]}"; do IFS=: read -r n t _ <<< "$a"; echo "${n}:${t}"; done | paste -sd,)
python3 "${FINAL_DIR}/22b_seg_donor_copies.py" --hmm_dir "$D20" --ref_dir "$D19" \
  --pool_summary "${D21}/pool_summary.tsv" --assemblies "$ASMS" --plants "$KMER_PLANTS" \
  --out "${D21}/seg_donor_copies.tsv" --selftest
