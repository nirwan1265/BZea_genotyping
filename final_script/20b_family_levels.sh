#!/bin/bash
#SBATCH --job-name=kmer_famlev
#SBATCH --output=logs/20b_family_levels.%j.out
#SBATCH --error=logs/20b_family_levels.%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --mem=32G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# 20b — Family DONOR levels for round ROUND, from round ROUND-1 (20b_family_levels.py).
#   Reads the previous round's step-20 .npz and step-21 rho_loglik.tsv; writes ${D20B}/<SAMPLE>.npy
#   (leave-one-out family level per pool) and summary.tsv (per family).
#   Round 1 (after the assembly run 20-21):  sbatch --export=ALL,ROUND=1 20b_family_levels.sh
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"
need python3
(( ROUND >= 1 )) || die "set ROUND >= 1 (sbatch --export=ALL,ROUND=1 20b_family_levels.sh)"
P=$(( ROUND - 1 ))
if (( P == 0 )); then PREV_HMM=${KMER_DIR}/hmm; PREV_OUT=${KMER_DIR}/introgression
else PREV_HMM=${KMER_DIR}/hmm_fam${P}; PREV_OUT=${KMER_DIR}/introgression_fam${P}; fi
[[ -s ${PREV_OUT}/rho_loglik.tsv ]] || die "missing ${PREV_OUT}/rho_loglik.tsv (run round ${P} steps 20-21 first)"
ASMS=$(for a in "${DONOR_ASSEMBLIES[@]}"; do IFS=: read -r n t _ <<< "$a"; echo "${n}:${t}"; done | paste -sd,)
python3 "${FINAL_DIR}/20b_family_levels.py" --hmm_dir "$PREV_HMM" --rho_file "${PREV_OUT}/rho_loglik.tsv" \
  --meta "$SAMPLE_META" --ref_dir "$D19" --assemblies "$ASMS" --out_dir "$D20B" --noisy_sigma "$KMER_NOISY_SIGMA"
