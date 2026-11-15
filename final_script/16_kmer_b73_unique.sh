#!/bin/bash
#SBATCH --job-name=kmer_b73
#SBATCH --output=logs/16_kmer_b73.%j.out
#SBATCH --error=logs/16_kmer_b73.%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --mem=48G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# 16 — B73 reference k-mers for the k-mer windows (run once).
#
#   1. ${KMER_DIR}/B73_chr1-10.fa : chr1-10 of the B73 v5 reference (windows are made on this)
#   2. ${KMER_DIR}/B73            : all 31-mers of the whole B73 genome (KMC)
#   3. ${KMER_DIR}/B73_uniq       : the 31-mers that occur exactly once in B73 (single-copy).
#      Repeats (most of the maize genome) are present in every pool, donor or not, and hide
#      the signal; only single-copy k-mers are used. (Pilot: 680 M of 817 M k-mers.)
#   Skips anything that already exists.
#
#   Submit with the k-mer env active:  sbatch 16_kmer_b73_unique.sh
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"
need kmc kmc_tools samtools

mkdir -p "$KMER_DIR"; cd "$KMER_DIR"
TMP=${KMER_DIR}/tmp_b73; mkdir -p "$TMP"
trap 'rm -rf "$TMP"' EXIT

[[ -s B73_chr1-10.fa ]] || samtools faidx "$REF" "${CHROMS[@]}" > B73_chr1-10.fa
[[ -s B73.kmc_pre ]] || kmc -k"$KMER_K" -m40 -t"${SLURM_CPUS_PER_TASK:-8}" -ci1 -fm "$REF" B73 "$TMP"
[[ -s B73_uniq.kmc_pre ]] || kmc_tools transform B73 -ci1 -cx1 reduce B73_uniq
kmc_tools info B73_uniq
echo "WROTE: ${KMER_DIR}/B73_chr1-10.fa, B73.kmc_*, B73_uniq.kmc_*"
