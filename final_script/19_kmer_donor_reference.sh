#!/bin/bash
#SBATCH --job-name=kmer_donor_ref
#SBATCH --output=logs/19_kmer_donor_ref.%A_%a.out
#SBATCH --error=logs/19_kmer_donor_ref.%A_%a.err
#SBATCH --array=1-9
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=06:00:00
#SBATCH --mem=48G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# 19 — Expected DONOR level per window, from the teosinte genome assemblies.
#
#   Task 1: B73 itself -> how many B73 single-copy k-mers each 100-kb window has
#           (${D19}/B73.tsv.gz; observed_kmers = all of them).
#   Task 2..: one assembly from DONOR_ASSEMBLIES (config.sh) -> how many of those k-mers
#           the assembly has (${D19}/<NAME>.tsv.gz).
#   Per window, d = assembly observed_kmers / B73 observed_kmers = share of B73 single-copy
#   k-mers a fixed donor tract of that taxon keeps. This is the DONOR level used by the HMM
#   (step 20), instead of a fitted d. Same windows and k-mers as step 17.
#
#   Submit after 16 (k-mer env active):  sbatch 19_kmer_donor_reference.sh
#   (array 1-9 = B73 + the 8 assemblies; change --array if DONOR_ASSEMBLIES changes)
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"
need kmc kmc_tools kcftools awk gzip

IDX=${SLURM_ARRAY_TASK_ID:-${1:?usage: bash $0 <task index>}}
CPUS=${SLURM_CPUS_PER_TASK:-8}
cd "$KMER_DIR"
[[ -s B73_uniq.kmc_pre && -s B73_chr1-10.fa ]] || die "missing B73_uniq / B73_chr1-10.fa in ${KMER_DIR} (run 16)"
mkdir -p "$D19"

if (( IDX == 1 )); then
  NAME=B73; FA=""
else
  (( IDX - 2 < ${#DONOR_ASSEMBLIES[@]} )) || { echo "task ${IDX}: no assembly"; exit 0; }
  IFS=: read -r NAME _ FA <<< "${DONOR_ASSEMBLIES[$(( IDX - 2 ))]}"
  [[ -s $FA ]] || die "missing assembly ${FA}"
fi
OUT=${D19}/${NAME}.tsv.gz
[[ -s $OUT ]] && { echo "${NAME}: done already"; exit 0; }
T=${KMER_DIR}/tmp/ref_${NAME}; rm -rf "$T"; mkdir -p "$T"

if [[ $NAME == B73 ]]; then
  kmc_tools transform B73_uniq dump "${T}/u.txt"                      # all B73 single-copy k-mers
else
  kmc -k"$KMER_K" -m40 -t"$CPUS" -ci1 -fm "$FA" "${T}/asm" "$T" > "${T}/kmc_asm.log" 2>&1
  kmc_tools simple "${T}/asm" B73_uniq intersect "${T}/u" -ocleft     # B73 single-copy k-mers in the assembly
  kmc_tools transform "${T}/u" dump "${T}/u.txt"
  rm -f "${T}"/asm.kmc_* "${T}"/u.kmc_*
fi
awk '{print ">k\n" $1}' "${T}/u.txt" > "${T}/u.fa"; rm -f "${T}/u.txt"
kmc -k"$KMER_K" -m40 -t"$CPUS" -ci1 -fm "${T}/u.fa" "${T}/u2" "$T" > "${T}/kmc_u2.log" 2>&1      # KMC2 format for kcftools
rm -f "${T}/u.fa"
kcftools getVariations -r B73_chr1-10.fa -k "${T}/u2" -o "${T}/${NAME}.kcf" -s "$NAME" -f window -w "$KMER_WINDOW" -t "$CPUS" > "${T}/kcftools.log" 2>&1
[[ -s ${T}/${NAME}.kcf ]] || { cat "${T}/kcftools.log"; die "kcftools wrote nothing for ${NAME}"; }
kcftools kcf2tsv -i "${T}/${NAME}.kcf" -o "${T}/${NAME}" -s "$NAME" > /dev/null 2>&1
[[ -s ${T}/${NAME}.${NAME}.tsv ]] || die "kcf2tsv wrote nothing for ${NAME}"
gzip -c "${T}/${NAME}.${NAME}.tsv" > "${OUT}.tmp" && mv "${OUT}.tmp" "$OUT"
rm -rf "$T"
echo "${NAME}: single-copy B73 k-mers in windows: $(zcat "$OUT" | awk -F'\t' 'NR > 1 {s += $7} END {printf "%.0f M", s / 1e6}')"
echo "WROTE: ${OUT}"
