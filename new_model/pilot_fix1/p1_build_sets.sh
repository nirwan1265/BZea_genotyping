#!/bin/bash
#SBATCH --job-name=p1_sets
#SBATCH --output=logs/p1_build_sets.%A_%a.out
#SBATCH --error=logs/p1_build_sets.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=08:00:00
#SBATCH --mem=48G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# pilot fix 1 / step 1 — informative B73 k-mers per taxon (one array task per taxon: 1 = Zv, 2 = Zx).
#   B73_inf_<T> = B73 single-copy k-mers (B73_uniq, step 16) that are in NONE of the taxon's assemblies
#   (kmc_tools kmers_subtract). Only these can tell DONOR from SEG; k-mers the teosinte shares with B73 are
#   present in every pool whatever its state and only add counting noise.
#   Also: B73_chr3_7.fa (chr3 + chr7 of B73_chr1-10.fa; the pilot is run on these two chromosomes) and the
#   per-window count of informative k-mers in B73 (B73_inf_<T>.tsv.gz, the denominator, like step 19's B73.tsv.gz).
#   RULE=all (default): k-mers in NONE of the taxon's assemblies -> sets/B73_inf_<T>
#   RULE=any: k-mers missing from AT LEAST ONE assembly (B73_uniq minus the k-mers shared by all of them)
#            -> sets/B73_inf_<T>_any. Use when one assembly looks like B73 in places (Zv chr3 75-79 Mb kept
#            only 34-505 k-mers per window under RULE=all): sbatch --export=ALL,RULE=any --array=1 p1_build_sets.sh
#   Output: ${P}/sets/ with P = ${KMER_DIR}/new_model/pilot_fix1
#   Run from new_model/pilot_fix1/ with the k-mer env active (kmc, kmc_tools, kcftools):
#     mkdir -p logs && sbatch --array=1-2 p1_build_sets.sh
# =============================================================================
set -euo pipefail
HERE=$(cd "${SLURM_SUBMIT_DIR:-$PWD}" && pwd)
source "${HERE}/../../final_script/config.sh"
need kmc kmc_tools kcftools awk gzip
CPUS=${SLURM_CPUS_PER_TASK:-8}
P=${KMER_DIR}/new_model/pilot_fix1; mkdir -p "${P}/sets" "${P}/tmp"
TAXA=(Zv Zx); T=${TAXA[$(( ${SLURM_ARRAY_TASK_ID:-1} - 1 ))]}
RULE=${RULE:-all}; [[ $RULE == all || $RULE == any ]] || die "RULE must be all or any"
TAG=""; [[ $RULE == any ]] && TAG=_any
SET=B73_inf_${T}${TAG}
cd "$KMER_DIR"
[[ -s B73_uniq.kmc_pre && -s B73_chr1-10.fa ]] || die "missing B73_uniq / B73_chr1-10.fa in ${KMER_DIR} (step 16)"
if [[ ! -s ${P}/B73_chr3_7.fa ]]; then
  awk '/^>/{p=($1==">chr3"||$1==">chr7")} p' B73_chr1-10.fa > "${P}/B73_chr3_7.fa.${T}.tmp"
  [[ -s ${P}/B73_chr3_7.fa ]] || mv "${P}/B73_chr3_7.fa.${T}.tmp" "${P}/B73_chr3_7.fa"
  rm -f "${P}/B73_chr3_7.fa.${T}.tmp"
fi
grep '^>' "${P}/B73_chr3_7.fa"
TMP=${P}/tmp/set_${T}${TAG}; rm -rf "$TMP"; mkdir -p "$TMP"
cur=${KMER_DIR}/B73_uniq; shared=""; i=0
for a in "${DONOR_ASSEMBLIES[@]}"; do
  IFS=: read -r NAME TAX FA <<< "$a"
  [[ $TAX == "$T" ]] || continue
  i=$((i + 1)); echo "${RULE}: assembly ${NAME}"
  kmc -k"$KMER_K" -m40 -t"$CPUS" -ci1 -fm "$FA" "${TMP}/asm${i}" "$TMP" > "${TMP}/kmc_asm_${i}.log" 2>&1
  if [[ $RULE == all ]]; then                               # remove every k-mer of every assembly
    kmc_tools simple "$cur" "${TMP}/asm${i}" kmers_subtract "${TMP}/s${i}"
    rm -f "${TMP}"/asm${i}.kmc_*; cur=${TMP}/s${i}
  else                                                       # keep the k-mers shared by all assemblies
    if [[ -z $shared ]]; then shared=${TMP}/asm${i}
    else kmc_tools simple "$shared" "${TMP}/asm${i}" intersect "${TMP}/sh${i}"; shared=${TMP}/sh${i}; fi
  fi
done
(( i > 0 )) || die "no assembly for ${T} in DONOR_ASSEMBLIES"
if [[ $RULE == any ]]; then kmc_tools simple "$cur" "$shared" kmers_subtract "${TMP}/s_any"; cur=${TMP}/s_any; fi
mv "${cur}.kmc_pre" "${P}/sets/${SET}.kmc_pre"; mv "${cur}.kmc_suf" "${P}/sets/${SET}.kmc_suf"
rm -f "${TMP}"/asm*.kmc_* "${TMP}"/sh*.kmc_*
kmc_tools transform "${P}/sets/${SET}" dump "${TMP}/u.txt"
echo "${T} (${RULE}): $(wc -l < "${TMP}/u.txt") informative k-mers (of $(kmc_tools info "${KMER_DIR}/B73_uniq" 2>/dev/null | awk '/total k-mers/{print $NF}') B73 single-copy)"
awk '{print ">k\n" $1}' "${TMP}/u.txt" > "${TMP}/u.fa"; rm -f "${TMP}/u.txt"
kmc -k"$KMER_K" -m40 -t"$CPUS" -ci1 -fm "${TMP}/u.fa" "${TMP}/u2" "$TMP" > "${TMP}/kmc_u2.log" 2>&1; rm -f "${TMP}/u.fa"
kcftools getVariations -r "${P}/B73_chr3_7.fa" -k "${TMP}/u2" -o "${TMP}/B73.kcf" -s B73 -f window -w "$KMER_WINDOW" -t "$CPUS" > "${TMP}/kcftools.log" 2>&1
kcftools kcf2tsv -i "${TMP}/B73.kcf" -o "${TMP}/B73" -s B73 > /dev/null 2>&1
gzip -c "${TMP}/B73.B73.tsv" > "${P}/sets/${SET}.tsv.gz"
rm -rf "$TMP"
echo "wrote ${P}/sets/${SET}.{kmc_pre,kmc_suf,tsv.gz}"
