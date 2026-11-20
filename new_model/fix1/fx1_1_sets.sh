#!/bin/bash
#SBATCH --job-name=fx1_sets
#SBATCH --output=logs/fx1_1_sets.%A_%a.out
#SBATCH --error=logs/fx1_1_sets.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=10:00:00
#SBATCH --mem=48G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# fix1 / 1 — informative B73 k-mers per taxon, whole genome (array 1-5: Zd Zh Zl Zv Zx).
#   B73_inf_<T> = B73 single-copy k-mers (B73_uniq, step 16) in NONE of the taxon's assemblies (RULE=all of the
#   pilot). Zv and Zx are copied from the pilot (same rule, the k-mer sets were already whole-genome).
#   Then the per-window count of these k-mers in B73 on chr1-10 (the denominator, like step 19's B73.tsv.gz).
#   Output: ${FX}/sets/B73_inf_<T>.{kmc_pre,kmc_suf,tsv.gz}
# =============================================================================
set -euo pipefail
HERE=$(cd "${SLURM_SUBMIT_DIR:-$PWD}" && pwd)
source "${HERE}/../../final_script/config.sh"; source "${HERE}/fx1_config.sh"; use_kmer_env
need kmc kmc_tools kcftools awk gzip
CPUS=${SLURM_CPUS_PER_TASK:-8}
T=${FX_TAXA[$(( ${SLURM_ARRAY_TASK_ID:-1} - 1 ))]}
SET=${FX}/sets/B73_inf_${T}
mkdir -p "${FX}/sets" "${FX}/tmp"
cd "$KMER_DIR"
[[ -s B73_uniq.kmc_pre && -s B73_chr1-10.fa ]] || die "missing B73_uniq / B73_chr1-10.fa in ${KMER_DIR} (step 16)"
[[ -s ${SET}.tsv.gz ]] && { echo "${T}: done already"; exit 0; }
TMP=${FX}/tmp/set_${T}; rm -rf "$TMP"; mkdir -p "$TMP"
if [[ ! -s ${SET}.kmc_pre ]]; then
  if [[ -s ${PILOT_SETS}/B73_inf_${T}.kmc_pre ]]; then
    cp "${PILOT_SETS}/B73_inf_${T}.kmc_pre" "${SET}.kmc_pre.tmp"; cp "${PILOT_SETS}/B73_inf_${T}.kmc_suf" "${SET}.kmc_suf"
    mv "${SET}.kmc_pre.tmp" "${SET}.kmc_pre"; echo "${T}: k-mer set copied from the pilot"
  else
    cur=${KMER_DIR}/B73_uniq; i=0
    for a in "${DONOR_ASSEMBLIES[@]}"; do
      IFS=: read -r NAME TAX FA <<< "$a"; [[ $TAX == "$T" ]] || continue
      i=$((i + 1)); echo "assembly ${NAME}"
      kmc -k"$KMER_K" -m40 -t"$CPUS" -ci1 -fm "$FA" "${TMP}/asm${i}" "$TMP" > "${TMP}/kmc_asm_${i}.log" 2>&1
      kmc_tools simple "$cur" "${TMP}/asm${i}" kmers_subtract "${TMP}/s${i}"
      rm -f "${TMP}"/asm${i}.kmc_*; cur=${TMP}/s${i}
    done
    (( i > 0 )) || die "no assembly for ${T} in DONOR_ASSEMBLIES"
    mv "${cur}.kmc_suf" "${SET}.kmc_suf"; mv "${cur}.kmc_pre" "${SET}.kmc_pre"
  fi
fi
kmc_tools transform "$SET" dump "${TMP}/u.txt"
echo "${T}: $(wc -l < "${TMP}/u.txt") informative k-mers"
awk '{print ">k\n" $1}' "${TMP}/u.txt" > "${TMP}/u.fa"; rm -f "${TMP}/u.txt"
kmc -k"$KMER_K" -m40 -t"$CPUS" -ci1 -fm "${TMP}/u.fa" "${TMP}/u2" "$TMP" > "${TMP}/kmc_u2.log" 2>&1; rm -f "${TMP}/u.fa"
kcftools getVariations -r "${KMER_DIR}/B73_chr1-10.fa" -k "${TMP}/u2" -o "${TMP}/B73.kcf" -s B73 -f window -w "$KMER_WINDOW" -t "$CPUS" > "${TMP}/kcftools.log" 2>&1
kcftools kcf2tsv -i "${TMP}/B73.kcf" -o "${TMP}/B73" -s B73 > /dev/null 2>&1
[[ -s ${TMP}/B73.B73.tsv ]] || die "kcftools wrote nothing (see ${TMP}/kcftools.log)"
gzip -c "${TMP}/B73.B73.tsv" > "${SET}.tsv.gz.tmp" && mv "${SET}.tsv.gz.tmp" "${SET}.tsv.gz"
rm -rf "$TMP"
echo "${T}: wrote ${SET}.tsv.gz ($(( $(zcat "${SET}.tsv.gz" | wc -l) - 1 )) windows)"
