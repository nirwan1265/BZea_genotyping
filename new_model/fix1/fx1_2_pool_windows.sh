#!/bin/bash
#SBATCH --job-name=fx1_pools
#SBATCH --output=logs/fx1_2_pools.%A_%a.out
#SBATCH --error=logs/fx1_2_pools.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --mem=48G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# fix1 / 2 — informative-k-mer windows for every pool, whole genome (FX_POOLS_PER_TASK pools per task).
#   Pools: the B73 controls (all 5 taxa's sets) then the donor pools (their taxon's set), see fx_pools in fx1_config.sh.
#   Per pool, as step 17: count its reads with KMC, keep the informative k-mers (kmc_tools intersect with
#   ${FX}/sets/B73_inf_<T>), re-count in KMC2 format, kcftools windows on B73 chr1-10.
#   Output: ${FX}/windows/<T>/<SAMPLE>.tsv.gz. Pools already done are skipped (safe to resubmit).
#   Temporary disk ~25 GB per running pool (${FX}/tmp, deleted after each pool).
# =============================================================================
set -euo pipefail
HERE=$(cd "${SLURM_SUBMIT_DIR:-$PWD}" && pwd)
source "${HERE}/../../final_script/config.sh"; source "${HERE}/fx1_config.sh"
need kmc kmc_tools kcftools awk gzip
CPUS=${SLURM_CPUS_PER_TASK:-8}
IDX=${SLURM_ARRAY_TASK_ID:-${1:?usage: bash $0 <task index>}}
for T in "${FX_TAXA[@]}"; do [[ -s ${FX}/sets/B73_inf_${T}.kmc_pre ]] || die "missing ${FX}/sets/B73_inf_${T} (fx1_1)"; done
mapfile -t POOLS < <(fx_pools)
FIRST=$(( (IDX - 1) * FX_POOLS_PER_TASK ))
(( FIRST < ${#POOLS[@]} )) || { echo "task ${IDX}: no pools (only ${#POOLS[@]})"; exit 0; }

one_pool() {
  local S=$1 G=$2 TMP=${FX}/tmp/pool_${1} reads plate R1 R2 T out todo=0 sets
  if [[ $G == control ]]; then sets=("${FX_TAXA[@]}"); else sets=("$G"); fi
  for T in "${sets[@]}"; do [[ -s ${FX}/windows/${T}/${S}.tsv.gz ]] || todo=1; done
  (( todo )) || { echo "${S}: done already"; return 0; }
  rm -rf "$TMP"; mkdir -p "$TMP" || return 1
  if   [[ -s ${KMER_DIR}/kmc/${S}.kmc_pre ]]; then reads=${KMER_DIR}/kmc/${S}
  elif [[ -s ${KMER_DIR}/${S}.kmc_pre ]];     then reads=${KMER_DIR}/${S}
  else
    plate=${S%%_*}; plate=${plate#PN}
    R1=${FASTQ_DIR}/filtered_S${plate}/${S}/${S}_paired_R1.fq.gz; R2=${R1/_R1/_R2}
    if [[ ! -s $R1 || ! -s $R2 ]]; then echo "${S}: MISSING FASTQ ${R1}"; rm -rf "$TMP"; return 0; fi
    printf '%s\n%s\n' "$R1" "$R2" > "${TMP}/fqlist"
    kmc -k"$KMER_K" -m40 -t"$CPUS" -ci1 -p9 -fq @"${TMP}/fqlist" "${TMP}/reads" "$TMP" > "${TMP}/kmc_reads.log" 2>&1 || return 1
    reads=${TMP}/reads
  fi
  for T in "${sets[@]}"; do
    out=${FX}/windows/${T}/${S}.tsv.gz; mkdir -p "${FX}/windows/${T}"
    [[ -s $out ]] && continue
    kmc_tools simple "$reads" "${FX}/sets/B73_inf_${T}" intersect "${TMP}/u" -ocleft || return 1
    kmc_tools transform "${TMP}/u" dump "${TMP}/u.txt" || return 1
    awk '{print ">k\n" $1}' "${TMP}/u.txt" > "${TMP}/u.fa" || return 1
    rm -f "${TMP}/u.txt" "${TMP}"/u.kmc_*
    kmc -k"$KMER_K" -m40 -t"$CPUS" -ci1 -fm "${TMP}/u.fa" "${TMP}/u2" "$TMP" > "${TMP}/kmc_u2.log" 2>&1 || return 1
    rm -f "${TMP}/u.fa"
    kcftools getVariations -r "${KMER_DIR}/B73_chr1-10.fa" -k "${TMP}/u2" -o "${TMP}/${S}.kcf" -s "$S" -f window -w "$KMER_WINDOW" -t "$CPUS" > "${TMP}/kcftools.log" 2>&1 || return 1
    kcftools kcf2tsv -i "${TMP}/${S}.kcf" -o "${TMP}/${S}" -s "$S" > /dev/null 2>&1 || return 1
    [[ -s ${TMP}/${S}.${S}.tsv ]] || { echo "${S} (${T}): kcf2tsv wrote nothing"; return 1; }
    gzip -c "${TMP}/${S}.${S}.tsv" > "${out}.tmp" && mv "${out}.tmp" "$out" || return 1
    rm -f "${TMP}"/u2.kmc_* "${TMP}/${S}".*
    echo "${S} (${T}): $(( $(zcat "$out" | wc -l) - 1 )) windows"
  done
  rm -rf "$TMP"
}

for (( i = FIRST; i < FIRST + FX_POOLS_PER_TASK && i < ${#POOLS[@]}; i++ )); do
  IFS=$'\t' read -r S G <<< "${POOLS[$i]}"
  one_pool "$S" "$G" || echo "${S}: FAILED (see ${FX}/tmp/pool_${S})"
done
