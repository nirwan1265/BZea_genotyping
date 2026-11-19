#!/bin/bash
#SBATCH --job-name=p1_pools
#SBATCH --output=logs/p1_pool_windows.%A_%a.out
#SBATCH --error=logs/p1_pool_windows.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=06:00:00
#SBATCH --mem=48G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# pilot fix 1 / step 2 — informative-k-mer windows on chr3 + chr7 for the pilot pools and B73 controls.
#   One array task per pool: tasks 1-60 = pools.tsv (this folder), 61-70 = the first 10 B73 controls of
#   <D11>/pools.tsv. Per pool: count its reads with KMC (as step 17; reuses ${KMER_DIR}/kmc/<SAMPLE> if
#   it exists), intersect with B73_inf_<TAXON> (controls: both Zv and Zx), windows with kcftools on chr3 + chr7.
#   Output: ${P}/windows/<TAXON>/<SAMPLE>.tsv.gz. Temporary disk ~25 GB per running task (deleted after).
#   Run from new_model/pilot_fix1/ after step 1, k-mer env active:  sbatch --array=1-70%100 p1_pool_windows.sh
# =============================================================================
set -euo pipefail
HERE=$(cd "${SLURM_SUBMIT_DIR:-$PWD}" && pwd)
source "${HERE}/../../final_script/config.sh"
need kmc kmc_tools kcftools awk gzip
CPUS=${SLURM_CPUS_PER_TASK:-8}
P=${KMER_DIR}/new_model/pilot_fix1
IDX=${SLURM_ARRAY_TASK_ID:-${1:?usage: bash $0 <task index>}}
mapfile -t PILOT < <(awk -F'\t' 'NR > 1 {print $1 "\t" $2}' "${HERE}/pools.tsv")
mapfile -t CTRL < <(awk -F'\t' 'NR > 1 && $2 == "control" {print $1}' "${D11}/pools.tsv" | head -10)
if (( IDX <= ${#PILOT[@]} )); then
  IFS=$'\t' read -r S TX <<< "${PILOT[$((IDX - 1))]}"; SETS=("$TX")
else
  k=$(( IDX - ${#PILOT[@]} - 1 )); (( k < ${#CTRL[@]} )) || { echo "task ${IDX}: nothing"; exit 0; }
  S=${CTRL[$k]}; SETS=(Zv Zx)
fi
todo=0; for T in "${SETS[@]}"; do [[ -s ${P}/windows/${T}/${S}.tsv.gz ]] || todo=1; done
(( todo )) || { echo "${S}: done already"; exit 0; }
for T in "${SETS[@]}"; do [[ -s ${P}/sets/B73_inf_${T}.kmc_pre ]] || die "missing ${P}/sets/B73_inf_${T} (run p1_build_sets.sh)"; done
TMP=${P}/tmp/pool_${S}; rm -rf "$TMP"; mkdir -p "$TMP"
if   [[ -s ${KMER_DIR}/kmc/${S}.kmc_pre ]]; then reads=${KMER_DIR}/kmc/${S}
elif [[ -s ${KMER_DIR}/${S}.kmc_pre ]];     then reads=${KMER_DIR}/${S}
else
  plate=${S%%_*}; plate=${plate#PN}
  R1=${FASTQ_DIR}/filtered_S${plate}/${S}/${S}_paired_R1.fq.gz; R2=${R1/_R1/_R2}
  [[ -s $R1 && -s $R2 ]] || die "${S}: missing FASTQ ${R1}"
  printf '%s\n%s\n' "$R1" "$R2" > "${TMP}/fqlist"
  kmc -k"$KMER_K" -m40 -t"$CPUS" -ci1 -p9 -fq @"${TMP}/fqlist" "${TMP}/reads" "$TMP" > "${TMP}/kmc_reads.log" 2>&1
  reads=${TMP}/reads
fi
for T in "${SETS[@]}"; do
  mkdir -p "${P}/windows/${T}"; out=${P}/windows/${T}/${S}.tsv.gz
  [[ -s $out ]] && continue
  kmc_tools simple "$reads" "${P}/sets/B73_inf_${T}" intersect "${TMP}/u" -ocleft
  kmc_tools transform "${TMP}/u" dump "${TMP}/u.txt"
  awk '{print ">k\n" $1}' "${TMP}/u.txt" > "${TMP}/u.fa"; rm -f "${TMP}/u.txt" "${TMP}"/u.kmc_*
  kmc -k"$KMER_K" -m40 -t"$CPUS" -ci1 -fm "${TMP}/u.fa" "${TMP}/u2" "$TMP" > "${TMP}/kmc_u2.log" 2>&1; rm -f "${TMP}/u.fa"
  kcftools getVariations -r "${P}/B73_chr3_7.fa" -k "${TMP}/u2" -o "${TMP}/${S}.kcf" -s "$S" -f window -w "$KMER_WINDOW" -t "$CPUS" > "${TMP}/kcftools.log" 2>&1
  kcftools kcf2tsv -i "${TMP}/${S}.kcf" -o "${TMP}/${S}" -s "$S" > /dev/null 2>&1
  gzip -c "${TMP}/${S}.${S}.tsv" > "${out}.tmp" && mv "${out}.tmp" "$out"
  rm -f "${TMP}"/u2.kmc_* "${TMP}/${S}".*
  echo "${S} (${T}): $(( $(zcat "$out" | wc -l) - 1 )) windows"
done
rm -rf "$TMP"
