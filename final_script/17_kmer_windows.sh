#!/bin/bash
#SBATCH --job-name=kmer_win
#SBATCH --output=logs/17_kmer_win.%A_%a.out
#SBATCH --error=logs/17_kmer_win.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=10:00:00
#SBATCH --mem=48G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# 17 — k-mer windows for every pool (Slurm array, KMER_POOLS_PER_TASK pools per task).
#
#   Pools: ${D11}/pools.tsv, ROLE control or donor (same pools as steps 13-15).
#   Per pool:
#     1. count 31-mers in its trimmed paired reads (KMC, -ci1: at 0.8x most k-mers are seen once)
#        (reuses an existing database in ${KMER_DIR}/kmc/<SAMPLE> or ${KMER_DIR}/<SAMPLE>)
#     2. keep k-mers that are single-copy in B73 (kmc_tools intersect with B73_uniq)
#     3. KCFtools reads only KMC2-format databases and kmc_tools writes KMC1, so the kept
#        k-mers are dumped and re-counted with kmc (same k-mers, KMC2 format)
#     4. kcftools getVariations on B73 chr1-10, ${KMER_WINDOW}-bp windows
#   Output: ${D16}/<SAMPLE>.tsv.gz (kcftools kcf2tsv table; observed_kmers = B73 single-copy
#           k-mers of the window seen in the pool). Pools already done are skipped; outputs of
#           the 44-pool test batch (${KMER_DIR}/batch) are reused.
#   Temporary disk: ~25 GB per running pool in ${KMER_DIR}/tmp (deleted after each pool).
#
#   Submit after 16, with the k-mer env active:
#     N=$(awk -F'\t' 'NR>1 && ($2=="control"||$2=="donor")' <D11>/pools.tsv | wc -l)
#     sbatch --array=1-$(( (N + 3) / 4 ))%10 17_kmer_windows.sh
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"
need kmc kmc_tools kcftools awk gzip

IDX=${SLURM_ARRAY_TASK_ID:-${1:?usage: bash $0 <task index>}}
CPUS=${SLURM_CPUS_PER_TASK:-8}
cd "$KMER_DIR"
[[ -s B73_uniq.kmc_pre && -s B73_chr1-10.fa ]] || die "missing B73_uniq / B73_chr1-10.fa in ${KMER_DIR} (run 16)"
[[ -s ${D11}/pools.tsv ]] || die "missing ${D11}/pools.tsv (run 11)"
mkdir -p "$D16" "${KMER_DIR}/tmp" "${KMER_DIR}/kmc"

mapfile -t POOLS < <(awk -F'\t' 'NR > 1 && ($2 == "control" || $2 == "donor") {print $1}' "${D11}/pools.tsv")
FIRST=$(( (IDX - 1) * KMER_POOLS_PER_TASK ))
(( FIRST < ${#POOLS[@]} )) || { echo "task ${IDX}: no pools (only ${#POOLS[@]})"; exit 0; }

# every command checked with "|| return 1": set -e does not act inside a function called with ||
one_pool() {
  local S=$1 out=${D16}/${1}.tsv.gz T=${KMER_DIR}/tmp/${1} reads plate R1 R2
  [[ -s $out ]] && { echo "${S}: done already"; return 0; }
  if [[ -s ${KMER_DIR}/batch/${S}.uniq.scores.${S}.tsv ]]; then          # from the test batch
    gzip -c "${KMER_DIR}/batch/${S}.uniq.scores.${S}.tsv" > "$out" || return 1
    echo "${S}: copied from test batch"; return 0
  fi
  rm -rf "$T"; mkdir -p "$T" || return 1
  if   [[ -s ${KMER_DIR}/kmc/${S}.kmc_pre ]]; then reads=${KMER_DIR}/kmc/${S}
  elif [[ -s ${KMER_DIR}/${S}.kmc_pre ]];     then reads=${KMER_DIR}/${S}
  else
    plate=${S%%_*}; plate=${plate#PN}
    R1=${FASTQ_DIR}/filtered_S${plate}/${S}/${S}_paired_R1.fq.gz; R2=${R1/_R1/_R2}
    if [[ ! -s $R1 || ! -s $R2 ]]; then echo "${S}: MISSING FASTQ ${R1}"; rm -rf "$T"; return 0; fi
    printf '%s\n%s\n' "$R1" "$R2" > "${T}/fqlist"
    kmc -k"$KMER_K" -m40 -t"$CPUS" -ci1 -p9 -fq @"${T}/fqlist" "${T}/reads" "$T" > "${T}/kmc_reads.log" 2>&1 || return 1
    reads=${T}/reads
  fi
  kmc_tools simple "$reads" B73_uniq intersect "${T}/u" -ocleft || return 1
  kmc_tools transform "${T}/u" dump "${T}/u.txt" || return 1
  awk '{print ">k\n" $1}' "${T}/u.txt" > "${T}/u.fa" || return 1
  rm -f "${T}/u.txt" "${T}"/u.kmc_*
  kmc -k"$KMER_K" -m40 -t"$CPUS" -ci1 -fm "${T}/u.fa" "${T}/u2" "$T" > "${T}/kmc_u2.log" 2>&1 || return 1
  rm -f "${T}/u.fa"
  kcftools getVariations -r B73_chr1-10.fa -k "${T}/u2" -o "${T}/${S}.kcf" -s "$S" -f window -w "$KMER_WINDOW" -t "$CPUS" > "${T}/kcftools.log" 2>&1 || return 1
  [[ -s ${T}/${S}.kcf ]] || { echo "${S}: kcftools wrote nothing"; cat "${T}/kcftools.log"; return 1; }
  kcftools kcf2tsv -i "${T}/${S}.kcf" -o "${T}/${S}" -s "$S" > /dev/null 2>&1 || return 1
  [[ -s ${T}/${S}.${S}.tsv ]] || { echo "${S}: kcf2tsv wrote nothing"; return 1; }
  gzip -c "${T}/${S}.${S}.tsv" > "${out}.tmp" && mv "${out}.tmp" "$out" || return 1
  if [[ $KMER_KEEP_KMC == 1 && $reads == "${T}/reads" ]]; then mv "${T}"/reads.kmc_* "${KMER_DIR}/kmc/"; fi
  rm -rf "$T"
  echo "${S}: done ($(( $(zcat "$out" | wc -l) - 1 )) windows)"
}

for (( i = FIRST; i < FIRST + KMER_POOLS_PER_TASK && i < ${#POOLS[@]}; i++ )); do
  one_pool "${POOLS[$i]}" || echo "${POOLS[$i]}: FAILED (see ${KMER_DIR}/tmp/${POOLS[$i]})"
done
