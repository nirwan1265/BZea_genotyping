#!/bin/bash
# =============================================================================
# 18 — Check step 17 (run on the login node: bash 18_kmer_check.sh).
#   Lists pools without output, MISSING FASTQ / FAILED lines from the logs, and writes
#   ${KMER_DIR}/kmer_pool_qc.tsv: SAMPLE ROLE GROUP WINDOWS OBSERVED_KMERS_M
#   (pools with few observed k-mers are low depth and give noisy calls).
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"

QC=${KMER_DIR}/kmer_pool_qc.tsv
printf 'SAMPLE\tROLE\tGROUP\tWINDOWS\tOBSERVED_KMERS_M\n' > "$QC"
n=0; miss=0
while IFS=$'\t' read -r S role group _; do
  n=$(( n + 1 ))
  f=${D16}/${S}.tsv.gz
  if [[ ! -s $f ]]; then miss=$(( miss + 1 )); echo "no output: ${S}"; continue; fi
  zcat "$f" | awk -F'\t' -v s="$S" -v r="$role" -v g="$group" 'NR > 1 {w++; o += $7} END {printf "%s\t%s\t%s\t%d\t%.1f\n", s, r, g, w, o / 1e6}' >> "$QC"
done < <(awk -F'\t' 'NR > 1 && ($2 == "control" || $2 == "donor")' "${D11}/pools.tsv")
echo "pools: ${n} | with output: $(( n - miss )) | missing: ${miss}"
grep -h -E "MISSING FASTQ|FAILED" logs/17_kmer_win.*.out 2>/dev/null | sort -u | head -50 || true
echo "observed single-copy k-mers per pool (M): $(awk -F'\t' 'NR > 1 {print $5}' "$QC" | sort -g | awk '{a[NR] = $1} END {if (NR) printf "min %s | median %s | max %s", a[1], a[int((NR + 1) / 2)], a[NR]}')"
echo "WROTE: ${QC}"
