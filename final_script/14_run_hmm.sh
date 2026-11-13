#!/bin/bash
#SBATCH --job-name=hmm
#SBATCH --array=1-143%50
#SBATCH --output=logs/14_hmm.%A_%a.out
#SBATCH --error=logs/14_hmm.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mem=8G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# 14 — Run the HMM (14_hmm_pooled.py, adapted from hmm_viterbi_bc2s3.py) per pool.
#      Each array task runs SAMPLES_PER_TASK pools (default 10 -> 1,424 pools = 143 tasks).
#
#   Pools: B73 controls + donor pools from ${D11}/pools.tsv (excluded pools left out).
#   Input: per-pool GL files from step 13 (sequencing error is inside the GLs).
#
#   Output per pool in ${D14}: <SAMPLE>.statepath.tsv.gz, <SAMPLE>.tracts.bed.gz, <SAMPLE>.summary.tsv
#   Submit after 13:  sbatch 14_run_hmm.sh
#   (if the number of pools changes: sbatch --array=1-<ceil(pools/10)>%50 14_run_hmm.sh)
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"
need python3

IDX=${SLURM_ARRAY_TASK_ID:-${1:?usage: bash $0 <task index>}}
[[ -s ${D11}/pools.tsv ]] || die "missing ${D11}/pools.tsv (run 11)"
mkdir -p "$D14"

mapfile -t ALL < <(awk -F'\t' 'NR > 1 && ($2 == "control" || $2 == "donor") {print $1}' "${D11}/pools.tsv")
NT=$(( (${#ALL[@]} + SAMPLES_PER_TASK - 1) / SAMPLES_PER_TASK ))
(( IDX <= NT )) || { echo "task ${IDX} > ${NT} tasks needed for ${#ALL[@]} pools; nothing to do"; exit 0; }
FIRST=$(( (IDX - 1) * SAMPLES_PER_TASK ))

CHR_LIST=$(IFS=,; echo "${CHROMS[*]}")
echo "task ${IDX}/${NT}: pools $((FIRST + 1))-$(( FIRST + SAMPLES_PER_TASK < ${#ALL[@]} ? FIRST + SAMPLES_PER_TASK : ${#ALL[@]} )) of ${#ALL[@]} | rho ${HMM_RHO} | ${HMM_DECODE} R=${HMM_RIGIDITY}"

for S in "${ALL[@]:FIRST:SAMPLES_PER_TASK}"; do
  if [[ -s ${D14}/${S}.summary.tsv ]]; then echo "${S}: done already"; continue; fi
  python3 "${FINAL_DIR}/14_hmm_pooled.py" --sample "$S" --in_dir "$D13" --map_dir "$MAP_DIR" \
    --fai "${REF}.fai" --out_prefix "${D14}/${S}" --chroms "$CHR_LIST" \
    --rho "$HMM_RHO" --decode "$HMM_DECODE" --rigidity "$HMM_RIGIDITY"
  echo "${S}: $(tail -n 1 "${D14}/${S}.summary.tsv" | cut -f2,6,7)"
done
echo "DONE"
