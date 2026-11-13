#!/bin/bash
#SBATCH --job-name=hmm_family
#SBATCH --array=1-143%50
#SBATCH --output=logs/14b_hmm_family.%A_%a.out
#SBATCH --error=logs/14b_hmm_family.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mem=8G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# 14b — Family-informed HMM run (same 14_hmm_pooled.py, plus the BC1 family's donor
#       alleles). Run alongside step 14 to compare.
#
#   At each marker, c (chance the donor carries ALT) = the family's PALT from step 13c
#   where known, else teosinte panel AC/AN x kappa (as in step 14). Non-panel PASS sites
#   are used where the family's PALT is known. Pools whose family has no donor column
#   (or listed in sample_lists/no_family.txt) run exactly like step 14.
#   Inputs: ${D13B} (GLs at all PASS sites), ${D13C} (family PALT), pedigree in SAMPLE_META
#   Output per pool in ${D14B}: <SAMPLE>.statepath.tsv.gz (+ C, C_SOURCE), .tracts.bed.gz, .summary.tsv
#   Submit after 13b and 13c:  sbatch 14b_run_hmm_family.sh
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"
need python3

IDX=${SLURM_ARRAY_TASK_ID:-${1:?usage: bash $0 <task index>}}
[[ -s ${D11}/pools.tsv ]] || die "missing ${D11}/pools.tsv (run 11)"
mkdir -p "$D14B"

mapfile -t ALL < <(awk -F'\t' 'NR > 1 && ($2 == "control" || $2 == "donor") {print $1}' "${D11}/pools.tsv")
NT=$(( (${#ALL[@]} + SAMPLES_PER_TASK - 1) / SAMPLES_PER_TASK ))
(( IDX <= NT )) || { echo "task ${IDX} > ${NT} tasks needed for ${#ALL[@]} pools; nothing to do"; exit 0; }
FIRST=$(( (IDX - 1) * SAMPLES_PER_TASK ))

# pool -> BC1 family (first two pedigree fields), from the metadata (CRLF line endings)
declare -A FAMILY
while IFS=$'\t' read -r sid geno; do
  FAMILY[$sid]=$(echo "$geno" | cut -d_ -f1,2)
done < <(tr -d '\r' < "$SAMPLE_META" | awk -F'\t' 'NR > 1 {print $1 "\t" $2}')
# pools to run without family info (e.g. flagged CHECK by 15b): sample_lists/no_family.txt
if [[ -s ${SAMPLE_LISTS}/no_family.txt ]]; then
  while read -r sid; do [[ -n $sid ]] && FAMILY[$sid]=none; done < <(tr -d '\r' < "${SAMPLE_LISTS}/no_family.txt")
fi

CHR_LIST=$(IFS=,; echo "${CHROMS[*]}")
echo "task ${IDX}/${NT}: pools $((FIRST + 1))-$(( FIRST + SAMPLES_PER_TASK < ${#ALL[@]} ? FIRST + SAMPLES_PER_TASK : ${#ALL[@]} )) of ${#ALL[@]} | rho ${HMM_RHO} | ${HMM_DECODE} R=${HMM_RIGIDITY}"

for S in "${ALL[@]:FIRST:SAMPLES_PER_TASK}"; do
  if [[ -s ${D14B}/${S}.summary.tsv ]]; then echo "${S}: done already"; continue; fi
  FAM=${FAMILY[$S]:-none}
  python3 "${FINAL_DIR}/14_hmm_pooled.py" --sample "$S" --in_dir "$D13B" --map_dir "$MAP_DIR" \
    --fai "${REF}.fai" --out_prefix "${D14B}/${S}" --chroms "$CHR_LIST" \
    --donor_dir "$D13C" --family "$FAM" \
    --rho "$HMM_RHO" --decode "$HMM_DECODE" --rigidity "$HMM_RIGIDITY"
  echo "${S} (${FAM}): $(tail -n 1 "${D14B}/${S}.summary.tsv" | cut -f2,6,7,9,13)"
done
echo "DONE"
