#!/bin/bash
#SBATCH --job-name=kmer_hmm_new
#SBATCH --output=logs/20n_kmer_hmm.%A_%a.out
#SBATCH --error=logs/20n_kmer_hmm.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH --mem=8G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# new_model / 20n — the step-20 HMM with DONOR levels from 20c (20n_kmer_hmm.py = copy of step 20 + --level_dir).
#   Same pools, rho grid, scale grid and d_spread as the current run (SPREAD=1). Output: ${DNEW}/hmm/<SAMPLE>.npz
#   Run from new_model/ after 20c (env with numpy + scipy):
#     N=$(awk -F'\t' 'NR>1 && ($2=="control"||$2=="donor") && $3!="Zm"' <D11>/pools.tsv | wc -l)
#     sbatch --array=1-$(( (N + 29) / 30 )) 20n_run_kmer_hmm.sh
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/../final_script/config.sh"
need python3
python3 -c "import numpy, scipy" 2> /dev/null || die "python3 here has no numpy/scipy"
DNEW=${DNEW:-${KMER_DIR}/new_model}
[[ -n $(ls "${DNEW}/levels/"*.npy 2> /dev/null) ]] || die "no levels in ${DNEW}/levels (run 20c first)"
IDX=${SLURM_ARRAY_TASK_ID:-${1:?usage: bash $0 <task index>}}
OUT=${DNEW}/hmm; mkdir -p "$OUT"
mapfile -t POOLS < <(awk -F'\t' -v skip=" ${KMER_SKIP_GROUPS} " 'NR > 1 && ($2 == "control" || $2 == "donor") && index(skip, " " $3 " ") == 0 {print $1}' "${D11}/pools.tsv")
FIRST=$(( (IDX - 1) * KMER_HMM_POOLS_PER_TASK ))
(( FIRST < ${#POOLS[@]} )) || { echo "task ${IDX}: no pools"; exit 0; }
TODO=()
for (( i = FIRST; i < FIRST + KMER_HMM_POOLS_PER_TASK && i < ${#POOLS[@]}; i++ )); do
  [[ -s ${OUT}/${POOLS[$i]}.npz ]] || TODO+=("${POOLS[$i]}")
done
(( ${#TODO[@]} )) || { echo "task ${IDX}: all done already"; exit 0; }
ASMS=$(for a in "${DONOR_ASSEMBLIES[@]}"; do IFS=: read -r n t _ <<< "$a"; echo "${n}:${t}"; done | paste -sd,)
python3 "${SLURM_SUBMIT_DIR:-$PWD}/20n_kmer_hmm.py" --d_spread --level_dir "${DNEW}/levels" \
  --samples "$(IFS=,; echo "${TODO[*]}")" --pools_tsv "${D11}/pools.tsv" --win_dir "$D16" --ref_dir "$D19" \
  --assemblies "$ASMS" --map_dir "$MAP_DIR" --out_dir "$OUT" --plants "$KMER_PLANTS"
