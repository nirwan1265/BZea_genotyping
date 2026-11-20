#!/bin/bash
#SBATCH --job-name=fx1_hmm
#SBATCH --output=logs/fx1_4_hmm.%A_%a.out
#SBATCH --error=logs/fx1_4_hmm.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00:00
#SBATCH --mem=12G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# fix1 / 4 — the HMM with the switched windows (fx1_kmer_hmm.py = 20n + --fix_dir). Same pools, rho grid, scale
#   grid, d_spread and assembly levels as the current run (SPREAD=1) in all other windows.
#   Array size: ceil(pools / KMER_HMM_POOLS_PER_TASK). Output: ${FX}/hmm/<SAMPLE>.npz
# =============================================================================
set -euo pipefail
HERE=$(cd "${SLURM_SUBMIT_DIR:-$PWD}" && pwd)
source "${HERE}/../../final_script/config.sh"; source "${HERE}/fx1_config.sh"; use_py_env
FINAL_DIR=$(cd "${HERE}/../../final_script" && pwd); MAP_DIR=${FINAL_DIR}/genetic_maps
python3 -c "import numpy, scipy" 2> /dev/null || die "python3 in ${PY_ENV} has no numpy/scipy"
[[ -n $(ls "${FX}/choice/"*.npz 2> /dev/null) ]] || die "no choices in ${FX}/choice (run fx1_3 first)"
IDX=${SLURM_ARRAY_TASK_ID:-${1:?usage: bash $0 <task index>}}
OUT=${FX}/hmm; mkdir -p "$OUT"
mapfile -t POOLS < <(awk -F'\t' -v skip=" ${KMER_SKIP_GROUPS} " 'NR > 1 && ($2 == "control" || $2 == "donor") && index(skip, " " $3 " ") == 0 {print $1}' "${D11}/pools.tsv")
FIRST=$(( (IDX - 1) * KMER_HMM_POOLS_PER_TASK ))
(( FIRST < ${#POOLS[@]} )) || { echo "task ${IDX}: no pools"; exit 0; }
TODO=()
for (( i = FIRST; i < FIRST + KMER_HMM_POOLS_PER_TASK && i < ${#POOLS[@]}; i++ )); do
  [[ -s ${OUT}/${POOLS[$i]}.npz ]] || TODO+=("${POOLS[$i]}")
done
(( ${#TODO[@]} )) || { echo "task ${IDX}: all done already"; exit 0; }
ASMS=$(for a in "${DONOR_ASSEMBLIES[@]}"; do IFS=: read -r n t _ <<< "$a"; echo "${n}:${t}"; done | paste -sd,)
python3 "${HERE}/fx1_kmer_hmm.py" --d_spread --fix_dir "${FX}/choice" --new_win_dir "${FX}/windows" --new_ref_dir "${FX}/sets" \
  --min_new "$FX_MIN_NEW" --samples "$(IFS=,; echo "${TODO[*]}")" --pools_tsv "${D11}/pools.tsv" --win_dir "$D16" --ref_dir "$D19" \
  --assemblies "$ASMS" --map_dir "$MAP_DIR" --out_dir "$OUT" --plants "$KMER_PLANTS"
