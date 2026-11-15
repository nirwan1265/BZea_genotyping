#!/bin/bash
#SBATCH --job-name=kmer_hmm
#SBATCH --output=logs/20_kmer_hmm.%A_%a.out
#SBATCH --error=logs/20_kmer_hmm.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH --mem=8G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# 20 — HMM on k-mer windows, KMER_HMM_POOLS_PER_TASK pools per array task (20_kmer_hmm.py).
#   Each pool is run over the whole rho grid (log-likelihood saved); step 21 picks rho.
#   Inputs: ${D16} (step 17), ${D19} (step 19), ${D11}/pools.tsv, genetic_maps/
#   Output: ${D20}/<SAMPLE>.npz
#   Needs python3 with numpy + scipy (check: python3 -c "import numpy, scipy").
#   Submit after 17 and 19:
#     N=$(awk -F'\t' 'NR>1 && ($2=="control"||$2=="donor") && $3!="Zm"' <D11>/pools.tsv | wc -l)
#     sbatch --array=1-$(( (N + 29) / 30 )) 20_run_kmer_hmm.sh
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"
need python3
python3 -c "import numpy, scipy" 2> /dev/null || die "python3 here has no numpy/scipy"

IDX=${SLURM_ARRAY_TASK_ID:-${1:?usage: bash $0 <task index>}}
mkdir -p "$D20"
# pools: controls + donors, minus taxa in KMER_SKIP_GROUPS
mapfile -t POOLS < <(awk -F'\t' -v skip=" ${KMER_SKIP_GROUPS} " 'NR > 1 && ($2 == "control" || $2 == "donor") && index(skip, " " $3 " ") == 0 {print $1}' "${D11}/pools.tsv")
FIRST=$(( (IDX - 1) * KMER_HMM_POOLS_PER_TASK ))
(( FIRST < ${#POOLS[@]} )) || { echo "task ${IDX}: no pools"; exit 0; }
TODO=()
for (( i = FIRST; i < FIRST + KMER_HMM_POOLS_PER_TASK && i < ${#POOLS[@]}; i++ )); do
  [[ -s ${D20}/${POOLS[$i]}.npz ]] || TODO+=("${POOLS[$i]}")
done
(( ${#TODO[@]} )) || { echo "task ${IDX}: all done already"; exit 0; }
ASMS=$(for a in "${DONOR_ASSEMBLIES[@]}"; do IFS=: read -r n t _ <<< "$a"; echo "${n}:${t}"; done | paste -sd,)
python3 "${FINAL_DIR}/20_kmer_hmm.py" --samples "$(IFS=,; echo "${TODO[*]}")" \
  --pools_tsv "${D11}/pools.tsv" --win_dir "$D16" --ref_dir "$D19" --assemblies "$ASMS" \
  --map_dir "$MAP_DIR" --out_dir "$D20" --plants "$KMER_PLANTS"
