#!/bin/bash
#SBATCH --job-name=tau_hmm
#SBATCH --output=logs/tau_1_hmm.%A_%a.out
#SBATCH --error=logs/tau_1_hmm.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --mem=12G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# tau_run / 1 — rerun the current model's HMM (fix1, with_levels) on the pools of tau_samples.tsv, one pool per
#   array task, with a floor on tau: every taxon's assembly tau is raised to at least TAU_MULT x the tau of
#   TAU_FLOOR (default Zv, the widest). Zh and Zl (one assembly, tau = 0) get Zv's tau; Zd and Zx get max(own, Zv);
#   Zv is unchanged, so its pools must come out the same as the current run (a check that the rerun is exact).
#   Everything else (pools, levels, switched windows, rho grid, scale grid) is the current run's.
#   Output: ${TAU_OUT}/hmm/<SAMPLE>.npz     (TAU_OUT = ${KMER_DIR}/new_model/tau_run/<TAU_NAME>)
# Submit (login node, python env with numpy/scipy active):
#   cd /share/maize/ntanduk/BZea_genotyping/new_model/tau_run && mkdir -p logs
#   sbatch --array=1-25%100 tau_1_hmm.sh
# Half the floor instead (a second run next to the first):
#   TAU_NAME=zv_half TAU_MULT=0.5 sbatch --export=ALL,TAU_NAME,TAU_MULT --array=1-25%100 tau_1_hmm.sh
# =============================================================================
set -euo pipefail
HERE=$(cd "${SLURM_SUBMIT_DIR:-$PWD}" && pwd)
unset MODEL
source "${HERE}/../../final_script/config.sh"
export FX_RUN=with_levels FX_LEVEL_DIR=${KMER_DIR}/new_model/levels
source "${HERE}/../fix1/fx1_config.sh"
FINAL_DIR=$(cd "${HERE}/../../final_script" && pwd); MAP_DIR=${FINAL_DIR}/genetic_maps
TAU_FLOOR=${TAU_FLOOR:-Zv}; TAU_MULT=${TAU_MULT:-1}; TAU_NAME=${TAU_NAME:-zv_floor}
TAU_OUT=${KMER_DIR}/new_model/tau_run/${TAU_NAME}
python -c "import numpy, scipy" 2> /dev/null || die "python here has no numpy/scipy: activate your python env before submitting"
IDX=${SLURM_ARRAY_TASK_ID:-${1:?usage: bash $0 <task index>}}
S=$(awk -F'\t' -v i=$((IDX + 1)) 'NR == i {print $1}' "${HERE}/tau_samples.tsv")
[[ -n $S ]] || { echo "task ${IDX}: no pool"; exit 0; }
mkdir -p "${TAU_OUT}/hmm"
[[ -s ${TAU_OUT}/hmm/${S}.npz ]] && { echo "${S}: done already"; exit 0; }
ASMS=$(for a in "${DONOR_ASSEMBLIES[@]}"; do IFS=: read -r n t _ <<< "$a"; echo "${n}:${t}"; done | paste -sd,)
echo "pool ${S} | floor ${TAU_MULT} x ${TAU_FLOOR} tau | out ${TAU_OUT}"
python "${HERE}/tau_kmer_hmm.py" --d_spread --tau_floor_taxon "$TAU_FLOOR" --tau_floor_mult "$TAU_MULT" --tau_floor_for all \
  --level_dir "$FX_LEVEL_DIR" --fix_dir "${FX}/choice" --new_win_dir "${FX}/windows" --new_ref_dir "${FX}/sets" \
  --min_new "$FX_MIN_NEW" --samples "$S" --pools_tsv "${D11}/pools.tsv" --win_dir "$D16" --ref_dir "$D19" \
  --assemblies "$ASMS" --map_dir "$MAP_DIR" --out_dir "${TAU_OUT}/hmm" --plants "$KMER_PLANTS"
