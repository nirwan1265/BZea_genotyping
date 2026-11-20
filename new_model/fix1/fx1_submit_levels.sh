#!/bin/bash
# =============================================================================
# fix1 — second run of steps 4-6 only (no recount): fix 1 in the switched windows AND new_model/20c's taxon
#   DONOR levels (learned from all pools) in the windows that keep the old counts. 20c alone fixed the chr6 21-23 and
#   chr1 297 Mb hotspots (48 -> 11, 28 -> 14 pools), which fix 1 does not touch. Outputs: ${FX}/with_levels/
#     cd /share/maize/ntanduk/BZea_genotyping/new_model/fix1 && bash fx1_submit_levels.sh
# =============================================================================
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); cd "$HERE"
source ../../final_script/config.sh
export FX_RUN=with_levels FX_LEVEL_DIR=${KMER_DIR}/new_model/levels
source ./fx1_config.sh
[[ -n $(ls "${FX_LEVEL_DIR}"/*.npy 2> /dev/null) ]] || die "no 20c levels in ${FX_LEVEL_DIR}"
[[ -n $(ls "${FX}/choice/"*.npz 2> /dev/null) ]] || die "no fix-1 choices in ${FX}/choice"
mkdir -p logs "$FXR"
NH=$(awk -F'\t' -v skip=" ${KMER_SKIP_GROUPS} " 'NR > 1 && ($2 == "control" || $2 == "donor") && index(skip, " " $3 " ") == 0' "${D11}/pools.tsv" | wc -l)
NH=$(( (NH + KMER_HMM_POOLS_PER_TASK - 1) / KMER_HMM_POOLS_PER_TASK ))
J4=$(sbatch --parsable --export=ALL,FX_RUN,FX_LEVEL_DIR --array=1-${NH}%100 fx1_4_run_hmm.sh)
J5=$(sbatch --parsable --export=ALL,FX_RUN,FX_LEVEL_DIR --dependency=afterok:${J4} fx1_5_collect.sh)
J6=$(sbatch --parsable --export=ALL,FX_RUN,FX_LEVEL_DIR --dependency=afterok:${J5} fx1_6_checks.sh)
echo "submitted: hmm ${J4} -> calls ${J5} -> checks ${J6}; outputs in ${FXR}"
