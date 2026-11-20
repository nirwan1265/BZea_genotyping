#!/bin/bash
# =============================================================================
# fix1 — submit the whole build as one chain of Slurm jobs (run on the login node, from new_model/fix1/):
#     cd /share/maize/ntanduk/BZea_genotyping/new_model/fix1 && bash fx1_submit.sh
#   1 sets (5 taxa) -> 2 pool windows (all pools, %100) -> 3 choose windows (5 taxa) -> 4 HMM (%100) -> 5 calls -> 6 checks
#   Each job loads its own env (fx1_config.sh), so no conda activate is needed. Step 2 continues even if a few pools
#   fail (they keep the old counts in the HMM); every other step needs the one before it to succeed.
#   Safe to rerun: finished pieces are skipped.
# =============================================================================
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); cd "$HERE"
source ../../final_script/config.sh; source ./fx1_config.sh
mkdir -p logs "$FX"
NP=$(fx_pools | wc -l)
NT=$(( (NP + FX_POOLS_PER_TASK - 1) / FX_POOLS_PER_TASK ))
NH=$(awk -F'\t' -v skip=" ${KMER_SKIP_GROUPS} " 'NR > 1 && ($2 == "control" || $2 == "donor") && index(skip, " " $3 " ") == 0' "${D11}/pools.tsv" | wc -l)
NH=$(( (NH + KMER_HMM_POOLS_PER_TASK - 1) / KMER_HMM_POOLS_PER_TASK ))
echo "pools to count: ${NP} (${NT} tasks); HMM tasks: ${NH}; outputs in ${FX}"
J1=$(sbatch --parsable --array=1-5 fx1_1_sets.sh)
J2=$(sbatch --parsable --dependency=afterok:${J1} --array=1-${NT}%100 fx1_2_pool_windows.sh)
J3=$(sbatch --parsable --dependency=afterany:${J2} --array=1-5 fx1_3_choose.sh)
J4=$(sbatch --parsable --dependency=afterok:${J3} --array=1-${NH}%100 fx1_4_run_hmm.sh)
J5=$(sbatch --parsable --dependency=afterok:${J4} fx1_5_collect.sh)
J6=$(sbatch --parsable --dependency=afterok:${J5} fx1_6_checks.sh)
echo "submitted: sets ${J1} -> pools ${J2} -> choose ${J3} -> hmm ${J4} -> calls ${J5} -> checks ${J6}"
echo "${J1} ${J2} ${J3} ${J4} ${J5} ${J6}" > logs/fx1_jobs.txt
