#!/bin/bash
#SBATCH --job-name=fx1_checks
#SBATCH --output=logs/fx1_6_checks.%j.out
#SBATCH --error=logs/fx1_6_checks.%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --mem=24G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# fix1 / 6 — the checks of the new calls, as for the other new_model runs:
#   targeted sequencing (final_script/34c + 34d, unchanged) and old vs new (new_model/compare_old_new.py).
#   Output: ${FX}/targseq/, ${FX}/targseq_tracts/, ${FX}/compare_old_new.txt (also printed in the log)
# =============================================================================
set -euo pipefail
HERE=$(cd "${SLURM_SUBMIT_DIR:-$PWD}" && pwd)
source "${HERE}/../../final_script/config.sh"; source "${HERE}/fx1_config.sh"; use_py_env
FINAL_DIR=$(cd "${HERE}/../../final_script" && pwd); REPO=$(cd "${HERE}/../.." && pwd)
mkdir -p "${FXR}/targseq" "${FXR}/targseq_tracts"
python3 "${FINAL_DIR}/34c_targseq_compare.py" --vcf_dir "${D21}/targseq/vcf" --tracts "${FXR}/introgression/tracts_all.tsv.gz" \
  --pool_summary "${FXR}/introgression/pool_summary.tsv" --names "${REPO}/data/phenotypes/raw/name_change_final.txt" \
  --out "${FXR}/targseq" > "${FXR}/targseq/log.txt"
python3 "${FINAL_DIR}/34d_targseq_tracts.py" --sites "${FXR}/targseq/sites_used.tsv.gz" --tracts "${FXR}/introgression/tracts_all.tsv.gz" \
  --out "${FXR}/targseq_tracts" > "${FXR}/targseq_tracts/log.txt"
python3 "${HERE}/../compare_old_new.py" --old "${D21}/tracts_all.tsv.gz" --new "${FXR}/introgression/tracts_all.tsv.gz" \
  --pool_summary "${D21}/pool_summary.tsv" --targseq_old "${REPO}/data/introgression_spread/targseq/tracts_targseq.tsv" \
  --targseq_new "${FXR}/targseq_tracts/tracts_targseq.tsv" | tee "${FXR}/compare_old_new.txt"
echo; echo "== which windows switched (fx1_3)"; cat "${FX}"/choice/summary_*.tsv
