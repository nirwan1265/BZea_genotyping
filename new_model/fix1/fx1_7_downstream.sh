#!/bin/bash
# =============================================================================
# fix1 / 7 — everything downstream of the calls, on the new model (MODEL=fix1, set in final_script/config.sh):
#   22  genome shares, per chromosome           22b SEG donor copies (all SEG windows, and the 5-Mb core check)
#   23  introgression viewer                     24  result figures (after 22)
#   25-27, 29.2  k-mer GWAS + permutations + joint QTL, taxon and mexicana races (submit_assoc_all.sh)
#   the same per year (submit_year_check.sh) and for predicted leaf N 2025 (submit_trait_run.sh)
# All outputs go next to the new calls: ${KMER_DIR}/new_model/fix1/with_levels/introgression/ (the current run's
# folders are not touched). Run on the login node with the GWAS env active (R jobs inherit it); the python jobs of
# 22-24 get the machine_learning env from fx1_config.sh:
#   conda activate /usr/local/usrapps/maize/ntanduk/bzea_gwas
#   cd /share/maize/ntanduk/BZea_genotyping/new_model/fix1 && bash fx1_7_downstream.sh
# =============================================================================
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
PY_ENV=$(grep '^PY_ENV=' "${HERE}/fx1_config.sh" | cut -d= -f2 | cut -d' ' -f1)
cd "${HERE}/../../final_script"
export MODEL=fix1
source config.sh
mkdir -p logs
[[ -s ${D21}/tracts_all.tsv.gz ]] || die "no calls in ${D21} (fx1_submit_levels.sh first)"
Rscript -e 'library(GridLMM); library(qtl2)' > /dev/null 2>&1 || die "GridLMM/qtl2 not found: activate the GWAS conda env first"
PYX="--export=ALL,PATH=${PY_ENV}/bin:${PATH}"
echo "model ${MODEL}: calls ${D21}"
J22=$(sbatch --parsable $PYX 22_kmer_stats.sh)
J22B=$(sbatch --parsable $PYX 22b_seg_donor_copies.sh)
J22C=$(sbatch --parsable "${PYX},CORE_MB=5,TRIM_MB=1" 22b_seg_donor_copies.sh)
J23=$(sbatch --parsable $PYX 23_introgression_viewer.sh)
J24=$(sbatch --parsable $PYX --dependency=afterok:${J22} 24_kmer_figures.sh)
echo "22 ${J22} | 22b ${J22B} ${J22C} | 23 ${J23} | 24 ${J24}"
bash submit_assoc_all.sh
bash submit_year_check.sh
bash submit_trait_run.sh
echo "follow: squeue -u \$USER; outputs in ${D21}"
