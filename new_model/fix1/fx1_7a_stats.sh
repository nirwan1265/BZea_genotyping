#!/bin/bash
# =============================================================================
# fix1 / 7a — summaries of the new calls (MODEL=fix1): 22 genome shares, 22b SEG donor copies (all SEG windows and
#   the 5-Mb core check), 23 introgression viewer, 24 result figures. Outputs next to the calls:
#   ${KMER_DIR}/new_model/fix1/with_levels/introgression/. Jobs use the env active when you submit (numpy, scipy, pandas,
#   matplotlib). Run on the login node:  cd new_model/fix1 && bash fx1_7a_stats.sh
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")/../../final_script"
export MODEL=fix1
source config.sh
mkdir -p logs
[[ -s ${D21}/tracts_all.tsv.gz ]] || die "no calls in ${D21} (fx1_submit_levels.sh first)"
J22=$(sbatch --parsable 22_kmer_stats.sh)
J22B=$(sbatch --parsable 22b_seg_donor_copies.sh)
J22C=$(sbatch --parsable --export=ALL,CORE_MB=5,TRIM_MB=1 22b_seg_donor_copies.sh)
J23=$(sbatch --parsable 23_introgression_viewer.sh)
J24=$(sbatch --parsable --dependency=afterok:${J22} 24_kmer_figures.sh)
echo "calls ${D21}: 22 ${J22} | 22b ${J22B} ${J22C} | 23 ${J23} | 24 ${J24}"
