#!/bin/bash
# =============================================================================
# submit_counts.sh — submit steps 07 -> 08 -> 09 -> 10 as one Slurm chain.
#   Each step starts only if every task of the previous step succeeded (afterok).
#
#   Run from inside final_script/ on a login node:
#       bash submit_counts.sh test   # 07, then only part 1 of 08 (check run time and memory first)
#       bash submit_counts.sh rest   # after the test: parts 2-200 of 08, then 09 and 10
#       bash submit_counts.sh        # everything in one go, no test
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs
mode=${1:-all}

case $mode in
  test)
    j07=$(sbatch --parsable 07_make_site_tables.sh)
    j08=$(sbatch --parsable --dependency=afterok:"${j07}" --array=1 08_angsd_counts_by_part.sh)
    echo "07 sites          job ${j07}"
    echo "08 counts part 1  job ${j08}"
    echo "When it finishes, check logs/08_counts.${j08}_1.out (ends with DONE) and"
    echo "  sacct -j ${j08} --format=JobID,Elapsed,MaxRSS,State"
    echo "then run: bash submit_counts.sh rest"
    exit 0 ;;
  rest)
    j08=$(sbatch --parsable --array=2-200%40 08_angsd_counts_by_part.sh) ;;
  all)
    j07=$(sbatch --parsable 07_make_site_tables.sh)
    echo "07 sites          job ${j07}"
    j08=$(sbatch --parsable --dependency=afterok:"${j07}" 08_angsd_counts_by_part.sh) ;;
  *)
    echo "usage: bash submit_counts.sh [test|rest]" >&2; exit 1 ;;
esac

echo "08 counts         job ${j08}"
j09=$(sbatch --parsable --dependency=afterok:"${j08}" 09_concat_counts_by_chr.sh)
echo "09 concat         job ${j09}"
j10=$(sbatch --parsable --dependency=afterok:"${j09}" 10_counts_qc.sh)
echo "10 QC             job ${j10}"
echo "Watch: squeue -u \$USER    If any 08 part fails: bash 08_check_unfinished.sh"
