#!/bin/bash
#SBATCH --job-name=collect_hmm_family
#SBATCH --output=logs/15b_collect_hmm_family.%j.out
#SBATCH --error=logs/15b_collect_hmm_family.%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --mem=16G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# 15b — Collect the family-informed HMM results (step 14b) (wrapper for 15_collect_hmm.R).
#   Output in ${D15B}: pool_summary.tsv, tracts_all.tsv.gz, family_check.tsv
#
#   family_check.tsv compares each pool with step 15 (panel-only run). If a pool's family
#   label is wrong, its DONOR tracts turn into SEG and FRAC_DONOR drops (simulation: from
#   0.26-0.30 to 0). Pools flagged CHECK: confirm the pedigree; to run them without family
#   info, list them in sample_lists/no_family.txt, delete their ${D14B}/<SAMPLE>.* files and
#   resubmit 14b.
#   Submit after 14b (and after 15 for the check):  sbatch 15b_collect_hmm_family.sh
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"
need Rscript
Rscript "${FINAL_DIR}/15_collect_hmm.R" "$D14B" "${D11}/pools.tsv" "$D15B"

# family check: FRAC_DONOR in the family run vs the panel-only run (step 15)
V1=${D15}/pool_summary.tsv
if [[ -s $V1 ]]; then
  awk -F'\t' -v OFS='\t' '
    FNR == 1 { for (i = 1; i <= NF; i++) h[$i] = i; next }
    FILENAME == ARGV[1] { d1[$h["SAMPLE"]] = $h["FRAC_DONOR"]; s1[$h["SAMPLE"]] = $h["DONOR_SHARE"]; next }
    !hdr { print "SAMPLE", "GROUP", "FAMILY", "FAMILY_MARKERS", "FRAC_DONOR_PANEL", "FRAC_DONOR_FAMILY", "SHARE_PANEL", "SHARE_FAMILY", "FLAG"; hdr = 1 }
    ($h["SAMPLE"] in d1) && $h["FRAC_DONOR"] != "NA" && d1[$h["SAMPLE"]] != "NA" {
      x = $h["SAMPLE"]
      flag = (d1[x] > 0.02 && $h["FRAC_DONOR"] < 0.5 * d1[x]) ? "CHECK" : "ok"
      print x, $h["GROUP"], $h["FAMILY"], $h["FAMILY_MARKERS"], d1[x], $h["FRAC_DONOR"], s1[x], $h["DONOR_SHARE"], flag
    }' "$V1" "${D15B}/pool_summary.tsv" > "${D15B}/family_check.tsv"
  echo "family check: $(awk -F'\t' '$NF == "CHECK"' "${D15B}/family_check.tsv" | wc -l) pools flagged CHECK -> ${D15B}/family_check.tsv"
else
  echo "no ${V1} yet (run 15 first) -> family check skipped"
fi
