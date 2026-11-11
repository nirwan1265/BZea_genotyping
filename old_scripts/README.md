# old_scripts

This folder keeps the earlier genotyping and introgression pipeline for reference. Nothing here is used by the current pipeline in [`../final_script/`](../final_script/README.md).

- `HMM_introgression/` contains the earlier pipeline:
  - ANGSD discovery (`0.0`–`0.5`)
  - BCF filter (`0.6`)
  - per-sample GL extraction (`1`, `1.2`)
  - the 3-state GL-based HMM (`hmm_viterbi_bc2s3.py`, `hmm_viterbi_rigidity.py`)
  - Optuna tuning
  - R post-processing

  Its own documentation is in [`HMM_introgression/README.md`](HMM_introgression/README.md).

## Why it was replaced

These points are covered in more detail in `final_script/README.md`.

- **Pools treated as single plants.** Each sample is a pool of plants from one family, but the HMM modeled it as one diploid plant (RR/RH/HH).
- **Alleles not oriented to B73.** ANGSD's BCF alleles are major/minor, not B73 REF/ALT. Wherever the major allele was not the B73 base, "RR" meant donor homozygous.
- **Donor-specific SNPs removed.** The `0.6` filter (population-wide AF 0.05–0.95) removed most SNPs found in only one donor's families.
- **HMM bugs:**
  - Chromosomes ran through one forward–backward pass, with a sticky transition at chromosome boundaries.
  - `thetas` started as `[0.0]`, which shifts every transition by one marker.
- **Inconsistent parameters.** The README, `best_params.txt` and `2_run_hmm_bySample.sh` each list different values.

Steps `0.0`–`0.5` map to `final_script/01`–`06`.
