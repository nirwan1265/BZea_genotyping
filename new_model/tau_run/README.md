# tau_run — does a floor on τ change the calls?

**Question.** τ_w is how unsure the model is about a taxon's DONOR level in a window (main README, 2.6.1). It is measured from how much the taxon's two assemblies disagree. *huehuetenangensis* (Zh) and *luxurians* (Zl) have one assembly each, so their τ is 0 outside the switched windows, as if their DONOR level were known exactly. If that matters, giving them a realistic τ should move some windows between SEG and DONOR.

**Test.** Rerun the current model's HMM (fix1, `with_levels`) on 25 pools with every taxon's assembly τ raised to at least *parviglumis*' τ (the widest; the B73 controls already use it). Nothing else changes.
- Zh, Zl: τ 0 → Zv's τ. Zd, Zx: max(own, Zv's). **Zv: unchanged**, so its 5 pools must come out the same as the current run, a check that the rerun is exact.
- Switched windows (fix1) keep the τ learned from our pools.
- Pools (`tau_samples.tsv`): per taxon the 5 pools with the most clearly called SEG/DONOR tracts in the targeted sequencing (three bad samples left out); Zh has only 3 such pools, so 2 more Zh pools with the most DONOR are added (no targeted data).

**What to look at.** How many windows change state (and which way), and whether the targeted sequencing agrees more or less with the new calls in the same pools.

## Run

```bash
# 1. cluster (login node), python env with numpy/scipy active
cd /share/maize/ntanduk/BZea_genotyping && git pull && cd new_model/tau_run && mkdir -p logs
sbatch --array=1-25%100 tau_1_hmm.sh
# 2. cluster, when the 25 tasks are done
source ../../final_script/config.sh
python tau_2_collect.py --samples tau_samples.tsv --old_hmm /share/maize/ntanduk/kmer/new_model/fix1/with_levels/hmm \
  --new_hmm /share/maize/ntanduk/kmer/new_model/tau_run/zv_floor/hmm \
  --rho_tsv /share/maize/ntanduk/kmer/new_model/fix1/with_levels/introgression/rho_loglik.tsv \
  --ref_dir $D19 --out /share/maize/ntanduk/kmer/new_model/tau_run/zv_floor/compare
# 3. laptop, from the repo root
mkdir -p data/tau_run/zv_floor
rsync -av ntanduk@login.hpc.ncsu.edu:/share/maize/ntanduk/kmer/new_model/tau_run/zv_floor/compare/ data/tau_run/zv_floor/
python new_model/tau_run/tau_3_targseq.py --run_dir data/tau_run/zv_floor
```

Files: `tau_kmer_hmm.py` (copy of `new_model/fix1/fx1_kmer_hmm.py` plus `--tau_floor_taxon`, `--tau_floor_for`, `--tau_floor_mult`), `tau_1_hmm.sh` (one pool per array task), `tau_2_collect.py` (old vs new per pool: `pools.tsv`, `changed_regions.tsv`, both runs' tracts), `tau_3_targseq.py` (step 34d on both, per taxon: `targseq_by_taxon.tsv`).
