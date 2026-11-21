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

## Scale runs: is the per-pool scale a a problem?

**Question.** a is one number per pool that multiplies every DONOR level of the pool (DONOR expected at a × d_w). It is not measured: the HMM is run with each a of a grid (0.5, 0.6, …, 1.4, 1.6, 1.8) and the best fit is kept. In the current run, 20% of *diploperennis* pools (45 of 227) end at an edge of the grid (0.5 or 1.8), against 1–7% in the other taxa. An edge means the best a may lie outside the grid, or that a is not pinned down by the data and trades off with the SEG/DONOR call (in the τ run, `PN15_SID1348` went from a = 1.8 with a chr3 DONOR tract to a = 0.5 with the same stretch SEG).

**Test.** Two reruns of the current model (current τ, no floor) on 35 pools (`scale_samples.tsv`: the 25 pools above plus 10 Zd pools at an edge, 5 at 0.5 and 5 at 1.8, those with the most SEG + DONOR):
- `scale_wide`: grid 0.3–2.5. If edge pools settle inside and their calls hardly change, the edges are harmless.
- `scale_1`: a fixed at 1. The current model learns each window's DONOR level from our pools, which covers most of what a was for; if a = 1 agrees with the targeted sequencing as well or better, a may not be needed.

```bash
# 1. cluster (login node), python env active
cd /share/maize/ntanduk/BZea_genotyping && git pull && cd new_model/tau_run && mkdir -p logs
TAU_NAME=scale_wide TAU_FLOOR=none TAU_SAMPLES=scale_samples.tsv SCALE_GRID=0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.6,1.8,2.0,2.2,2.5 \
  sbatch --export=ALL,TAU_NAME,TAU_FLOOR,TAU_SAMPLES,SCALE_GRID --array=1-35%100 tau_1_hmm.sh
TAU_NAME=scale_1 TAU_FLOOR=none TAU_SAMPLES=scale_samples.tsv SCALE_GRID=1 \
  sbatch --export=ALL,TAU_NAME,TAU_FLOOR,TAU_SAMPLES,SCALE_GRID --array=1-35%100 tau_1_hmm.sh
# 2. cluster, when both are done
source ../../final_script/config.sh
for R in scale_wide scale_1; do
  python tau_2_collect.py --samples scale_samples.tsv --old_hmm /share/maize/ntanduk/kmer/new_model/fix1/with_levels/hmm \
    --new_hmm /share/maize/ntanduk/kmer/new_model/tau_run/$R/hmm \
    --rho_tsv /share/maize/ntanduk/kmer/new_model/fix1/with_levels/introgression/rho_loglik.tsv \
    --ref_dir $D19 --out /share/maize/ntanduk/kmer/new_model/tau_run/$R/compare
done
# 3. laptop, from the repo root
for R in scale_wide scale_1; do
  mkdir -p data/tau_run/$R
  rsync -av ntanduk@login.hpc.ncsu.edu:/share/maize/ntanduk/kmer/new_model/tau_run/$R/compare/ data/tau_run/$R/
  python new_model/tau_run/tau_3_targseq.py --run_dir data/tau_run/$R
done
```
