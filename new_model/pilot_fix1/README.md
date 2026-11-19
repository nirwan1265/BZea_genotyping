# Pilot of fix 1: count only the informative k-mers

**The problem.** In windows where the teosinte keeps most of the B73 k-mers (chr7 54–65 Mb and chr3 75–78 Mb keep 61–78% in *parviglumis* and *mexicana*), DONOR and SEG score almost the same, and the HMM calls SEG islands inside DONOR blocks (main README 2.9.1, point 7 and "contrast"). Each window's score counts all B73 single-copy k-mers, but the k-mers the teosinte shares with B73 are present in every pool, whatever its state: they add counting noise and no signal.

**The idea.** Count only the **informative** k-mers: B73 single-copy k-mers that are in **none** of the taxon's assemblies. The signal (k-mers going missing in DONOR and SEG) stays; the noise of the shared k-mers goes. In theory the contrast rises by up to 1/√(1 − d_w): about 2× where the teosinte keeps 75% (chr7), hardly at all where it keeps 10%. Part of the noise is per window and does not shrink, so the real gain is smaller; this pilot measures it before anything is built for the whole population.

**The pilot.** 60 pools of *parviglumis* and *mexicana* (`pools.tsv`), chosen at the two hotspots: 30 with a SEG island between two DONOR tracts there (the error; 11 of them have targeted sequencing), 14 with DONOR there (to measure the DONOR level), 8 with SEG there, and 8 more *parviglumis* pools with targeted sequencing; plus 10 B73 controls. Only chr3 and chr7 are counted.

| Step | Script | Env | What it does |
|---|---|---|---|
| 1 | `p1_build_sets.sh` (array 1–2: Zv, Zx) | k-mer env (kmc, kmc_tools, kcftools) | informative k-mers per taxon = B73_uniq minus every k-mer of the taxon's two assemblies (`kmc_tools kmers_subtract`); chr3 + chr7 reference; per-window count of informative k-mers in B73 |
| 2 | `p1_pool_windows.sh` (array 1–70) | k-mer env | per pool: count reads with KMC (as step 17), keep the informative k-mers of its taxon (controls: both taxa), windows on chr3 + chr7 |
| 3 | `p1_contrast.py` | numpy, scipy, pandas, matplotlib (machine_learning) | scores old (all k-mers) and new (informative only) with step 20's normalisation; contrast per window; whether the SEG-island windows sit closer to the DONOR or the SEG level; figure |

**Test on simulated pools** (DONOR keeps 76% of the B73 k-mers at a hotspot, binomial counting noise and 2% per-window noise): the contrast at the hotspot rises from 15.5 to 26.8 (1.7×; theory 2.0×), and the rest changes little (38.6 to 46.0, where the teosinte keeps 30%).

## Run (from `new_model/pilot_fix1/` on the cluster)

```bash
cd /share/maize/ntanduk/BZea_genotyping/new_model/pilot_fix1 && mkdir -p logs
conda activate /rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk/kmer_env      # kmc, kmc_tools, kcftools
J=$(sbatch --parsable --array=1-2 p1_build_sets.sh)                             # ~1-3 h
sbatch --dependency=afterok:$J --array=1-70%10 p1_pool_windows.sh               # ~30-60 min per pool, 10 at a time

# when all 70 are done:
conda activate /rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk/machine_learning
eval "$(bash -c 'source ../../final_script/config.sh; echo D16=$D16 D19=$D19 D20=$D20 D21=$D21 KMER_DIR=$KMER_DIR')"
python p1_contrast.py --pilot_dir $KMER_DIR/new_model/pilot_fix1 --pools pools.tsv --d16 $D16 --d19 $D19 \
  --hmm_dir $D20 --tracts $D21/tracts_all.tsv.gz --out $KMER_DIR/new_model/pilot_fix1/result
```

Outputs: `result/contrast_old_new.tsv` (printed too) and `result/island_examples.png`.

**How to read the result.**
- `contrast_new` clearly above `contrast_old` at chr7 54–65 and chr3 75–78 Mb (e.g. 1.5× or more): the fix works; build it for all pools.
- `island_windows_closer_to_DONOR` rising from old to new: the SEG islands move towards the DONOR level, i.e. the informative k-mers see them as DONOR, as the targeted sequencing does.
- Little change: the noise in those windows is mostly per window (not counting noise), and only counting both sides (teosinte k-mers too, fix 2) or deeper sequencing would help.

Checks before step 2: `grep "informative k-mers" logs/p1_build_sets.*.out` should give a count per taxon (a share of the ~680 M B73 single-copy k-mers), and `zcat $KMER_DIR/new_model/pilot_fix1/sets/B73_inf_Zv.tsv.gz | head` should show chr3/chr7 windows with counts.
