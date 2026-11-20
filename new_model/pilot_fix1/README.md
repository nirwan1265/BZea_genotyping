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
sbatch --dependency=afterok:$J --array=1-70%100 p1_pool_windows.sh              # ~30-60 min per pool; ~25 GB temp disk each

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

## Result of the pilot (RULE=all), 19 Nov 2026

| Hotspot | Contrast old → new | SEG-island windows closer to the DONOR level, old → new |
|---|---|---|
| *parviglumis* chr7 54–65 Mb | 0.46 → **3.16** | 49% → **93%** |
| *mexicana* chr7 54–65 Mb | 1.31 → 1.94 | 15% → **87%** |
| *mexicana* chr3 75–78 Mb | 2.45 → 3.36 | 6% → **71%** |
| *parviglumis* chr3 75–78 Mb | 1.85 → no estimate | – |
| rest of chr3 + chr7 | ~3.0 → 3.1–3.4 | – |

Where the informative k-mers exist, the SEG islands at the hotspots read as DONOR, as the targeted sequencing says (figure `../../data/new_model/island_examples.png`: the long chr7 islands sit flat at the DONOR level in the new counts). The old counts are useless at chr7 in *parviglumis*: the teosinte keeps 90% of the B73 k-mers there, so DONOR and B73 are only 0.1 apart.

**The gap:** under RULE=all, *parviglumis* chr3 75–79 Mb keeps only 34–505 informative k-mers per window (*mexicana*: 3,100–13,900), and parts of chr7 55–78 Mb too, so those windows drop out (≥ 200 needed). One or both *parviglumis* assemblies carry nearly all the B73 k-mers there, while our *parviglumis* donors keep only 67%: the assemblies wrongly mark as shared the very k-mers our donors lack.

## Result of RULE=any for *parviglumis*, 19 Nov 2026

481,508,596 informative k-mers (RULE=all: 292,099,170). *parviglumis* chr3 75–79 Mb now has 9,000–18,000 per window (RULE=all: 34–505).

| Hotspot (*parviglumis*) | Contrast old → RULE=all → RULE=any | Island windows closer to DONOR, old → all → any |
|---|---|---|
| chr7 54–65 Mb | 0.46 → 3.16 → 2.10 | 49% → 93% → 70% |
| chr3 75–78 Mb | 1.85 → no estimate → 2.97 | 9% → – → 2% |
| rest of chr3 + chr7 | 3.12 → – → 3.50 | – |

Per pool (`result_any/pool_scores.tsv`; score 0 = DONOR level, 1 = B73):
- **chr7:** most island pools move to 0–0.1 in the new counts, including PN13_SID1159, which the targeted sequencing calls BB. Pools with a targeted BB call score 0 to −0.1; the one AB pool scores 0.76. The fix works here. RULE=all is better than RULE=any, because some of the extra k-mers are carried by our donors.
- **chr3:** the three island pools score 0.55, 0.51 and 1.14 in **both** old and new counts. The pools called DONOR score about 0, and the SEG pools (targeted AB) about 0.75–0.87. A second k-mer set with higher contrast gives the same scores, so the islands are not a contrast problem: the reads put them between DONOR and SEG, or at B73. The one with targeted sequencing (PN16_SID1505) is BB. Targeted calls at *mexicana* chr3 islands are mixed (AA, BB and AB). Fix 1 does not resolve chr3.
- Two pools look odd in both counts: PN15_SID1390 (above B73 at chr7, called DONOR at chr3 but scoring 0.64) and PN12_SID1086 (called DONOR at chr7, scoring 0.65).

## RULE=any: how it was run

RULE=any keeps the k-mers missing from **at least one** assembly (B73_uniq minus the k-mers all assemblies share). More k-mers survive, including those our donors lack; some will be present in our donors, which is handled because the DONOR level is learned from our own pools (20c) in the full build.

```bash
cd /share/maize/ntanduk/BZea_genotyping/new_model/pilot_fix1
conda activate /rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk/kmer_env
J=$(sbatch --parsable --export=ALL,RULE=any --array=1 p1_build_sets.sh)                       # Zv only
sbatch --dependency=afterok:$J --export=ALL,RULE=any,TAXA=Zv --array=1-70%100 p1_pool_windows.sh   # Zv pools + controls
# when done (machine_learning env):
python p1_contrast.py --pilot_dir $KMER_DIR/new_model/pilot_fix1 --pools pools.tsv --d16 $D16 --d19 $D19 \
  --hmm_dir $D20 --tracts $D21/tracts_all.tsv.gz --out $KMER_DIR/new_model/pilot_fix1/result_any --rule any --taxa Zv
```
