# new_model: DONOR levels learned from all pools

A test of a new model, kept apart from the main pipeline. **Nothing in `final_script/` is changed**: the scripts here read the current run's outputs and write to their own folder on the cluster (`${KMER_DIR}/new_model`, called `DNEW` below).

## Why

The targeted sequencing (main README, section 2.9.1, point 7) found one systematic error of the current HMM: in *parviglumis* and *mexicana*, SEG tracts that lie between two DONOR tracts are really DONOR in about half of the cases (47% BB, against 1.4% for SEG tracts with no DONOR neighbour). The errors sit at the same windows in many pools (chr6 21–23 Mb in 48 pools, chr1 297 Mb, chr7 54–65 Mb, chr3 75–78 Mb).

The cause is the DONOR level of a window (d_w), which comes from teosinte assemblies of *other* plants. Where our donors keep more of the B73 k-mers than those assemblies, a DONOR window reads too high and looks like SEG. Two things that do **not** fix it:
- **A minimum tract length** (RTIGER's rigidity R): the false SEG islands are 1–10+ Mb long, as long as real tracts, so a length rule would fix only a few of them and would delete the 30% of our tracts that are shorter than 3 Mb.
- **Step 20b** (family levels): it learns levels only from pools *called* DONOR, so it never sees the windows that are wrongly called SEG.

## What changes

`20c_pool_levels.py` learns each window's DONOR level from **all** pools of a taxon at once, without using any calls. This is the idea behind RTIGER's shared emission fit (one set of emission parameters learned from all samples), applied to the parameter that is wrong here. In each window, the scores r of all pools of the taxon are a mixture of three groups that share the window's level d:

| Group | Score r |
|---|---|
| B73 | ≈ 1 |
| SEG | d + (1 − d)·g(f), f = k/12, k ~ Binomial(12, ½) (as in step 20) |
| DONOR | ≈ d (times the pool's fitted scale) |

d is chosen by maximum likelihood on a grid (0.02–0.95), with the three group sizes of the window fitted by EM at each d. The estimate is shrunk towards the assembly level by how much evidence the window has (level = (n·d̂ + K·d_asm)/(n + K), n = expected number of pools carrying teosinte there, K = 5), so windows with few such pools keep the assembly level.

`20n_kmer_hmm.py` is a **copy** of step 20 with one option added, `--level_dir`, which reads these levels instead of the assembly mean. Everything else (scale grid, d_spread, SEG emission, transitions, rho grid, decoding) is the same as the current run.

**Test on simulated pools** (200 pools, a 40-window hotspot where the true level is 0.18 above the assembly level): the fitted level is within 0.02 of the truth in the hotspot (the assembly level is 0.18 off), and unchanged elsewhere (median difference 0.000). Exactly the 40 hotspot windows were raised.

## Run (from this folder, on the cluster)

```bash
cd new_model && mkdir -p logs
# 1. levels per taxon (any env with numpy, e.g. machine_learning)
sbatch --array=1-5 20c_pool_levels.sh
# 2. HMM with the new levels (env with numpy + scipy); same array size as step 20
N=$(awk -F'\t' 'NR>1 && ($2=="control"||$2=="donor") && $3!="Zm"' <D11>/pools.tsv | wc -l)
sbatch --array=1-$(( (N + 29) / 30 )) 20n_run_kmer_hmm.sh
# 3. pick rho and write the calls
sbatch 21n_kmer_collect.sh
# 4. targeted-sequencing check of the new calls (scripts of final_script, unchanged)
D21=/share/maize/ntanduk/kmer/introgression_spread; DNEW=/share/maize/ntanduk/kmer/new_model
python ../final_script/34c_targseq_compare.py --vcf_dir $D21/targseq/vcf --tracts $DNEW/introgression/tracts_all.tsv.gz \
  --pool_summary $DNEW/introgression/pool_summary.tsv --names ../data/phenotypes/raw/name_change_final.txt --out $DNEW/targseq
python ../final_script/34d_targseq_tracts.py --sites $DNEW/targseq/sites_used.tsv.gz \
  --tracts $DNEW/introgression/tracts_all.tsv.gz --out $DNEW/targseq_tracts
# 5. old vs new
python compare_old_new.py --old $D21/tracts_all.tsv.gz --new $DNEW/introgression/tracts_all.tsv.gz \
  --pool_summary $D21/pool_summary.tsv \
  --targseq_old ../data/introgression_spread/targseq/tracts_targseq.tsv --targseq_new $DNEW/targseq_tracts/tracts_targseq.tsv
```

**What success looks like:** SEG tracts between two DONOR tracts stop being BB in the targeted sequencing (47% now), the hotspot counts drop, and the rest of the validation table (B73 → AA, DONOR → BB, SEG → AB) does not get worse. Only then would the new levels replace the current ones.

## Result of the taxon-level fit (20c → 20n → 21n), 19 Nov 2026

| Check | Assembly levels (current) | Taxon levels (20c) |
|---|---|---|
| Genome share DONOR / SEG | 4.94% / 6.22% | 4.88% / 6.22% |
| SEG between two DONOR tracts (% of genome) | 0.38% | 0.32% |
| Hotspots: chr6 21–23 / chr1 297 / chr7 54–65 / chr3 75–78 Mb (pools) | 48 / 28 / 28 / 19 | 11 / 14 / 23 / 16 |
| Targeted seq.: DONOR tracts that are BB | 90.7% | 93.7% |
| Targeted seq.: SEG tracts that are AA / AB / BB | 17.9 / 74.3 / 7.8% | 6.2 / 81.8 / 12.1% |
| Targeted seq.: B73 tracts that are AA | 97.1% | 96.2% |
| Targeted seq.: SEG between two DONOR tracts that are BB | 47% (38 tracts) | 71% (38 tracts) |

Better DONOR calls, fewer SEG calls where the plant has no teosinte, and two hotspots mostly gone; but the
SEG-between-DONOR error stays (the remaining ones are BB even more often), and chr7 and chr3 barely move. So the
error is not one level per taxon: within a taxon, families carry different teosinte haplotypes. Next: one level
per accession. The genome shares hardly move, so the DONOR deficit (main README 2.8.1) does not depend on how the
DONOR level is set.

## Per accession (20d → 20na → 21na)

`20d_accession_levels.py` fits the same three-group model on the pools of **one accession** (84 accessions; median
11 pools; in a window often 1–3 carry teosinte). To stay stable with so few pools, the group sizes get 5 pseudo-pools
at the taxon-wide mix, and the level is shrunk towards the accession's taxon level from 20c:
level = (n·d̂ + K·d_taxon)/(n + K), K = 2. Accessions with fewer than 4 pools keep the taxon level.

Simulation (two accessions of 12 pools; one has a hotspot 0.18 above the taxon level): the error in the hotspot falls
from 0.18 to 0.07 with K = 2 (0.05 with K = 1, 0.09 with K = 3), the other accession is unchanged, and 90% of the
other windows stay within ±0.047 (K = 1: ±0.058; K = 3: ±0.039). So the per-accession fit corrects about 60% of an
accession-specific shift and adds little noise.

`20n_kmer_hmm.py` takes `--acc_level_dir` (with `--meta`): a pool uses its accession's level if there is one, else
the taxon level. Outputs go to `${DNEW}/acc/` so the taxon-level run is kept.

```bash
cd new_model && conda activate /rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk/machine_learning
sbatch --array=1-5 20d_accession_levels.sh                       # accession levels (needs 20c's levels)
D11=$(bash -c 'source ../final_script/config.sh; echo $D11')
N=$(awk -F'\t' 'NR>1 && ($2=="control"||$2=="donor") && $3!="Zm"' $D11/pools.tsv | wc -l)
J=$(sbatch --parsable --array=1-$(( (N + 29) / 30 )) 20na_run_kmer_hmm.sh)   # after 20d has finished
sbatch --dependency=afterok:$J 21na_kmer_collect.sh
# then the checks, with $DNEW/acc in place of $DNEW:
D21=/share/maize/ntanduk/kmer/introgression_spread; DA=/share/maize/ntanduk/kmer/new_model/acc
mkdir -p $DA/targseq $DA/targseq_tracts
python ../final_script/34c_targseq_compare.py --vcf_dir $D21/targseq/vcf --tracts $DA/introgression/tracts_all.tsv.gz \
  --pool_summary $DA/introgression/pool_summary.tsv --names ../data/phenotypes/raw/name_change_final.txt --out $DA/targseq > $DA/targseq/log.txt
python ../final_script/34d_targseq_tracts.py --sites $DA/targseq/sites_used.tsv.gz --tracts $DA/introgression/tracts_all.tsv.gz \
  --out $DA/targseq_tracts > $DA/targseq_tracts/log.txt
python compare_old_new.py --old $D21/tracts_all.tsv.gz --new $DA/introgression/tracts_all.tsv.gz --pool_summary $D21/pool_summary.tsv \
  --targseq_old ../data/introgression_spread/targseq/tracts_targseq.tsv --targseq_new $DA/targseq_tracts/tracts_targseq.tsv | tee $DA/compare_old_new.txt
```

Files: `20c_pool_levels.py` / `.sh` (taxon levels), `20d_accession_levels.py` / `.sh` (accession levels), `20n_kmer_hmm.py` (HMM, copy of step 20 with `--level_dir` / `--acc_level_dir`), `20n_run_kmer_hmm.sh` / `20na_run_kmer_hmm.sh` (taxon / accession runs), `21n_kmer_collect.sh` / `21na_kmer_collect.sh` (calls), `compare_old_new.py` (old vs new).
