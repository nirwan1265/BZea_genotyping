# fix1: count only informative k-mers, where it helps (full build)

The pilot (`../pilot_fix1/`) showed that at *parviglumis* and *mexicana* chr7 54–65 Mb, counting only the **informative** k-mers (B73 single-copy k-mers missing from all of the taxon's assemblies) raises the DONOR-vs-SEG contrast from 0.5–1.3 to 1.9–3.2. With these counts the SEG islands there read as DONOR, as the targeted sequencing says. Away from the hotspots the contrast hardly changed, and at *parviglumis* chr3 75–78 Mb nothing changed. So this build switches to the new counts **only in the windows where they raise the contrast**, and keeps the old counts, assembly levels and everything else of the current run elsewhere. Nothing in `final_script/` or the current run is changed. All outputs go to `${KMER_DIR}/new_model/fix1` (`FX`).

| Step | Script | What it does |
|---|---|---|
| 1 | `fx1_1_sets.sh` (array 1–5: Zd Zh Zl Zv Zx) | informative k-mers per taxon (Zv, Zx copied from the pilot) and their count per window in B73 (denominator) |
| 2 | `fx1_2_pool_windows.sh` (2 pools per task, `%100`) | every pool, as step 17: count reads (KMC), keep its taxon's informative k-mers, windows on chr1–10. The 10 B73 controls are counted for all 5 taxa |
| 3 | `fx1_3_choose.sh` → `fx1_3_choose.py` (array 1–5) | per taxon and window: old and new scores of every pool (step 20's normalisation), noise from the controls and pools called B73, DONOR level from pools called DONOR (tracts ≥ 2 Mb), contrast old vs new. A window switches if B73 has ≥ 200 informative k-mers there and the new contrast is ≥ 1.2× the old one. For switched windows, the DONOR level and its spread between our donors come from our own pools (the assemblies have none of these k-mers by construction) |
| 4 | `fx1_4_run_hmm.sh` → `fx1_kmer_hmm.py` | the HMM (copy of `../20n_kmer_hmm.py` plus `--fix_dir`): in switched windows a pool and the controls use the new counts, level and spread; elsewhere it is exactly step 20 |
| 5 | `fx1_5_collect.sh` | rho and calls (`final_script/21_kmer_collect.py`, unchanged) → `${FX}/introgression/` |
| 6 | `fx1_6_checks.sh` | targeted sequencing (`34c`, `34d`, unchanged) and old vs new (`../compare_old_new.py`) → `${FX}/compare_old_new.txt` |

One difference from step 20, in switched windows only: step 20 drops a window when the controls see less than 20% of the genome-median count. Informative k-mers are few where the teosinte keeps most B73 k-mers, which is exactly at the hotspots, so that rule would drop the windows the fix is for. In switched windows the rule compares the **share** of the window's k-mers that the controls see with the typical share instead.

**Test on simulated pools** (40 *parviglumis*-like pools, 10 controls, chr3 + chr7):
- **The setup:** a chr7 hotspot where the assemblies keep 80% of the B73 k-mers and our donors keep 90%, with overdispersed counting noise. The current calls have a false SEG island there in half of the DONOR pools.
- **The old model** calls the hotspot SEG in all 10 DONOR pools tested (contrast 1.3).
- **This build** switches 118 of the 120 hotspot windows (contrast 4.7) and calls all 10 DONOR there.
- **No side effects:** SEG and B73 pools are unchanged, and the windows called correctly genome-wide rise from 97.1% to 100% in the DONOR pools.

## Run (one command, on the login node; each job loads its own env)

```bash
cd /share/maize/ntanduk/BZea_genotyping && git pull
cd new_model/fix1 && bash fx1_submit.sh
```

It prints the number of pools and the job ids (also in `logs/fx1_jobs.txt`). Expected time: step 1 is 1–3 h. Step 2 is about 1,400 pools at 30–60 min each, 100 at a time, so roughly 6–10 h. Steps 3–6 take about 1–2 h together. Step 2 needs about 25 GB of temporary disk per running pool in `${FX}/tmp`, so about 2.5 TB with 100 running.

Watch: `squeue -u $USER`; failed pools: `grep -l FAILED logs/fx1_2_pools.*.out`; resubmitting is safe (finished pieces are skipped): `bash fx1_submit.sh`.

**In the morning:** `cat ${FX}/compare_old_new.txt` (old vs new: genome shares, SEG-between-DONOR, hotspot counts, targeted-sequencing table) and `cat ${FX}/choice/summary_*.tsv` (windows switched per taxon, chromosome and hotspot; contrast before and after).

**What success looks like:** sandwiched SEG tracts that are BB in the targeted sequencing drop from 47%, the chr7 hotspot count drops, and B73 → AA, DONOR → BB and SEG → AB do not get worse. *parviglumis* chr3 75–78 Mb is not expected to change (pilot: the islands there read between DONOR and SEG in every count).
