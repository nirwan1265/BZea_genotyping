# final_script — ANGSD site discovery and per-pool read counts

This is a clean, numbered rewrite of the BZea genotyping steps. Each sequenced sample is a **pool of several BC2S3 plants from one plot (family)**, so the pipeline keeps **read counts per pool** instead of diploid genotype likelihoods. Every SNP is polarized to the **B73 reference base**.

- **Steps 01–06** document how SNP discovery was run. They have already been run, and their output lives in `.../angsd_full_snps/4_per_chr_renamed`.
- **Steps 07–10** are new. They produce the read-count data used for introgression calling.

All paths are set in `config.sh`, and every script reads them from there.

## Run

All jobs are Slurm (`sbatch`). Submit them from inside `final_script/`. Slurm needs the `logs/` folder to exist before a job starts, or the job fails without writing a log.

```bash
cd final_script
mkdir -p logs

# ---- read counts (run these now) ----
bash submit_counts.sh test     # 07 (10 tasks) + part 1 of 08; check its log and sacct
bash submit_counts.sh rest     # 08 parts 2-200, then 09 and 10, chained with --dependency=afterok

# if any 08 part was killed: list them and get the resubmit command
bash 08_check_unfinished.sh    # prints e.g.: sbatch --array=37,112 08_angsd_counts_by_part.sh

# the same steps by hand:
sbatch 07_make_site_tables.sh                                    # array 1-10, one per chromosome
sbatch --dependency=afterok:<07 jobid> 08_angsd_counts_by_part.sh   # array 1-200 (10 chr x 20 parts)
sbatch --dependency=afterok:<08 jobid> 09_concat_counts_by_chr.sh   # array 1-10
sbatch --dependency=afterok:<09 jobid> 10_counts_qc.sh

# ---- site filters + teosinte panel (after 10) ----
sbatch 11_filter_sites.sh                       # array 1-10
sbatch --dependency=afterok:<11 jobid> 12_signal_check.sh

# ---- introgression calling (after 11) ----
sbatch 13_extract_pool_gl.sh                    # array 1-10
sbatch --dependency=afterok:<13 jobid> 14_run_hmm.sh      # array 1-143 (10 pools each)
sbatch --dependency=afterok:<14 jobid> 15_collect_hmm.sh

# ---- site discovery (ALREADY RUN; only needed to redo it) ----
sbatch 01_angsd_discovery_by_part.sh            # array 1-200
bash 02_find_failed_parts.sh                    # prints the sbatch line for 03
sbatch --array=1-N%60 03_angsd_rerun_failed_split.sh
sbatch 04_merge_subparts.sh
sbatch 05_concat_parts_by_chr.sh                # array 1-10
sbatch 06_rename_samples.sh                     # array 1-10
# then set DISCOVERY_BCF_DIR=${D06} in config.sh
```

Options given on the `sbatch` command line override the `#SBATCH` lines. For example, `sbatch --array=5 --mem=64G 08_angsd_counts_by_part.sh` reruns part 5 with more memory.

Every job skips work that already finished (a `.done` marker), so resubmitting the same command is safe. To test one part interactively, pass the index: `bash 08_angsd_counts_by_part.sh 1`.

**Resources** (in each script's `#SBATCH` lines):

| Steps | CPUs | Memory | Time |
|---|---|---|---|
| ANGSD (01, 03, 08) | 4 | 16 / 48 / 32 G | 3 days 11 h |
| Others | 1 | 8–20 G | 2 h – 10 days |

No partition is set. If your account needs one, uncomment `##SBATCH --partition=...` in each script.

## Steps

| # | Script | Jobs | What it does | Main output |
|---|--------|------|--------------|-------------|
| 01 | `01_angsd_discovery_by_part.sh` | 200 | ANGSD SNP discovery + GL per part (`-GL 1 -doMajorMinor 1 -doMaf 1 -doBcf 1`) | `01_angsd_parts/*.bcf` |
| 02 | `02_find_failed_parts.sh` | login node | Finds parts without `[ALL done]` and splits each into 2 sub-intervals | `03_angsd_rerun/failed_subintervals.tsv` |
| 03 | `03_angsd_rerun_failed_split.sh` | N | Reruns those sub-intervals with 48 GB | `03_angsd_rerun/*.bcf` |
| 04 | `04_merge_subparts.sh` | 1 | Joins sub-intervals back into whole parts | `01_angsd_parts/*.bcf` |
| 05 | `05_concat_parts_by_chr.sh` | 10 | Joins 20 parts per chromosome; stops if any part is missing | `05_per_chr/BZea_chrN.full.bcf` |
| 06 | `06_rename_samples.sh` | 10 | BAM paths → `PN*_SID*`, after checking sample order | `06_per_chr_renamed/BZea_chrN.full.renamed.bcf` |
| 07 | `07_make_site_tables.sh` + `07_polarize_to_b73.py` | 10 | All biallelic ANGSD SNPs, re-polarized to B73; ANGSD `-sites` files; B73 control list; VCF header | `07_sites/` |
| 08 | `08_angsd_counts_by_part.sh` + `08_counts_to_ad.awk` | 200 | `angsd -doCounts 1 -dumpCounts 4 -sites`, with the same read filters as discovery → B73/ALT counts per pool | `08_counts_parts/*.AD.bcf` |
| 08 | `08_check_unfinished.sh` | login node | Lists unfinished 08 parts and prints the resubmit command | — |
| 09 | `09_concat_counts_by_chr.sh` | 10 | Joins parts per chromosome; builds site and pool summaries | `09_counts_by_chr/` |
| 10 | `10_counts_qc.sh` + `10_counts_qc.R` | 1 | Per-pool QC: reads, coverage, ALT fraction against the B73 controls, error proxy | `10_qc/` |
| 11 | `11_filter_sites.sh` (bcftools + awk) | 10 | Site filters (clean B73 controls, depth, third base, panel allele mismatch); splits sites into in / not in the 227-plant teosinte panel; reads per pool on each set | `11_sites_filtered/` |
| 12 | `12_signal_check.sh` + `12_signal_check.R` | 1 | Donor signal per taxon on teosinte-panel sites vs other sites: are panel sites alone enough? | `12_signal_check/` |
| 13 | `13_extract_pool_gl.sh` (bcftools + awk) | 10 | Per-pool GL files for the HMM: ANGSD GLs from the discovery BCF at PASS teosinte-panel sites, flipped to B73 orientation | `13_pool_gl/` |
| 14 | `14_run_hmm.sh` + `14_hmm_pooled.py` | 143 | Nirwan's BC2S3 HMM adapted to pools: GL emissions for B73 / SEG / DONOR, NAM genetic-map transitions, per-pool κ by EM | `14_hmm/` |
| 15 | `15_collect_hmm.sh` + `.R` | 1 | Collects all pools: summary table and all tracts | `15_introgression/` |

## Final outputs (`09_counts_by_chr/`)

**`BZea_chrN.AD.bcf`** holds the read counts.
- `REF` is the B73 base and `ALT` is the other allele.
- `FORMAT/AD` = (B73 reads, ALT reads) per pool. A pool with no reads at a site is `0,0`.
- INFO fields:
  - `DP`: all reads at the site.
  - `OTH`: reads carrying a third base.
  - `NCOV`: pools with any B73 or ALT read.
  - `NALT`: pools with any ALT read.
  - `SWAP`: ANGSD's major allele was not B73.

```bash
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' BZea_chr1.AD.bcf | head
bcftools view -s PN9_SID857,PN9_SID854 BZea_chr1.AD.bcf | ...     # subset pools, e.g. by donor accession
```

**`BZea_chrN.sitesum.tsv.gz`** has one row per site. Columns: `CHROM POS B73 ALT DP OTH REF_READS ALT_READS NCOV NALT CTRL_REF CTRL_ALT CTRL_NCOV`. The `CTRL_*` columns are reads from the 13 B73 control pools; an ALT read there flags a reference error, paralog or seed-lot difference. These are the inputs for step-2 site filtering.

**QC decisions:** see [`QC_NOTES.md`](QC_NOTES.md) and the lists in `sample_lists/`.

**`BZea_chrN.sampleqc.tsv`** has per-pool totals, and **`BZea_chrN.siteqc.tsv`** a one-line chromosome summary.

QC in `10_qc/`:
- `sample_qc.tsv` (one row per pool)
- `taxon_summary.tsv`
- `site_qc.tsv`
- `qc_plots.png`

The B73 controls give the background ALT fraction; donor pools should sit above it.

## Notes on decisions and problems found in the earlier pipeline

1. **ANGSD's BCF alleles are major/minor, not REF/ALT.** ANGSD writes `major,minor` (`abcWriteBcf.cpp`), and its GLs are ordered major/major, major/minor, minor/minor. Wherever the major allele is not the B73 base, the old GL "RR" meant *donor* homozygous. Step 07 re-polarizes every SNP to the B73 FASTA base and logs how many were flipped (`07_sites/BZea_chrN.polarize.log`).

2. **The old `0.6_filter_bcf.sh` is not used.** Its population-wide `AF 0.05–0.95` cut removes donor-specific SNPs. Example: a SNP carried only by parviglumis donors has expected AF ≈ 394/1510 × 0.125 ≈ 0.03 across all pools. Its `F_MISSING` and `AF` were also computed from fields ANGSD fills differently from bcftools. Counts are therefore taken at **all** biallelic ANGSD SNPs. Site filtering moves to step 2 and uses the B73 controls, depth and donor groups.

3. **Counting uses the same read filters as discovery**, set in `READ_FILTERS` in `config.sh`. It drops `-minInd` and `-setMinDepthInd`, so every read that passed the filters is counted.

4. **The killed parts were rerun with stricter settings** (`SNP_pval 1e-10`, `minInd 20`, `setMinDepthInd 4`; the main run used `1e-8`, `5`, `2`). Those parts have a sparser, stricter SNP set than the rest of the genome. `config.sh` records the values that were actually used. Set the `RERUN_*` values equal to the main ones if discovery is ever redone.

5. **Metadata issues** (`data/name_change/name_change_final.txt`):
   - The file has Windows (CRLF) line endings. Step 07 and the QC strip them.
   - `PN13_SID1225` is listed twice, as `Zd.0020_P3_P1_P2.5.1.1` and as `PurpleCheck`. QC labels it `conflicting_metadata`.
   - The metadata has 13 B73 controls, but `scripts/HMM_introgression/b73.list` has 12 (it is missing `PN1_SID85`). Step 07 takes the controls from the metadata.

6. **PoolParty2 is not used.** It starts from FASTQ and realigns, and it sums barcoded individuals into population pools. It also requires ≥10 reads per pool, which at ~0.8x per pool removes nearly every site. Its outputs (allele frequencies, F<sub>ST</sub>, CMH, FLK) compare populations rather than giving per-line ancestry. Its useful part, per-sample coverage QC, is covered by step 10.

## Old → new scripts

| Old (`scripts/HMM_introgression/`) | New |
|---|---|
| `0.0_angsd_genotyping.sh` (manual `CHR=` edit per chromosome) | `01` (one 200-job array) |
| `0.1_find_failed_runs.sh` (grep LSF logs) | `02` (checks ANGSD's `[ALL done]`) |
| `0.2_run_failed_parts.sh` | `03` |
| `0.3_combine_sub_parts.sh` | `04` |
| `0.4_combine_parts_chrom.sh` (manual chr loop, globbed parts) | `05` (computed part names, gap check) |
| `0.5_rename_headers_files.sh` (hand-made rename list) | `06` (builds the list and checks the order) |
| `0.6_filter_bcf.sh`, `1_extract_GL_emissions.sh`, `1.2_filter_GL_emissions.sh` | replaced by `07`–`10` |

## Size and time

These figures come from a test of the awk converter at 1,500 pools with ~0.8 reads per pool per site:
- **Conversion speed:** about 1.5 ms per site with mawk and 2.2 ms with gawk, on top of the ANGSD run.
- **Storage:** the AD BCF takes about 0.5 kB per site, roughly 4.6 GB per 10 M sites.
- **Raw counts:** ANGSD's A/C/G/T table is deleted after conversion unless `KEEP_RAW_COUNTS=1`.

**Open files.** ANGSD keeps every BAM open at once, about 1,535 files. Steps 01, 03 and 08 raise the open-file limit to the node's hard limit. If that is still too low, 08 splits the BAM list into batches that fit, runs ANGSD on each batch, and merges the counts back in the original column order (`08_merge_batches.py`). The result is identical to a single run. To force smaller batches, set `BAMS_PER_RUN=500` in `config.sh`.

Memory is set to 32 GB per part in 08. If parts get killed, resubmit them with `sbatch --array=<ids> --mem=64G 08_angsd_counts_by_part.sh`, or set `NPARTS=40` in `config.sh` and change `--array` in 01 and 08 to `1-400%40`.

## Requirements

- angsd
- bcftools and bgzip (htslib)
- python3 (standard library only)
- Rscript (base R only)
- awk (mawk or gawk)

Load them where `config.sh` says.

## Steps 11–12: site filters, teosinte panel, signal check

These steps use only `bcftools`, `awk` and base R. They work on the pools we trust: the 10 clean B73 controls and the donor pools not listed in `sample_lists/exclude_all.txt`. The reasons for each exclusion are in [`QC_NOTES.md`](QC_NOTES.md).

**Site filters (step 11).** A site is removed if any of these is true:
- **Control ALT reads:** the clean B73 controls have an ALT read there (`bcftools view -S controls | bcftools query -i 'SUM(FMT/AD[*:1])>0'`). A real B73 plant cannot carry ALT, so ALT reads point to a repeat, a paralog or a reference error.
- **Depth:** depth is more than 2× or less than 0.25× the chromosome median (`MAX_DP_MULT`, `MIN_DP_MULT`).
- **Third base:** more than 5% of the site's reads show a base that is neither B73 nor ALT (`MAX_OTH_FRAC`).
- **Panel mismatch:** the teosinte panel has the position, but with different alleles (`bcftools isec -c all` minus `-c none`).

**Outputs per chromosome:**
- `BZea_<chr>.PASS.bcf`: all sites that passed, trusted pools only.
- `BZea_<chr>.PASS.panel.bcf`: passed sites that the teosinte panel also has, with the same alleles (`bcftools isec -c none`). Each site gets `INFO/PANEL_AC` and `PANEL_AN`: how many of the 454 teosinte chromosomes carry ALT.
- `BZea_<chr>.PASS.nonpanel.bcf`: passed sites that are not in the panel.
- `BZea_<chr>.pool_reads.tsv`: B73 and ALT reads per pool on panel and non-panel sites.
- `BZea_<chr>.ctrl_background.tsv`: the B73-control ALT rate, measured fairly. A control's ALT count on sites it helped filter is 0 by construction, so the controls are split in two halves, and each half is counted on sites filtered using only the other half.
- `BZea_<chr>.summary.tsv`: how many sites each filter removed.

**Signal check (step 12).** For a typical pool of each taxon, step 12 reports, per 1 Mb:
- reads on panel sites
- expected ALT reads in a B73 region
- expected ALT reads in a fully donor region
- `z_fixed`, how clearly those two separate

`z_fixed ≥ 10` on panel sites means the panel sites alone are enough for that taxon. Where it isn't, that taxon also gets its non-panel PASS sites in the introgression caller.

## Steps 13–15: introgression calling (Nirwan's HMM, adapted to pools)

The HMM is `old_scripts/HMM_introgression/hmm_viterbi_bc2s3.py`, adapted in `14_hmm_pooled.py`. These parts are unchanged from the original:
- per-marker HMM
- NAM genetic map (`genetic_maps/`): Haldane recombination fraction, scaled by `rho`
- sticky/stationary transition matrix
- forward–backward posteriors
- posterior + rigidity (hysteresis) or Viterbi decoding, with optional min-run cleanup

**Input (step 13).** The ANGSD GLs come from the existing discovery BCFs, with no ANGSD rerun. Only the PASS teosinte-panel sites from step 11 are used. ANGSD stores alleles as major,minor, so at sites where the major allele is not the B73 base, the three GLs are flipped. `GL = 0,0,0` (no reads in that pool) is left out.

**What changed in the HMM:**

1. **States are for a pool (one BC2S3 family), not one plant.** The GLs are read as the likelihood of the reads if the pool's ALT allele fraction is 0, ½ or 1. The likelihood at any other fraction comes from the curve through those three points:

   L(p) = L0 (2p−1)(p−1) + L1·4p(1−p) + L2·p(2p−1)

   This is exact for sites with 1–2 reads, which is almost all sites at ~0.8x.

   | State | Meaning | Emission |
   |---|---|---|
   | **B73** | family fixed B73 | L(0) |
   | **SEG** | family segregating | c · mean of L(¼), L(½), L(¾) + (1−c) · L(0) |
   | **DONOR** | family fixed donor | c · L(1) + (1−c) · L(0) |

   `c` is the chance the donor carries the ALT allele at the site: κ × the teosinte panel's AC/AN. Sequencing error is already inside the GLs. The old GL adjustments (`eta_hh_from_rh`, `eta_rr_from_rh`, `rh_penalty`) are dropped.

2. **κ per pool, learned by EM.** κ is how often this pool's donor carries the panel's ALT alleles. It is near 1 for donors like the panel (parviglumis, mexicana) and lower for donors the panel represents poorly (e.g. diploperennis). Without κ, their fixed-donor regions look segregating.

   At one read per site, "SEG at κ" and "DONOR at κ/2" look alike. So EM runs from several starts (κ = 1, 0.5, 0.25) and keeps the start with the best log-likelihood; sites with 2+ reads tell the two apart.

3. **Each chromosome is its own chain.** The original ran chr1–chr10 as one chain.

4. **`thetas[t]` is the recombination fraction between markers t−1 and t.** The original was shifted by one marker.

5. **Tracts cover each chromosome end to end.** Boundaries sit at the midpoint between markers, so tract lengths add up to the chromosome length.

**Outputs:**
- In `14_hmm/`, per pool:
  - `<SAMPLE>.statepath.tsv.gz`: per marker: GLs, panel AC/AN, state, P(B73), P(SEG), P(DONOR), theta.
  - `<SAMPLE>.tracts.bed.gz`
  - `<SAMPLE>.summary.tsv`: fractions B73/SEG/DONOR, donor share, number of tracts, κ.
- In `15_introgression/`:
  - `pool_summary.tsv`
  - `tracts_all.tsv.gz`

**Simulation check:** 54 pools (Zx, Zv, and Zd with a donor that carries panel ALT alleles at 0.35× the panel frequency), 100 Mb, ~0.6 reads per pool per site, ANGSD-style GLs, 5% of sites with a swapped major allele.
- Marker accuracy: B73 controls 100%, Zx 99.5%, Zv 99.6%, Zd 99.3%.
- κ came out at 1.0 for Zx and Zv, and 0.36 for Zd (true value 0.35).
- Donor share per pool vs truth: correlation 1.00, mean error 0.0015.

**Run time:** about 9 s per pool in the simulation (~54k markers). Real pools have ~1 M markers, so expect roughly 3–4 min per pool and about 35 min per array task of 10 pools.

**BC1 sisters (later):** a 10x sister from the same F1 shows exactly which ALT alleles that F1's donor carries. It can replace `c = κ × panel AC/AN` for that family.
