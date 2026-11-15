# final_script — BZea introgression calling

Each sequenced sample is a **pool of 6 BC2S3 plants from one plot (family)**, sequenced at about 0.8x, with no parental sequencing. B73 is the recurrent parent. The donors are teosinte: mexicana (Zx), parviglumis (Zv), diploperennis (Zd), luxurians (Zl) and huehuetenangensis (Zh). Landraces (Zm) are set aside for now.

**Introgressions are called from k-mers (steps 16–21).** The pool's reads are compared with B73 directly, without mapping. An HMM then runs along the genome on 100-kb windows. The SNP route (steps 01–15) is kept as a cross-check, and its outputs (sites, read counts, QC) are still used. All paths and settings are in `config.sh`, and every script reads them from there.

Contents:
1. [Run order](#1-run-order)
2. [Why k-mers, not mapped SNPs](#2-why-k-mers-not-mapped-snps)
3. [k-mer pipeline, steps 16–19](#3-k-mer-pipeline-steps-1619)
4. [The HMM, steps 20–21: model and every parameter](#4-the-hmm-steps-2021-model-and-every-parameter)
5. [Results: all 1,395 pools](#5-results-all-1395-pools)
6. [Validation](#6-validation)
7. [Open issues and next steps](#7-open-issues-and-next-steps)
8. [SNP route, steps 01–15 (cross-check)](#8-snp-route-steps-0115-cross-check)
9. [Requirements](#9-requirements)

---

## 1. Run order

All jobs are Slurm. Submit them from `final_script/` (`mkdir -p logs` first). Activate the k-mer conda env (kmc, kmc_tools, kcftools, samtools) before submitting steps 16–21; jobs inherit it. Every step skips work that is already done, so resubmitting is safe.

```bash
cd final_script; mkdir -p logs; source config.sh
# ---- k-mer introgression calling (main) ----
sbatch 16_kmer_b73_unique.sh                                   # once: B73 single-copy 31-mers
N=$(awk -F'\t' 'NR>1 && ($2=="control"||$2=="donor")' ${D11}/pools.tsv | wc -l)
sbatch --array=1-$(( (N + 3) / 4 ))%50 17_kmer_windows.sh       # k-mer windows per pool (4 pools per task)
bash 18_kmer_check.sh                                          # missing pools, k-mers per pool
sbatch 19_kmer_donor_reference.sh                              # expected DONOR level from teosinte assemblies
N=$(awk -F'\t' 'NR>1 && ($2=="control"||$2=="donor") && $3!="Zm"' ${D11}/pools.tsv | wc -l)
J=$(sbatch --parsable --array=1-$(( (N + 29) / 30 )) 20_run_kmer_hmm.sh)   # HMM, 30 pools per task
J2=$(sbatch --parsable --dependency=afterok:$J 21_kmer_collect.sh)   # pick rho, calls, tracts, figures
sbatch --dependency=afterok:$J2 22_kmer_stats.sh                # observed vs expected, family bootstrap
```

Steps 01–15 are in [section 8](#8-snp-route-steps-0115-cross-check).

---

## 2. Why k-mers, not mapped SNPs

The SNP HMM (step 14) gave too little donor, and more SEG than DONOR. Read counts in one known Zd DONOR tract show why: PN13_SID1199, chr2:22.1–25.3 Mb, compared with a B73 region of the same size and normalized with a B73 control.

| What is kept in the donor tract | Share of expected |
|---|---|
| Reads that map there at all | 34% |
| + MAPQ ≥ 30 | 30% |
| Bases, MAPQ 30, no `-C 50`, no BAQ | 25% |
| Bases with `-C 50` + BAQ (the ANGSD settings used) | **11%** |
| Same filters in a B73 region | ~96% |

Three problems follow, and they hit DONOR and SEG (het) regions hardest:

1. **Donor reads are lost (reference bias).** About 2/3 of donor reads don't map to the region; much of that is donor sequence with no counterpart in B73. `-C 50` and BAQ then remove another ~60% of the donor bases that do map: they lower the score of reads with mismatches, and a read carrying a teosinte allele has extra mismatches by definition. ALT reads are lost first, so **DONOR (BB) looks mixed and is called SEG, and SEG loses its ALT reads and is called B73.**
2. **Too few reads per site.** At 0.8x most SNPs have 0 or 1 read. One read shows one allele, so a het site looks homozygous, and a single error read looks like a het.
3. **Pools aren't diploid.** A pool's ALT fraction can be anything from 0 to 1 (e.g. 5/12 in a SEG region). Diploid genotype callers (0 / ½ / 1) mislabel these mixtures.

---

### 2.1 Prior work

Finding introgressions by counting reference k-mers in windows, without mapping, was developed in wheat: **IBSpy** (Uauy lab; Quiroz-Chávez 2023; einkorn introgressions in bread wheat, Ahmed et al. 2023) and **KCFtools** (Selvanayagam et al. 2026, *Bioinformatics*, doi:10.1093/bioinformatics/btag354), which we use for the window counts. k-mers are counted with KMC 3 (Kokot et al. 2017); DONOR levels come from PanAnd assemblies (Stitzer et al. 2025). What we added: low-coverage pools, the three-state pool model with per-window DONOR levels, and the HMM on the NAM map. Full references are in the [top-level README](../README.md).

### 2.2 How k-mers fix each problem

| Problem | How the k-mer approach avoids it |
|---|---|
| **1. Donor reads are lost** | No mapping. Raw reads are counted as k-mers and only the presence of B73 k-mers is checked, so donor reads can't be filtered away. |
| **2. Too few reads per site** | No per-site calls. Each 100-kb window adds up ~30,000 k-mers, so r is measured from thousands of positions at once. |
| **3. Pools aren't diploid** | The model is built for pools: DONOR = d, SEG = d + (1 − d)·g(f), averaged over Binomial(12, ½) donor fractions of a 6-plant family. |

**Telling DONOR from SEG.** It no longer depends on allele ratios at single sites, which reference bias distorts. It depends on the window's **level**: in a Zd example, DONOR ≈ 0.27 and SEG ≈ 0.69, with B73 = 1. The DONOR level is measured from teosinte genomes, not guessed.

**What we saw in our data:**
- **SEG matches the expectation:** 6.4% observed vs 6.25% expected (family bootstrap p = 0.26).
- **Plant-level het matches Mendel:** AB 3.2% vs 3.1%. So the k-mer calls don't inflate het, which is the classic failure.
- **Fewer broken tracts:** the SNP route broke tracts into many short SEG pieces (44–69% of its SEG tracts were under 1 Mb) and turned DONOR into SEG when its donor-allele model drifted. In two Zd pools it reported 17% donor where the k-mers find 3–7%.

---

## 3. k-mer pipeline, steps 16–19

**Idea.** Take a 100-kb window of B73. Count how many of its **single-copy** 31-mers (k-mers that occur exactly once in the whole B73 genome) appear in a pool's reads.
- **B73 regions:** the count is what a B73 pool of the same depth gives.
- **Donor regions:** every SNP, indel or missing sequence in the donor removes the B73 k-mers that overlap it, so the count drops.

Repeats are left out because they're present in every genome, donor or not. With repeats included, a fixed donor tract scored 0.91. With single-copy k-mers only, it scores 0.20.

| Step | Script | What it does | Output |
|---|---|---|---|
| 16 | `16_kmer_b73_unique.sh` | Counts all 31-mers in B73 v5 (KMC) and keeps those that occur once (`kmc_tools transform -ci1 -cx1`). That's 680 M of 817 M (83%). Also writes the chr1–10 FASTA. | `${KMER_DIR}/B73_uniq.kmc_*`, `B73_chr1-10.fa` |
| 17 | `17_kmer_windows.sh` | Per pool: (1) KMC count of the trimmed paired reads, `-ci1` (at 0.8x most k-mers are seen once); (2) keep k-mers single-copy in B73 (`kmc_tools simple … intersect`); (3) dump and recount with `kmc`, because KCFtools reads only the KMC2 format; (4) `kcftools getVariations` on 100-kb windows of chr1–10. | `${KMER_DIR}/windows/<SAMPLE>.tsv.gz` (`observed_kmers` per window) |
| 18 | `18_kmer_check.sh` | Pools without output, missing FASTQs, k-mers per pool | `${KMER_DIR}/kmer_pool_qc.tsv` |
| 19 | `19_kmer_donor_reference.sh` | Scores B73 and 8 teosinte genome assemblies (PanAnd) like pools. Per window, d = the share of B73 single-copy k-mers present in that taxon's genome: the level a fixed donor tract should show. | `${KMER_DIR}/donor_ref/<NAME>.tsv.gz` |

**Inputs:**
- Reads: `FASTQ_DIR/filtered_S<plate>/<SAMPLE>/<SAMPLE>_paired_R1/R2.fq.gz`
- Assemblies: `DONOR_ASSEMBLIES` in `config.sh`
- Outputs go to `KMER_DIR=/share/maize/ntanduk/kmer`.

**Disk:**
- Temporary: ~25 GB per running pool (17), deleted after each pool.
- Pool read databases: ~2–3 GB each, deleted unless `KMER_KEEP_KMC=1`.

**Expected DONOR level from the assemblies (step 19), whole genome:**

| Taxon | Assemblies | d (genome) | d per window, median (10–90%) |
|---|---|---|---|
| Zv | TIL01, TIL11 | 0.43, 0.44 | 0.38 (0.15–0.79) |
| Zx | TIL18, TIL25 | 0.38, 0.37 | 0.32 (0.13–0.67) |
| Zh | RIMHU001 | 0.34 | 0.30 (0.14–0.55) |
| Zd | Gigi, Momo | 0.26, 0.28 | 0.23 (0.13–0.42) |
| Zl | RIL003 | 0.26 | 0.22 (0.13–0.40) |

Zv shares the most single-copy k-mers with B73, and Zd and Zl the fewest, which fits the Zea phylogeny. In the 44-pool test, the level observed in pools' DONOR windows matched their taxon's assembly value. For example, Zh pools were at 0.22–0.32 against an assembly value of 0.26–0.31, with a window-by-window correlation of up to 0.8.

---

## 4. The HMM, steps 20–21: model and every parameter

**Code:**
- `22_kmer_stats.py` (step 22): observed vs expected tests, run by `22_kmer_stats.sh`.
- `20_kmer_hmm.py` (per pool), run by `20_run_kmer_hmm.sh` (30 pools per array task).
- `21_kmer_collect.py`, run by `21_kmer_collect.sh`.

The HMM itself is Nirwan's BC2S3 HMM. `20_kmer_hmm.py` imports these functions unchanged from `14_hmm_pooled.py`:
- `read_map` / `interp_cm`: NAM genetic map
- `haldane_theta`
- `build_A`: transitions
- `forward_backward`
- `decode_hysteresis`

**Only the emission is new.**

### 4.1 Observation per window

1. **Scale the pool:** N = the pool's single-copy k-mer count in the window ÷ the pool's genome-wide median.
2. **Build the reference:** a non-negative least-squares fit (`scipy.optimize.nnls`) of the B73 controls' scaled counts to this pool. A control is never part of its own reference. The fit is redone without windows below 1 − 2σ until it stops changing, so donor windows don't pull the reference down.
3. **Score:** r = N / reference, rescaled so the B73-like windows have median 1. **r = 1 means "looks like B73".**
4. **Skip:** windows with fewer than 2,000 B73 single-copy k-mers, or where the controls' median scaled count is below 0.2.

#### How r is calculated: a worked example

Take one 100-kb window. The numbers are rounded but realistic.

**Step 1: what B73 has.** B73's window contains **30,000** single-copy k-mers (step 16).

**Step 2: what a pure-B73 pool would find.** At 0.8x a pool doesn't see every k-mer, because many parts of the window aren't covered by any read. The 10 B73 control pools tell us how many a B73 pool of this depth finds in this window; say **13,000**. That's the expected count.

**Step 3: r = what our pool found ÷ the expected count.**

| Our pool in this window is… | k-mers found | r |
|---|---|---|
| B73 (all plants AA) | ~13,000 | 13,000 / 13,000 ≈ **1.0** |
| DONOR (all plants BB, Zd) | ~3,500 | 3,500 / 13,000 ≈ **0.27** |
| SEG (about half B73, half donor) | ~9,000 | ≈ **0.69** (below) |

**Where the DONOR level comes from (step 19).** We take a complete Zd genome assembly (Gigi, Momo) and count how many of the window's 30,000 B73 k-mers appear anywhere in it. Say 8,100: 8,100 / 30,000 = **d = 0.27**, so Zd keeps 27% of B73's k-mers here. A pool that is all Zd in this window can only find that 27%, so its expected r is 0.27.

d is measured like this for every window and every taxon. Genome-wide it's Zd 0.27, Zl 0.26, Zh 0.34, Zx 0.37 and Zv 0.43; Zv is closest to maize, so it keeps the most. In later rounds (step 20b), d comes from the family's own DONOR pools instead.

**Why SEG falls in between.** In a SEG region the pool is a mix of B73 and donor chromosomes. Split the 30,000 k-mers in two:
- **Shared with the donor (27%):** every chromosome carries them, so they're found as often as in a B73 pool (ratio 1.0).
- **Missing from the donor (73%):** only the B73 chromosomes carry them. With half the chromosomes B73, they get half the reads, so at this depth they're found about **57%** as often.

So r(SEG) = 0.27 × 1.0 + 0.73 × 0.57 ≈ **0.69**.

In general, r(SEG) = d + (1 − d) × g, where g = how often a k-mer is seen at the pool's B73 share of the reads ÷ how often at full depth. g comes from the standard coverage formula: the chance that a k-mer is seen at coverage c is 1 − e^(−c). The HMM doesn't assume exactly half donor. It averages over the possible donor fractions of a 6-plant pool (1/12 to 11/12), weighted by Binomial(12, ½), the 1:2:1 segregation summed over 6 plants.

**Then the HMM** compares each window's r with these three expected levels (B73 ≈ 1, SEG ≈ 0.6–0.7, DONOR ≈ d). It allows for each pool's noise and takes neighbouring windows into account, using the genetic map, to call the most likely state.

### 4.2 States and emissions

| State | Meaning in a pool of 6 plants | Emission for r |
|---|---|---|
| B73 | all plants B73/B73 | Normal(1, σ_w) |
| SEG | the family still segregates here: a mix of AA, AB and BB plants | mixture below |
| DONOR | all plants donor/donor (BB) | Normal(a·d_w, σ_w) |

- **d_w:** the expected DONOR level of the pool's taxon in that window: the mean of that taxon's assemblies (step 19).
- **a:** a per-pool scale that allows for the pool's donor plant differing from the assemblies. Fitted per pool.
- **SEG:** a family segregating 1:2:1 at a locus gives each plant 0, 1 or 2 donor copies, so the pool's donor copies ~ **Binomial(2n, ½)**, with n = 6 plants. The pool's donor fraction is f = k/12, k = 1–11, with binomial weights. k-mers the donor lacks come only from the B73 part of the pool, so SEG for fraction f gives Normal(d_w + (1 − d_w)·g(f), σ_w), where:
  - g(f) = (1 − e^(−(1−f)λ)) / (1 − e^(−λ)) (Lander–Waterman: the share of k-mers seen at coverage λ is 1 − e^(−λ))
  - λ = the pool's k-mer coverage, from the median share of single-copy k-mers it sees.
- **σ_w:** the pool's robust SD of r in B73 windows (1.4826 × MAD) × how much the controls disagree in that window relative to a typical window (clipped to 0.5–3×).
- **Outliers:** every state has a 2% chance that the window is junk (uniform on 0–2).

### 4.3 Transitions and decoding (unchanged from step 14)

- **Genetic distance:** cM interpolated at window midpoints from the NAM map (`genetic_maps/`), converted to θ with Haldane.
- **Switching:** probability s = ρ·θ (clamped to 1e-12–0.25). The new state is chosen in proportion to the prior π.
- **Chains:** each chromosome is its own chain.
- **Decoding:** forward–backward, then the most likely state per window.

### 4.4 ρ is fitted, R is not used

- **ρ:** step 20 runs every pool over a ρ grid and saves the log-likelihood. Step 21 sums it over donor pools that are not NOISY and uses the best ρ for every pool. This is the same as estimating the number of generations in local-ancestry HMMs (e.g. Ancestry_HMM).
  - Best: **ρ = 1.5**, from 1,292 pools; the 44-pool test also gave 1.5.
  - Alone, each taxon would pick ρ = 1 (Zd, Zh, Zl) or 2 (Zv, Zx).
- **R:** calls are made without rigidity. With R = 3 windows, the donor share differed by more than 0.5 percentage points in **0 of 1,385** pools, so R isn't needed once ρ is fitted. `DONOR_SHARE_R` in the summary keeps the R = 3 value for reference.

### 4.5 All parameters

| Parameter | Value | Where | Source / justification |
|---|---|---|---|
| k-mer size | 31 | `KMER_K` | standard (KMC, KCFtools) |
| Window | 100 kb | `KMER_WINDOW` | ~30,000 single-copy k-mers per window |
| Single-copy | count = 1 in B73 v5 | step 16 | removes repeats |
| KMC minimum count | 1 | steps 17, 19 | at 0.8x most k-mers are seen once |
| Reference fit | NNLS of B73 controls, refit without windows below 1 − 2σ | `--cut_sigma 2` | robust; no fixed cutoff |
| σ | 1.4826 × MAD, per pool; × controls' spread per window (0.5–3) | 20 | robust SD |
| Outliers | 2% | `--eps` | robust mixture |
| DONOR level | assembly d_w × a | step 19 + 20 | measured from teosinte genomes |
| a (per pool) | grid 0.5–1.8; fitted median 0.9–1.0 in every taxon | `--scale_grid` | donor plant ≠ assembly plant |
| Plants per pool | 6 | `KMER_PLANTS` | sampling design |
| SEG fractions | k/12, k ~ Binomial(12, ½) | 20 | Mendelian 1:2:1 summed over 6 plants |
| Prior π (and start) | B73 0.844, SEG 0.0625, DONOR 0.094 | `--prior` | expected state of a BC2S3 **family** = genotype of its BC2S2 parent (section 5.1) |
| ρ | grid 0.25–8; fitted **1.5** | `--rho_grid`, step 21 | maximum likelihood over all pools |
| R | none (R = 1) | — | R = 3 changed nothing |
| NOISY flag | σ > 0.15 | `KMER_NOISY_SIGMA` | left out of the ρ fit and the means |
| Skipped taxa | Zm | `KMER_SKIP_GROUPS` | no landrace assembly yet |
| Min. single-copy k-mers per window | 2,000 | `--min_b73_kmers` | drops near-empty windows |

### 4.6 Outputs (`${KMER_DIR}/introgression/`, from step 21)

- `pool_summary.tsv`: `SAMPLE ROLE GROUP SIGMA FLAG PARAM_KIND PARAM FRAC_B73 FRAC_SEG FRAC_DONOR DONOR_SHARE DONOR_SHARE_R N_TRACTS`. DONOR_SHARE = FRAC_DONOR + ½ FRAC_SEG.
- `tracts_all.tsv.gz`: `SAMPLE GROUP CHROM START END STATE`
- `rho_loglik.tsv`: log-likelihood by ρ (all pools and per taxon)
- `blocks_<GROUP>.pdf`: one page per pool. Blue B73 (AA), purple SEG (AB), red DONOR (BB).
- `family_summary.tsv`: family = first two pedigree fields, e.g. `Zd.0010_P2`. Built from `pool_summary.tsv` and `SAMPLE_META` (command in section 5.5).
- Per-pool HMM results for every ρ: `${KMER_DIR}/hmm/<SAMPLE>.npz`.

---

## 5. Results: all 1,395 pools

1,385 donor pools (Zm skipped) and 10 B73 controls. Run of 15 Nov 2026.

### 5.1 Expected for BC2S3: one plant vs one pool

Both views give the same **12.5% donor alleles**. They differ in how that 12.5% is split, because a pool is a whole family, not one plant.

**One BC2S3 plant** (top-level README). At BC2, 25% of the genome is het. Three rounds of selfing leave 1/8 of that het:

| B73/B73 | het | donor/donor | donor alleles |
|---|---|---|---|
| 85.94% | 3.13% | 10.94% | 10.94 + ½ × 3.13 = **12.5%** |

**One pool (6 plants of one BC2S3 family)**, which is what we sequence. A plot is the S3 progeny of one BC2S2 plant, so a region's state in the pool is that parent's genotype. After two selfs, 1/4 of the BC2 het regions are still het:

| Pool state | Parent (BC2S2) genotype | Expected | Plants in the pool |
|---|---|---|---|
| B73 | B73/B73 | 75% + 25% × 3/8 = **84.4%** | all AA |
| SEG | het | 25% × 1/4 = **6.25%** | mix of AA, AB, BB (1:2:1); pool is ~50% donor |
| DONOR | donor/donor | 25% × 3/8 = **9.4%** | all BB |
| donor alleles | | 9.4 + ½ × 6.25 = **12.5%** | |

So **SEG is not "het"**: a SEG region holds about 25% AA, 50% AB and 25% BB plants. The plant-level het fraction (3.13%) is the AB plants inside SEG regions: 6.25% × ½. Plant-level TT (10.94%) = DONOR (9.4%) + the BB plants inside SEG regions (6.25% × ¼). Our HMM calls pool states, so it's compared with the pool table.

### 5.2 Non-noisy pools (σ ≤ 0.15): the main results

| Taxon | Pools | Mean share % | Median share % | Mean SEG % | Mean DONOR % | Median tracts |
|---|---|---|---|---|---|---|
| B73 controls | 10 | 0.0 | 0.0 | 0.1 | 0.0 | 1 |
| Zd | 227 | 5.5 | 4.7 | 5.2 | 3.0 | 7 |
| Zh | 61 | 8.2 | 7.8 | 5.7 | 5.3 | 9 |
| Zl | 98 | 7.2 | 6.5 | 6.8 | 3.8 | 9 |
| Zv | 333 | 8.0 | 7.7 | 6.8 | 4.6 | 21 |
| Zx | 573 | 8.7 | 8.2 | 6.7 | 5.3 | 21 |
| **All donors** | **1,292** | **7.8** | **7.3** | **6.4** | **4.6** | **16** |

**Donor share across non-noisy donor pools:**

| 10% | 25% | Median | 75% | 90% |
|---|---|---|---|---|
| 2.0% | 4.2% | 7.3% | 10.5% | 14.2% |

- 199 pools reach at least 12.5%.
- 127 pools are below 2%.
- 23 pools are above 20%.

### 5.3 Noisy pools (σ > 0.15): flagged, not trusted

| Taxon | Pools | Mean share % | Median share % | Mean SEG % | Mean DONOR % | Median tracts |
|---|---|---|---|---|---|---|
| Zd | 22 | 6.5 | 5.3 | 7.6 | 2.7 | 7 |
| Zl | 20 | 7.0 | 6.3 | 7.2 | 3.3 | 10 |
| Zv | 47 | 9.5 | 9.4 | 6.7 | 6.2 | 20 |
| Zx | 4 | 17.0 | 16.3 | 19.8 | 7.1 | 64 |
| **All** | **93** | **8.6** | **8.7** | **7.6** | **4.8** | **13** |

**All 1,385 donor pools together:** mean share 7.9%, median 7.4%, SEG 6.5%, DONOR 4.6%, median 15 tracts.

**NOISY is a sequencing-plate effect.**
- **The count:** 81 of the 93 noisy pools are on three plates, **PN11, PN13 and PN15**.
- **The cause:** those plates have the fewest single-copy k-mers per pool, i.e. the least sequencing depth. Each pool is a 100-kb-window count, so less depth means a noisier score r and a larger σ.
- **The table:** all plates, from step 18 (`kmer_pool_qc.tsv`) and step 21 (`pool_summary.tsv`). Pools are those run in the HMM, including B73 controls, without Zm.

| Plate | Pools | NOISY | Median σ | Mean single-copy k-mers per pool | NOISY by taxon |
|---|---|---|---|---|---|
| PN2 | 3 | 0 | 0.118 | 189 M | — (rest of PN2 is Zm, skipped) |
| PN3 | 87 | 3 | 0.089 | 192 M | Zx 3 |
| PN4 | 93 | 0 | 0.081 | 211 M | — |
| PN5 | 92 | 0 | 0.097 | 170 M | — |
| PN6 | 94 | 2 | 0.112 | 145 M | Zx 1, Zd 1 |
| PN7 | 91 | 0 | 0.081 | 213 M | — |
| PN8 | 92 | 0 | 0.087 | 199 M | — |
| PN9 | 95 | 2 | 0.112 | 162 M | Zl 2 |
| PN10 | 94 | 0 | 0.086 | 245 M | — |
| **PN11** | 93 | **34** | **0.130** | **150 M** | Zl 18, Zv 16 |
| PN12 | 95 | 4 | 0.095 | 211 M | Zv 4 |
| **PN13** | 88 | **23** | **0.140** | **115 M** | Zd 18, Zv 5 |
| PN14 | 94 | 1 | 0.090 | 207 M | Zd 1 |
| **PN15** | 93 | **24** | **0.138** | **109 M** | Zv 22, Zd 2 |
| PN16 | 95 | 0 | 0.084 | 214 M | — |
| PN17 | 96 | 0 | 0.082 | 217 M | — |

**Recommendation: resequence PN11, PN13 and PN15.**
- **Depth:** PN13 and PN15 have about half the k-mers of a typical plate (109–115 M vs ~200 M), and PN11 about 75%. On PN11, 37% of pools are NOISY.
- **Taxa:** the noisy pools are concentrated in a few taxa:
  - Zl, 18 of its 20 noisy pools (PN11)
  - Zd, 18 of 22 (PN13)
  - Zv, 43 of 47 (PN11, PN13, PN15)
- **Families:** some are almost entirely noisy, for example Zv.0240_P2 (11 of 13 pools) and Zl.0030_P2 (7 of 10).

**Until they're resequenced,** their calls stay flagged and are left out of the ρ fit, the means and the tests.

**The list of the 93 noisy pools** (sample, plate, taxon, σ): `sample_lists/resequence_noisy_kmer.txt`. A whole-plate rerun of PN11, PN13 and PN15 would also cover their non-noisy pools, which are low depth too.

### 5.4 Are we getting what we expect?

- **SEG matches.** Observed SEG is 6.4%, against an expected 6.25%.
- **DONOR is about half the expected value.** Observed DONOR is 4.6% against 9.4%, so the donor share is 7.8% instead of 12.5%. The shortfall is almost entirely in **fixed** donor segments. Two explanations remain:
  1. **Biology:** selection against plants homozygous for teosinte segments during selfing (plant type, flowering, seed set). Some lines also lost almost all donor (127 pools below 2%).
  2. **Model:** a real DONOR block called SEG where the line's donor plant differs from the assemblies. This happens most in Zv and Zx; see section 7.
- **The low share isn't a mapping artifact.** The k-mer route never maps reads, and it gives the same shares as the SNP route for pools where the SNP route is reliable (section 6).
- **By taxon:** Zd is lowest (5.5%), and Zh, Zv and Zx highest (8–9%). Some families reach the expected 12.5%: Zx.0130_P1, Zx.0390_P2 and Zx.0560_P2.
- **Family means** range from 2.7% to 12.4% (98 families). See `family_summary.tsv`.

### 5.5 Plant genotypes and the statistical test (step 22)

**From pool states to plant genotypes.** A SEG region is a family segregating 1:2:1, so its plants are ¼ AA, ½ AB and ¼ BB. Each pool's state fractions therefore give the expected genotype mix of its plants:

- AA = B73 + ¼ SEG
- AB = ½ SEG
- BB = DONOR + ¼ SEG
- donor alleles = BB + ½ AB

These are expected plant fractions; single plants inside a pool can't be genotyped.

**The test.** Pools of one family share their donor segments, so they are not independent. The unit is therefore the **family**:
- 98 families, 1,292 non-noisy donor pools.
- The families are resampled with replacement, each with all its pools, 5,000 times (`22_kmer_stats.py`).
- 95% CI = the 2.5–97.5 percentiles.
- Two-sided p = 2 × the share of resamples on the far side of the expected value.

**Two kinds of measure:**
- **Amounts:** SEG, DONOR, AB, BB and donor alleles, as shares of the genome.
- **Ratios:** the share of the remaining donor that is still segregating, SEG/(SEG+DONOR), and its plant-level version, AB/(AB+BB) = the share of donor-carrying plants that are het. These ask how the donor that is present is split between het and homozygous, regardless of how much donor there is in total.

**All taxa (98 families, 1,292 pools):**

| Measure | Observed | 95% CI | Expected BC2S3 | p |
|---|---|---|---|---|
| SEG (pool) | 6.44% | 6.12–6.80 | 6.25% | 0.26 |
| AB (plant) | 3.22% | 3.06–3.40 | 3.13% | 0.26 |
| DONOR (pool) | 4.60% | 4.25–4.96 | 9.38% | < 0.0002 |
| BB (plant) | 6.21% | 5.81–6.61 | 10.94% | < 0.0002 |
| Donor alleles | 7.82% | 7.36–8.29 | 12.50% | < 0.0002 |
| **AB/(AB+BB)** | **34.1%** | 33.0–35.3 | 22.2% | < 0.0002 |
| SEG/(SEG+DONOR) | 58.3% | 56.7–60.0 | 40.0% | < 0.0002 |

**By taxon:**

| Taxon | Families | Pools | AB (plant) | BB (plant) | Donor alleles | AB/(AB+BB) |
|---|---|---|---|---|---|---|
| Zd | 6 | 227 | 2.58% | 4.26% | 5.55% | 37.8% (35.5–39.9) |
| Zh | 2 | 61 | 2.87% | 6.73% | 8.16% | 29.9% |
| Zl | 5 | 98 | 3.40% | 5.45% | 7.16% | 38.4% (36.4–41.3) |
| Zv | 40 | 333 | 3.40% | 6.30% | 8.00% | 35.1% (32.8–37.3) |
| Zx | 45 | 573 | 3.37% | 7.02% | 8.70% | 32.5% (31.2–33.9) |
| *Expected* | | | *3.13%* | *10.94%* | *12.50%* | *22.2%* |

BB, donor alleles and AB/(AB+BB) differ from the expectation at p < 0.0002 in every taxon. Zh (2 families), Zl (5) and Zd (6) have few families, so their intervals are less reliable. `stats_by_chrom.tsv` gives the same measures per chromosome.

**What it shows.**
- **Het is as expected.** AB (3.2%) and SEG (6.4%) don't differ from the BC2S3 expectation.
- **Homozygous donor is missing.** BB is 6.2% instead of 10.9%, and DONOR is 4.6% instead of 9.4%. That loss is the whole shortfall in donor alleles (7.8% vs 12.5%).
- **The remaining donor is mostly het.** 34% of donor-carrying plants are AB, against 22% expected. The same pattern holds in all five taxa.
- **This fits other work.** Lab mates' lines and TeoNAM reported more heterozygosity than expected after selection during selfing. Here, het is normal in absolute terms but strongly elevated relative to fixation.
- **Interpretation:** selection against plants homozygous for teosinte segments (BB) during the selfing generations. Recessive teosinte alleles that hurt plant type, flowering or seed set are exposed only in BB plants. Losing those plants removes fixed donor segments, while the same segments survive in heterozygotes.

**Why it is not a model artifact.**
1. **The prior works against it.** The HMM's prior π is the expected 84.4 / 6.25 / 9.4 split, which pulls calls toward the expectation. The departure happens in spite of it.
2. **SEG is not inflated.** If real DONOR were being called SEG, SEG would come out above 6.25%. It doesn't.
3. **The switching problem doesn't explain it.** SEG/DONOR switching occurs mainly in Zv and Zx, yet Zd and Zl, which have the least of it, show the strongest excess (38%).
4. **No false donor.** The B73 controls are clean, so B73 regions are not turned into donor.

**Next:** the per-chromosome breakdown (`stats_by_chrom.tsv`). If the BB deficit is concentrated in a few regions, those are candidate targets of selection.

### 5.6 Family table

```bash
source config.sh; PS=${D21}/pool_summary.tsv; OUT=${D21}/family_summary.tsv
tr -d '\r' < "$SAMPLE_META" | awk -F'\t' -v out="$OUT" 'NR==FNR { if (FNR > 1) { split($2, a, "_"); fam[$1] = a[1] "_" a[2] }; next }
  FNR == 1 { for (i = 1; i <= NF; i++) h[$i] = i; next }
  $h["ROLE"] == "donor" { f = ($h["SAMPLE"] in fam) ? fam[$h["SAMPLE"]] : "unknown"; n[f]++; grp[f] = $h["GROUP"]
    if ($h["FLAG"] == "NOISY") { noisy[f]++; next }
    k[f]++; sh[f] += $h["DONOR_SHARE"]; sg[f] += $h["FRAC_SEG"]; dn[f] += $h["FRAC_DONOR"]; tr[f] += $h["N_TRACTS"] }
  END { print "FAMILY\tTAXON\tPOOLS\tNOISY\tMEAN_SHARE\tMEAN_SEG\tMEAN_DONOR\tMEAN_TRACTS" > out
    for (f in n) if (k[f]) printf "%s\t%s\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.1f\n", f, grp[f], n[f], noisy[f], sh[f]/k[f], sg[f]/k[f], dn[f]/k[f], tr[f]/k[f] >> out }' - "$PS"
```

---

### 5.7 Family DONOR levels (step 20b, rounds)

**Why.** In round 0 the DONOR level comes from the teosinte assemblies. Our lines' donor plants differ from those assembly plants, and in Zv and Zx that causes SEG/DONOR switching inside tracts. Every pool of a BC1 family carries the same donor haplotype (one F1 gamete). So where pools of the family are called DONOR, their score r there is the family's real DONOR level.

**How (`20b_family_levels.py`).** For each pool and window, using the previous round's calls at its best ρ and only non-noisy pools:
- n = the number of the **other** pools of its family (leave-one-out) called DONOR in that window, and m = their median r.
- level = (n·m + K·prior) / (n + K), with K = 3.
- prior = the family's median fitted scale × the taxon assembly level.

So windows with few DONOR pools stay close to the assembly level. A pool's own calls never set its own level.

Step 20 then uses this level, with a narrow per-pool scale (0.8–1.3). Rounds repeat until the calls stop changing; usually 2–3.

**This tried a BC1-sister route first, and it did not work.** k-mer depth in BC1 plants at 10–17x didn't separate het from B73/B73 windows, even after normalizing by the pooled B73 controls (about 5x): one broad hump of depth ratios, 0.6–1.4. The expected het drop (15–35%) is about the size of the window-to-window depth noise.

**Run** (each round after the previous round's 21):
```bash
sbatch --export=ALL,ROUND=1 20b_family_levels.sh
J=$(sbatch --parsable --export=ALL,ROUND=1 --array=1-$(( (N + 29) / 30 )) --dependency=afterok:<20b jobid> 20_run_kmer_hmm.sh)
J2=$(sbatch --parsable --export=ALL,ROUND=1 --dependency=afterok:$J 21_kmer_collect.sh)
sbatch --export=ALL,ROUND=1 --dependency=afterok:$J2 22_kmer_stats.sh
```
Outputs: `KMER_DIR/family_levels_fam1/`, `hmm_fam1/`, `introgression_fam1/` (round 2: `…fam2`).

**Per chromosome (round 0, step 22):** AB/(AB+BB) is 31–37% on every chromosome (expected 22%; all p < 0.0002). The BB deficit is genome-wide, not a few regions:

| chr1 | chr2 | chr3 | chr4 | chr5 | chr6 | chr7 | chr8 | chr9 | chr10 |
|---|---|---|---|---|---|---|---|---|---|
| 34.9% | 33.4% | 31.8% | 35.7% | 35.6% | 30.9% | 34.5% | 36.7% | 34.3% | 34.1% |

## 6. Validation

- **B73 controls are clean.** All 10 have ≤ 0.15% non-B73 windows and 0–4 tiny tracts. This is a fair test: the k-mer route never uses the controls to filter sites, and each control is left out of its own reference.
  - The SNP route can't test this. Its site filter removed every site where the controls had an ALT read, so its controls are B73 by construction.
- **Known tracts are found.**
  - Zd PN13_SID1199 chr2:22–25 Mb: DONOR, r ≈ 0.20; the Zd assembly level there is 0.26.
  - Zv PN10_SID867 chr2:224–239 Mb: SEG, r ≈ 0.63. The reads confirm SEG: 55% of informative sites with 2+ reads have both B73 and ALT reads.
  - Boundaries agree with the SNP calls to within one 2-Mb bin.
- **Agreement with the SNP HMM,** in a 44-pool test (10 controls, 34 donors from all taxa):
  - The donor share is within ~1.5 points in 24 of 34 pools.
  - The disagreements have known SNP-side causes:
    - κ collapsing to 0.1–0.2 in two Zd pools (SNP 17% vs k-mer 3–7%)
    - Zm, which is too close to B73 for a teosinte panel
    - noisy pools
- **R:** no effect in any pool (section 4.4).

---

## 7. Open issues and next steps

1. **DONOR ↔ SEG switching inside tracts (Zv, Zx).** Zv and Zx have ~21 tracts per pool against 7–9 for the other taxa. Inside one long tract, some 5–10 Mb blocks are called SEG.
   - Example: PN16_SID1523, chr5 100–110 Mb (around the centromere). The pool's level is 0.53–0.76, while the mexicana assemblies predict 0.37–0.55.
   - Cause: the lines' donor plant differs locally from TIL18 and TIL25.
   - **Fix in progress:** take the expected DONOR level from the family's own **BC1 sisters** (10x, `/rsstu/users/r/rrellan/BZea/ZEAL/demultiplex/bc1`). An F1 carries one teosinte gamete, so all BC1 sisters and all BC2S3 lines of that family carry the same donor haplotype.
   - In a BC1's B73/donor windows, B73 k-mers the donor lacks appear at half depth, so d_w(family) = 2 × (depth ratio) − 1.
   - Pilot (`/share/maize/ntanduk/kmer/bc1_pilot.sh`, chr2, two BC1 plants) is running. It still needs the S_* → BC1 → family table.
   - Option if needed: estimate each family's d_w from its own DONOR pools (iteratively), with no BC1 data.
2. **Noisy plates (PN11, PN13, PN15):** resequence, or keep them flagged.
3. **Zm (landraces):** skipped. They need a maize-based DONOR level, e.g. from NAM founder assemblies.
4. **Shortfall in DONOR (4.6% vs 9.4%):** recheck after the BC1-based levels, then interpret biologically.

---

## 8. SNP route, steps 01–15 (cross-check)

These steps produced the SNP sites, per-pool read counts, QC, site filters and a GL-based HMM. Use them as a cross-check and for SNP-level work. For introgression calls, use section 5.

```bash
# site discovery (already run; only to redo it)
sbatch 01_angsd_discovery_by_part.sh; bash 02_find_failed_parts.sh; sbatch --array=1-N%60 03_angsd_rerun_failed_split.sh
sbatch 04_merge_subparts.sh; sbatch 05_concat_parts_by_chr.sh; sbatch 06_rename_samples.sh
# read counts
bash submit_counts.sh test; bash submit_counts.sh rest; bash 08_check_unfinished.sh
# site filters + teosinte panel
sbatch 11_filter_sites.sh; sbatch --dependency=afterok:<11> 12_signal_check.sh
# SNP HMM
sbatch 13_extract_pool_gl.sh; sbatch --dependency=afterok:<13> 14_run_hmm.sh; sbatch --dependency=afterok:<14> 15_collect_hmm.sh
```

| # | Script | What it does | Output |
|---|---|---|---|
| 01–06 | `01_…` to `06_…` | ANGSD SNP discovery + GLs per part (`-GL 1 -doMajorMinor 1 -doMaf 1 -doBcf 1`, `-minMapQ 30 -minQ 20 -C 50 -baq 1 -only_proper_pairs 1`); reruns of failed parts; joined and renamed per chromosome | `06_per_chr_renamed/` |
| 07 | `07_make_site_tables.sh` + `07_polarize_to_b73.py` | Biallelic SNPs re-polarized to the B73 base; ANGSD `-sites` files | `07_sites/` |
| 08–09 | `08_angsd_counts_by_part.sh` + `08_counts_to_ad.awk`, `09_…` | B73/ALT read counts per pool (`-doCounts 1 -dumpCounts 4 -sites`) | `09_counts_by_chr/BZea_chrN.AD.bcf` |
| 10 | `10_counts_qc.sh` + `.R` | Per-pool QC | `10_qc/` |
| 11 | `11_filter_sites.sh` | Removes sites with ALT in clean B73 controls, depth outside 0.25–2× the median, >5% third base, or a panel allele mismatch; splits sites into in / not in the 227-plant teosinte panel (Schnable 2023) | `11_sites_filtered/` (7.51 M PASS, 2.59 M in panel; `pools.tsv`) |
| 12 | `12_signal_check.sh` + `.R` | Donor signal per taxon on panel sites | `12_signal_check/` |
| 13–15 | `13_extract_pool_gl.sh`, `14_run_hmm.sh` + `14_hmm_pooled.py`, `15_collect_hmm.*` | GL-based HMM | `15_introgression/` |

**How the SNP HMM (step 14) works:**
- **States:** B73 / SEG / DONOR. The pool likelihood comes from interpolating the GLs: L(p) = L0(2p−1)(p−1) + L1·4p(1−p) + L2·p(2p−1).
- **Donor ALT chance:** c = κ × panel AF, with κ learned per pool by multi-start EM.
- **Transitions:** NAM map transitions, ρ = 0.2, rigidity 50 markers.

**Results and limits:**
- **Results:** median donor share Zd 5.4%, Zh 8.0%, Zl 6.7%, Zv 8.9%, Zx 8.6%.
- **Limit 1:** it loses donor reads to mapping (section 2).
- **Limit 2:** κ can collapse.
- **Limit 3:** its control test is circular.

**Dropped:** steps 13b–15b, which used the BC1 family donor alleles (`union/*.donor_alleles.vcf.gz`) as SNP-level PALT. The family alleles were wrong at 7–11% of sites even in good pools, and the run lowered DONOR in every taxon. The scripts are kept but not used.

**QC decisions:** [`QC_NOTES.md`](QC_NOTES.md); lists in `sample_lists/`.

---

## 9. Requirements

- **k-mer route:**
  - KMC 3 (`kmc`, `kmc_tools`), KCFtools and samtools, in the conda env `/rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk/kmer_env`
  - python3 with numpy, scipy and matplotlib (steps 20–21)
- **SNP route:** angsd, bcftools/htslib, awk, base R, python3 (standard library)
- **Slurm** for all jobs. No partition is set; uncomment `##SBATCH --partition=...` if needed.
