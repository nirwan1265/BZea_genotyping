# BZea: k-mer introgression calling and QTL mapping in maize × teosinte BC2S3 pools

This repository goes from raw paired-end reads of **pooled, low-coverage (~0.8×)** BC2S3 samples to **teosinte introgression calls** and **QTL** for flowering, height, ear and leaf traits. B73 was crossed to 5 teosinte taxa, backcrossed twice to B73 and selfed three times; each sample is a pool of 6 plants from one BC2S3 plot.

- **Introgressions are called from k-mers**, not from mapped SNPs: for each pool, which B73 single-copy 31-mers its reads contain, per 100-kb window, decoded with an HMM that models the pool: each window is **B73** (all 6 plants AA), **DONOR** (all BB) or **SEG** (segregating: the plants are a mix of AA, AB and BB).
- **QTL are mapped on that k-mer ancestry**: a genome scan (k-mer ancestry GWAS, GridLMM with leave-one-chromosome-out kinship) finds the peaks, and a joint QTL model (R/qtl2) fits them together, with an effect for each teosinte taxon. Each QTL interval is then fine-mapped (SuSiE).

> **This README describes the current model (adopted 20 November 2026).** It differs from the first model in two places, both in how the HMM is fed, not in the HMM itself: where the teosinte shares most of the B73 k-mers, only the **informative** k-mers are counted (B73 k-mers missing from every assembly of the taxon), and the DONOR level of each window is **learned from our own pools** instead of taken from assemblies of other plants (section [2.6.6](#266-the-current-model-informative-k-mers-and-donor-levels-from-our-pools)). Checked against targeted sequencing of the same plots (three bad samples left out), SEG calls that are heterozygous rose from 74% to 86% and DONOR calls that are homozygous teosinte from 91% to 95% (section 2.9.1); the genome-wide shares and the GWAS scans hardly change. The first model, with all its results, is kept in [`old_model_README.md`](old_model_README.md) (results in `data/introgression_spread/`).

> **Pipeline and commands from the k-mer and SNP steps on:** [`final_script/`](final_script/README.md) (the pre-processing scripts are not in the repository yet; see Section 1). Results tables: [`data/introgression_fix1/`](data/introgression_fix1/). Figures: [`figures/`](figures/).

## Contents

- [Words used in this README](#words-used-in-this-readme)
- [The BZea population](#bzea-breeding-scheme): breeding scheme and expected genotype proportions
- [Overview](#overview) and [software](#software-and-packages)
- [Section 1 — Pre-processing](#section-1--pre-processing): demultiplexing, trimming, alignment
- [Section 2 — Introgression calling](#section-2--introgression-calling): the k-mer HMM, results, validation, population structure
- [Section 3 — k-mer ancestry GWAS and joint QTL mapping](#section-3--k-mer-ancestry-gwas-and-joint-qtl-mapping): phenotypes, methods, results for 17 traits (11 in detail), *mexicana* races, candidate genes
- [Section 4 — References](#section-4--references)

## Words used in this README

- **Pool** or **sample:** one sequenced sample, i.e. one leaf each from 6 plants of one plot. Every result is per pool.
- **Window:** a 100-kb piece of the B73 genome. The HMM calls one state per pool per window.
- **Tract:** a run of neighbouring windows with the same state.
- **F1 family:** all pools that descend from one F1 plant (same accession and same first P in the name, see below). Used to group pools in the statistics.

## BZea Breeding Scheme

B73 was crossed with teosinte, backcrossed twice to B73 and then selfed three times:

**B73 × teosinte → F1 → BC1 (× B73) → BC2 (× B73) → BC2S1 → BC2S2 → BC2S3**

- **Selfing, ear to row:** a plant is selfed, its progeny are planted, several are selfed, **one is selected** and advanced to the next generation, and so on. So each BC2S3 plot is grown from the seed of one BC2S2 mother (one selfed ear).
- **Sequenced sample:** one leaf each was taken from 6 different BC2S3 plants of the plot and pooled into one sample. The 6 plants were not random: the best-looking plants of the plot were chosen. These 6 plants are siblings, not copies of one plant: each got its own two gametes from the BC2S2 mother. So where the mother is heterozygous, each plant is AA, AB or BB with chance ¼, ½, ¼.
- **Pedigree names:** in `Zv.0200_P3_P2_P1.1.1.1`, `Zv.0200` is the teosinte accession, `P3` the F1 plant, `P2` the BC1 plant, `P1` the BC2 plant, and `.1.1.1` the selfed line. Pools that share `Zv.0200_P3` descend from one F1 plant (an **F1 family**: they all carry the same teosinte haplotype); pools that share `Zv.0200_P3_P2_P1` come from the same BC2 plant.

## Expected genotype proportions (per locus)

**Notation, used everywhere in this README and in `final_script/`:**
- **A** = the B73 allele, **B** = the teosinte (donor) allele
- **AA** = homozygous B73, **AB** = heterozygous, **BB** = homozygous teosinte

**Assumptions:** B73 is AA, the teosinte donor contributes one B allele to the F1, and there is no selection or segregation distortion. (In practice one plant was selected at each selfing step, so a departure from these numbers can point to selection, but also to segregation distortion, lower survival or fertility of some genotypes, or errors in the calls; see [Results](#28-results).) These are the same for every donor taxon: Zd (*diploperennis*), Zh (*huehuetenangensis*), Zl (*luxurians*), Zv (*parviglumis*) and Zx (*mexicana*).

### Step 1 — Two backcrosses (BC2)

- F1 = B73 (AA) × teosinte (BB) → **100% AB**
- BC1 = F1 × AA → ½ AA, ½ AB
- BC2 = BC1 × AA → **¾ AA, ¼ AB**, 0 BB

Donor allele share at BC2: ¼ × ½ = **12.5%**. Selfing doesn't change allele shares, so 12.5% is the expectation at every later generation too.

### Step 2 — Selfing a heterozygote

Each self halves the heterozygotes; the rest become AA or BB in equal parts. Starting from AB, after *s* selfs:
- AB = (½)^s
- AA = BB = (1 − (½)^s) / 2

| Selfs | AA | AB | BB |
|---|---|---|---|
| 1 (BC2S1) | 1/4 | 1/2 | 1/4 |
| 2 (BC2S2) | 3/8 | 1/4 | 3/8 |
| 3 (BC2S3) | 7/16 | 1/8 | 7/16 |

### Two levels: plants and pools

Everything below is computed at two levels, and it helps to keep them apart:
- **Plant level:** the genotype of one BC2S3 plant at a locus, **AA**, **AB** or **BB**.
- **Pool level:** what the sequenced sample shows. A sample is a pool of 6 sibling BC2S3 plants, and the HMM gives each 100-kb window of a pool one of three **pool states**:
  - **B73:** all 6 plants are AA.
  - **DONOR:** all 6 plants are BB.
  - **SEG** (*segregating*): the plants differ. Their mother was AB, so each of the 6 plants is AA, AB or BB. A SEG pool is a **mix of plants**, not one heterozygous plant.

### Step 3 — One BC2S3 plant (three selfs)

Start from the BC2 numbers in Step 1: **¾ of the genome is AA** and **¼ is AB** (BC1 = ½ AA + ½ AB; crossing to B73 again turns the AA half into AA and splits the AB half into AA and AB, so BC2 = ½ + ¼ = ¾ AA and ¼ AB).

Then the BC2 plant is selfed three times (BC2S1, BC2S2, BC2S3):
- the **¾ that is AA** stays AA: selfing AA only gives AA.
- the **¼ that is AB** segregates, as in the **3 selfs** row of Step 2: of that ¼, 7/16 becomes AA, 1/8 stays AB and 7/16 becomes BB.

| Genotype | From the ¾ AA at BC2 | From the ¼ AB at BC2 | Expected |
|---|---|---|---|
| AA | ¾ | ¼ × 7/16 | 0.75 + 0.1094 = **85.94%** |
| AB | 0 | ¼ × 1/8 | **3.13%** |
| BB | 0 | ¼ × 7/16 | **10.94%** |
| Donor alleles | | BB + ½ × AB | 10.94 + 1.56 = **12.5%** |

### Step 4 — One pool (6 BC2S3 siblings)

The 6 plants in a pool are the selfed offspring of one **BC2S2 mother**, so the pool state at a locus is set by **her** genotype:
- mother **AA** → all her offspring are AA → pool state **B73**
- mother **BB** → all her offspring are BB → pool state **DONOR**
- mother **AB** → her offspring segregate 1 AA : 2 AB : 1 BB (Mendel, selfing a heterozygote), so each of the 6 plants is AA with chance ¼, AB with ½ and BB with ¼ → pool state **SEG**

The mother is a BC2S2 plant, i.e. **two** selfs after BC2. Same logic as Step 3, with the **2 selfs** row of Step 2: the ¾ of the genome that is AA at BC2 stays AA, and the ¼ that is AB at BC2 becomes 3/8 AA, 1/4 AB and 3/8 BB:

| BC2S2 mother | From the ¾ AA at BC2 | From the ¼ AB at BC2 | Expected share of the genome | The 6 plants | Pool state |
|---|---|---|---|---|---|
| AA | ¾ | ¼ × 3/8 | 0.75 + 0.0938 = **84.38%** | all AA | **B73** |
| AB | 0 | ¼ × 1/4 | **6.25%** | AA, AB or BB with chance ¼, ½, ¼ (1:2:1) | **SEG** |
| BB | 0 | ¼ × 3/8 | **9.38%** | all BB | **DONOR** |

### Step 5 — From pool states back to plants

Steps 3 and 4 describe the same plants, so one can be computed from the other. B73 pools contain only AA plants and DONOR pools only BB plants. SEG pools contain the 1:2:1 mix, so on average **¼ of a SEG region is AA in a plant, ½ is AB and ¼ is BB**. Adding up, with SEG = 6.25%:

| Plant genotype | From the pool states | Expected |
|---|---|---|
| AA | B73 + ¼ × SEG = 84.38 + 1.56 | **85.94%** |
| AB | ½ × SEG = ½ × 6.25 | **3.13%** |
| BB | DONOR + ¼ × SEG = 9.38 + 1.56 | **10.94%** |
| Donor alleles | DONOR + ½ × SEG = 9.38 + 3.13 | **12.5%** |

These are exactly the Step 3 numbers. This is how the HMM results are reported at both levels: the pool states as called (compared with Step 4), and converted to plant genotypes with this table (compared with Step 3). Heterozygosity (AB) can only be estimated this way, as half of SEG; no single plant is genotyped. So AB, BB and AA in this README are **Mendelian estimates from the HMM pool states**, not measured plant genotypes, and AB always equals half of SEG by construction. (The HMM's SEG state leaves out pools in which all 12 gametes happened to carry the same allele, k = 0 or 12; under Mendel that is 2 × ½¹² = 0.05% of SEG windows, small enough to ignore.)

### Implication for low-pass pooled genotyping

- **Het is small:** only 3.1% of a plant's genome is expected to be AB (Step 3). At 0.8x, single-site het calls are unreliable anyway, so we call pool states over 100-kb windows instead of genotypes at sites.
- **Most donor should be fixed:** of the 12.5% donor alleles, most should sit in DONOR regions (9.38%), where all 6 plants are BB. SEG regions (6.25%) hold the rest, as a 1:2:1 mix of plants (Step 5).
- **What we found:** het matches the expectation (SEG 6.23%, AB 3.12%), but DONOR is only about half of it (5.0% vs 9.4%). See [Results](#28-results).

---

## Overview

### 1. Pre-processing (Section 1)
- **Demultiplexing** by barcode (`sabre`) into per-sample paired FASTQs.
- **Trimming** of adapters, primers and low-quality bases (Trimmomatic).
- **Alignment** to B73 RefGen v5 (BWA-MEM), sorting, duplicate removal, indexing; used by the SNP route only.

### 2. Introgression calling (Section 2)
- **k-mer windows** (KMC 3 + KCFtools): for each pool, the B73 single-copy 31-mers seen in its reads, per 100-kb window; no read mapping (`final_script/` steps 16–19).
- **Informative k-mers** per taxon (B73 single-copy k-mers missing from every assembly of the taxon), counted where they separate DONOR from SEG better than all B73 k-mers (`new_model/fix1/`).
- **Expected DONOR level** per window: from the teosinte genome assemblies (MaizeGDB / PanAnd) to start with, then learned from all pools of each taxon (`new_model/20c_pool_levels.py`, `new_model/fix1/`).
- **HMM** with three states, B73 / SEG / DONOR, genetic-map transitions and pool-aware emissions (steps 20–21). Outputs: per-pool state per window, tracts, `pool_summary.tsv`, `family_summary.tsv`, chromosome paintings, interactive figures.
- **Population structure from k-mers** (sourmash; step 32).
- **SNP route** (ANGSD sites and read counts, GL-HMM; steps 01–15), kept as a cross-check.

### 3. k-mer ancestry GWAS and joint QTL mapping (Section 3)
- **Phenotypes:** SpATS spatially corrected BLUEs per year (2023, 2025), combined across years (step 28).
- **k-mer ancestry GWAS:** GridLMM, taxon as fixed effect, **LOCO kinship** from the ancestry dosage, 5% genome-wide threshold from 1,000 permutations (steps 25–27); rerun with *mexicana* split by race.
- **Joint QTL model:** the peaks fitted together in R/qtl2: effect of each QTL adjusted for the others, its effect in each taxon, and a test that the taxa differ (step 29).
- **Genes and GO terms** in each QTL interval, side-by-side figures (steps 30–31).
- **Fine-mapping** of each QTL interval with SuSiE: posterior inclusion probability (PIP) per window, credible sets, genes per QTL (steps 33, 33a, 33b).

## Software and packages

### Genotyping

- **sabre** — demultiplexing by barcode (FASTQ splitting)
- **Trimmomatic** — adapter/primer trimming + quality filtering
- **bwa-mem2** *(recommended)* or **bwa** — read alignment to reference genome
- **samtools** — BAM/SAM manipulation (sort, index, stats, depth)
- **picard** — duplicate marking
- **bcftools** — querying
- **ANGSD** — genotype likelihood estimation + SNP calling for low-pass data
- **vcftools**  — quick VCF filters / summaries
- **bgzip / tabix** *(htslib)* — compress + index VCF/BCF outputs

### Imputation and Variant handling

- **plink2** — variant QC, LD pruning, PCA, relationship matrices, LOCO kinship (`.rel/.rel.id`)
- **Beagle** (Java) — genotype imputation

### Introgression inference (k-mers + HMM)

- **KMC 3** (`kmc`, `kmc_tools`): k-mer counting
- **KCFtools**: k-mer counts per reference window
- **Python** (>= 3.8) with numpy, scipy, matplotlib: HMM and figures
- **R** (>= 4.1): QC and SNP-route summaries

### Phenotypes and QTL / association mapping

- **R** (>= 4.1) with **SpATS** (spatial BLUEs), dplyr, tidyr
- **GridLMM** R package (`github.com/deruncie/GridLMM`), data.table, ggplot2
- **R/qtl2** (joint QTL model) and **susieR** (fine-mapping)
- **Python** with pandas, numpy: inputs (dosage, LOCO kinship) and regions/genes

### Gene annotation resources

- **GFF3 annotation** for the reference build (we used v5 annotation - `Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3`)

---

## Section 1 — Pre-processing

The sequencing lanes hold many pooled samples. Pre-processing splits them into one pair of FASTQ files per sample, trims them, and (for the SNP route only) aligns them to B73. The two downstream routes use different products:

| Product | Used by |
|---|---|
| Trimmed paired reads, `filtered_S<plate>/<SAMPLE>/<SAMPLE>_paired_R{1,2}.fq.gz` | **k-mer route** (Section 2): k-mer counting (step 17) and k-mer sketches (step 32a). No alignment needed. |
| Sorted, deduplicated, indexed BAMs (`bam_list.txt`) | **SNP route** (cross-check): ANGSD site discovery and read counts (steps 01–12). |

Samples are named `PN<plate>_SID<number>` (e.g. `PN12_SID1106`); `name_change_final.txt` maps them to the pedigree names (Section 3).

### 1.1 Demultiplexing

Each lane's paired FASTQs (`..._R1_001.fastq.gz`, `..._R2_001.fastq.gz`) are split by the inline barcode at the start of read 1 with **sabre** (paired-end mode). The barcode file has two tab-separated columns, the barcode sequence and the sample name:

```txt
ACGTACGT    PN12_SID1106
TGCATGCA    PN12_SID1107
```

sabre writes one read pair per sample and puts reads whose barcode matches no sample into `no_bc_match_*_R1.fq.gz` / `_R2.fq.gz`; their share is a quick check of barcode quality per lane.

### 1.2 Trimming and quality filtering

**Trimmomatic** (paired-end mode) removes adapter and primer sequence and low-quality bases:

| Setting | What it does |
|---|---|
| `ILLUMINACLIP:<adapters.fa>:2:30:10` | removes adapter/primer matches (2 seed mismatches, palindrome score 30, simple score 10) |
| `LEADING:3`, `TRAILING:3` | trims bases below quality 3 at the read ends |
| `SLIDINGWINDOW:4:15` | cuts the read when the mean quality of 4 bases falls below 15 |
| `MINLEN:36` | drops reads shorter than 36 bp |

Only pairs in which both reads survive (`*_paired_R1.fq.gz`, `*_paired_R2.fq.gz`) are kept; these are the reads used by the k-mer route.

### 1.3 Alignment and BAM processing (SNP route)

The trimmed pairs are aligned to **B73 RefGen v5** with **BWA-MEM**. Each BAM is then coordinate-sorted, PCR duplicates are removed, and the BAM is indexed (SAMtools and Picard). The k-mer route does not use the BAMs; they feed the SNP route only, which is kept as a cross-check (Section 2, "The problem").

> **Scripts:** the pre-processing commands (sabre, Trimmomatic, BWA, BAM processing) are not yet in this repository; they will be added under `scripts/`. The outputs are at `/rsstu/users/r/rrellan/sara/BZea/filtered_S` (trimmed reads) and `/rsstu/users/r/rrellan/DOE_CAREER/BZea/dedup_filtered` (BAMs) on the NCSU cluster.

---

## Section 2 — Introgression calling

All code is in [`final_script/`](final_script/), numbered in run order. **The detailed analysis (every command, parameter, output, all statistics and the method history) is in [`final_script/README.md`](final_script/README.md).** This section explains what we do and why, and the main results.

### 2.1 What a sample is

Each sequenced sample is a **pool of 6 BC2S3 plants from one plot**, sequenced at about **0.8x**, with no parental sequencing. B73 is the recurrent parent. 1,385 donor pools were called, from 5 teosinte taxa (pools per race from `data/metadata_BZea.csv`):
- *mexicana* (Zx), 577 pools from four races:
  - **Chalco:** 14 accessions, 119 pools (the highland race, median 2,390 m)
  - **Mesa Central:** 16 accessions, 167 pools
  - **Durango:** 3 accessions, 146 pools
  - **Nobogame:** 2 accessions, 142 pools
  - 3 pools are from accessions with no race in the metadata.
- *parviglumis* (Zv): 380 pools from 36 accessions, all race Balsas.
- *diploperennis* (Zd): 249 pools from 4 accessions.
- *luxurians* (Zl): 118 pools from 5 accessions (Guatemala 4, Oaxaca 1).
- *huehuetenangensis* (Zh): 61 pools from 2 accessions.
- Landraces (Zm) are set aside for now.

The *mexicana* races are kept apart in Section 3 as an option (`SUBTAXON`), because they differ in origin and elevation.

Because a sample is several sibling plants, each genome region in a pool is in one of three states, set by the genotype of the plot's BC2S2 mother (B73, SEG or DONOR; see [Step 4](#step-4--one-pool-6-bc2s3-siblings)). The expected shares, at the pool level and for one plant:

| | Only B73 alleles | Both alleles | Only donor alleles | Donor alleles |
|---|---|---|---|---|
| **Pool state** (6 plants) | **B73** (all plants AA): 84.4% | **SEG** (plants AA, AB or BB, 1:2:1): 6.25% | **DONOR** (all plants BB): 9.4% | 12.5% |
| **One plant** | **AA:** 85.9% | **AB:** 3.1% | **BB:** 10.9% | 12.5% |

SEG is not the same as AB: a SEG pool holds a mix of AA, AB and BB plants, and only half of a SEG region is AB in an average plant (AB = ½ × 6.25% = 3.1%).

**What SEG means.** SEG is a state of the *pool*, not of a plant, and it is **not the same as "het"**.
- In a SEG region, the family's BC2S2 mother was AB, so its BC2S3 offspring segregate 1:2:1. **Each of the 6 plants in the pool can be AA, AB or BB.**
- On average a SEG pool holds ¼ AA, ½ AB and ¼ BB plants. The exact mix differs from pool to pool (e.g. 2 AA, 3 AB, 1 BB), and the model allows for that: the pool carries k donor copies out of 12 (6 plants × 2), with k ~ Binomial(12, ½).
- A SEG call says "this family is still segregating here". It does not say that any one plant is AB.
- **SEG vs "het" in SNP-based callers.** In a SEG region the pool's reads show both alleles, just like a het site in a single plant, so a diploid SNP genotyper would call it het. But in a pool, "both alleles present" only means the 6 plants are a mix: any combination of AA, AB and BB plants, as long as both alleles are there. That's why we call this state SEG, not het.
- To get plant genotypes from pool states: **AA = B73 + SEG/4, AB = SEG/2, BB = DONOR + SEG/4**.

**Can we tell which plants are AA, AB or BB in a SEG region?** No, not from this data. The 6 leaves are ground into one sample, so a read doesn't say which plant it came from.
- **What the pool can tell:** roughly how many donor copies it carries, k out of 12. In principle k can be estimated per SEG tract, because a 10-Mb tract has about 100 windows. Right now the HMM averages over all possible k and doesn't report one.
- **Even an exact k doesn't give the genotypes.** For example, 6/12 could be 6 AB plants, or 3 AA + 3 BB, or 2 AA + 2 AB + 2 BB.
- **The estimate is only as good as the family's DONOR level d**, because r(SEG) = d + (1 − d)·g(k/12) (section 2.5).
- **k can change inside a tract.** The 12 gametes recombine independently, so a plant can be AA at one end of a SEG tract and AB at the other.
- **Leaves are not equal.** Leaf size and DNA yield differ, so a plant's share of the reads isn't exactly 1/6. This adds noise to k.
- **To get plant genotypes:** genotype single plants, or seed of that family, e.g. with a few KASP markers inside the SEG tracts, or skim-sequence single plants from a subset of families.

### 2.2 The problem: mapping reads and calling SNPs at 0.8x in pools

We first did it the usual way (steps 01–15): map reads to B73, count each pool's reads at SNPs (ANGSD), then run an HMM. It found too little donor and too much het. Three things go wrong, and they hit donor (BB) and het regions hardest:

1. **Donor reads are lost (reference bias).** In a known Zd DONOR tract (PN13_SID1199, chr2:22.1–25.3 Mb), only **~34%** of the expected reads map there at all, because much teosinte sequence is missing from or rearranged against B73. The usual filters (`-C 50`, BAQ, MAPQ ≥ 30) then remove most of the rest: only **~11%** of expected donor bases survive, against ~96% in B73 regions. The filters remove reads with mismatches, and a read carrying a teosinte allele has extra mismatches by definition. So ALT alleles are lost first: **DONOR (BB) shows mixed reads and is called SEG, and SEG loses its ALT reads and is called B73.** This comes from mapping to one reference; tuning filters can't fix it.

   Read counts in that tract, compared with a B73 region of the same size and normalized with a B73 control sample:

   | What is kept in the donor tract | Share of the expected reads or bases |
   |---|---|
   | Reads that map there at all | 34% |
   | + MAPQ ≥ 30 | 30% |
   | Bases, MAPQ ≥ 30, no `-C 50`, no BAQ | 25% |
   | Bases with `-C 50` + BAQ (the ANGSD settings used) | **11%** |
   | Same filters in a B73 region | ~96% |

   The k-mer method (below) calls this tract DONOR, with r ≈ 0.20 (the Zd assemblies give 0.26 there).
2. **Too few reads per site.** At 0.8x most SNPs have 0 or 1 read. One read shows one allele, so a het site looks homozygous, and a single error read looks like a het. Every call leans on the HMM and imputation.
3. **Pools aren't diploid.** A pool's ALT fraction can be anything from 0 to 1 (e.g. 5/12 in a SEG region). Genotype callers built for diploids (0 / ½ / 1) mislabel these mixtures.

### 2.3 The solution: k-mers

Finding introgressions by counting reference k-mers in genome windows, without mapping reads, was developed in wheat: **IBSpy** (Uauy lab; Quiroz-Chávez 2023) was used to find einkorn introgressions in bread wheat (Ahmed et al. 2023), and **KCFtools** (Selvanayagam et al. 2026) is a fast tool for k-mer introgression screening and k-mer GWAS, which we use for the per-window counts. k-mers are counted with **KMC 3** (Kokot et al. 2017). The DONOR levels come from teosinte genome assemblies (PanAnd-1.0, Stitzer et al. 2025), downloaded from **MaizeGDB** (Portwood et al. 2019).

**What we added:** low-coverage *pools* instead of single plants, a three-state pool model (B73 / SEG / DONOR) with per-window DONOR levels (from assemblies, then learned from our own pools), counting only the informative k-mers where the teosinte shares most of B73's (section 2.6.6), and the HMM on the NAM genetic map.

### 2.4 How it works

![How k-mers find introgressions](final_script/figures/kmer_method.svg)

1. **B73 k-mer list (step 16).** The **B73 reference genome** (RefGen v5, `Zm-B73-REFERENCE-NAM-5.0`, the same reference used for mapping in the SNP route) is cut into every 31-letter word (k-mer), counted with KMC. Only the **single-copy** ones are kept: 680 M of 817 M. Repeats are dropped because they're in every genome, donor or not, and would hide the signal. Spread over the genome's 21,330 windows of 100 kb, that is **about 32,000 single-copy k-mers per window** on average (it varies from window to window).
2. **Look them up in each sample (step 17).** All k-mers in a sample's raw reads (the pool of 6 plants) are counted, with no mapping. For each 100-kb B73 window, KCFtools counts how many of the window's single-copy B73 k-mers the sample has.
   - **B73 region:** essentially all of them are present, limited only by sequencing depth.
   - **Donor region:** every teosinte SNP or indel removes the B73 k-mers that overlap it, so the count drops.
3. **Score each window (step 20).** At 0.8x even a pure-B73 sample finds only part of the window's k-mers, so the pool's count is not compared with the reference's 32,000. It is compared with what **our own sequenced B73 samples** find:
   - **B73 control samples:** 10 samples of B73 itself, sequenced by us the same way as the BC2S3 pools. In each window they show how many of the reference k-mers a pure-B73 sample actually finds, rescaled to the depth of the pool being scored. A control is never used to score itself.
   - **r = k-mers found in the pool ÷ k-mers the B73 controls find** (at the pool's depth).
   - So two kinds of B73 are used: the **B73 reference genome** decides *which* k-mers are counted (step 1), and the **B73 control samples** say *how many* of them a B73 sample at this depth should find.
   - B73 gives **r ≈ 1**, DONOR **r ≈ d**, and SEG falls between them.
   - d is the share of B73 k-mers a donor genome keeps. It depends on the taxon and on the window, and is measured per window from teosinte genome assemblies (PanAnd-1.0 from MaizeGDB, step 19). Genome-wide: Zd 0.27, Zl 0.26, Zh 0.34, Zx 0.37, Zv 0.43.
   - SEG: k-mers the donor lacks come only from the B73 plants in the pool. For 6 plants, the pool's donor fraction follows Binomial(12, ½).
4. **HMM along each chromosome (steps 20–21).** The windows' r values are read along each chromosome by a three-state hidden Markov model (B73 / SEG / DONOR). It uses the NAM genetic map, so a state change costs more where recombination is rare. Adjacent windows in the same state become tracts. **The model is explained step by step in [The HMM in detail](#26-the-hmm-in-detail) below.**
5. **Allow for our donor plants differing from the assembly plants (step 20).** Where the two assemblies of a taxon disagree on d, the model is told that the DONOR level there is uncertain, so it doesn't flip between SEG and DONOR on small differences (see τ_w under [Emissions](#261-emissions-how-well-does-r-fit-each-state)). In the current model, d is then learned from all our pools of the taxon, and where the teosinte keeps most of B73's k-mers only the informative ones are counted ([2.6.6](#266-the-current-model-informative-k-mers-and-donor-levels-from-our-pools)).
6. **Tests (step 22).** Observed vs expected BC2S3 proportions, with a bootstrap over F1 families. Pools of one F1 family share segments, so the F1 family is the unit.

### 2.5 How r is calculated: a worked example

An average window in a typical Zd pool. The numbers are genome-wide averages from our data, rounded; real windows vary around them.

**Step 1: which k-mers are counted (B73 reference).** The B73 reference has 680 M single-copy k-mers in 21,330 windows, so an average window has **~32,000** (step 16). This is the same for every pool and every taxon.

**Step 2: how many a B73 sample at this depth finds (B73 controls).** At 0.8x many parts of the window are not covered by any read, so a sample finds only part of the 32,000. A typical pool contains about 200 M of the 680 M single-copy k-mers, i.e. **about 30%** (plate averages range from 109 M to 245 M; section 2.11). So in this window our B73 control samples, rescaled to the pool's depth, find about 30% × 32,000 ≈ **9,600**. That's the expected count. It depends on the pool's depth, not on its taxon.

**Step 3: r = what our pool found ÷ the expected count.**

| Our pool in this window is… | k-mers found | r |
|---|---|---|
| B73 (all plants AA) | ~9,600 | 9,600 / 9,600 ≈ **1.0** |
| DONOR (all plants BB, Zd) | ~2,600 | 2,600 / 9,600 ≈ **0.27** |
| SEG (about half B73, half donor) | ~6,300 | 6,300 / 9,600 ≈ **0.66** (below) |

**Where the DONOR level comes from (step 19).** We take the complete Zd genome assemblies (Gigi and Momo, PanAnd-1.0 from MaizeGDB) and count how many of the window's ~32,000 B73 k-mers appear anywhere in them. Genome-wide it's 27% (Gigi 0.26, Momo 0.28), i.e. about 8,600 of 32,000: **d = 0.27**. A pool that is all Zd in this window can only find that 27%, so it finds 0.27 × 9,600 ≈ 2,600 and its expected r is 0.27. A real example: the Zd tract of PN13_SID1199 (chr2:22–25 Mb, section 2.2) has r ≈ 0.20, where the Zd assemblies give 0.26.

d depends on the taxon and on the window:
- **Taxon:** genome-wide Zd 0.27, Zl 0.26, Zh 0.34, Zx 0.37 and Zv 0.43. Zv is closest to maize, so it keeps the most.
- **Window:** for Zd, the per-window median is 0.23, and 80% of windows are between 0.13 and 0.42. The HMM uses each window's own d.

Sequencing the **actual donor parents** of our lines would have been better, because d would then be measured on the plants our lines came from. We don't have them, so we use the MaizeGDB assemblies of the same taxa.

**Why SEG falls in between.** In a SEG region the pool is a mix of B73 and donor chromosomes. Split the window's k-mers in two:
- **Shared with the donor (27%):** every chromosome carries them, so they're found as often as in a B73 sample (ratio 1.0).
- **Missing from the donor (73%):** only the B73 chromosomes carry them, about half the chromosomes, so they get about half the reads. Half the reads doesn't mean half the k-mers found, because at low depth many k-mers are missed either way:
  - A B73 sample at this depth finds **30%** of the k-mers. From the coverage formula (the chance a k-mer is seen at coverage c is 1 − e^(−c)), that means a k-mer coverage of c = 0.36.
  - With half the reads (c = 0.18), it finds 1 − e^(−0.18) = **16%**.
  - So these k-mers are found 16.3 / 30 = **0.54** times as often as in a B73 sample. That ratio is **g**.

So r(SEG) = 0.27 × 1.0 + 0.73 × 0.54 ≈ **0.66**.

In general, r(SEG) = d + (1 − d) × g. g depends on the pool's depth, so each pool gets its own. The HMM doesn't assume exactly half donor. It averages over the possible donor fractions of a 6-plant pool (1/12 to 11/12), weighted by Binomial(12, ½), the 1:2:1 segregation summed over 6 plants.

**Then the HMM** compares each window's r with these three expected levels (B73 ≈ 1, SEG ≈ 0.66, DONOR ≈ d). It allows for each pool's noise and takes neighbouring windows into account, using the genetic map, to call the most likely state.

### 2.6 The HMM in detail

A hidden Markov model (HMM) is used because we can't see a pool's true state directly. We only see a noisy score r in each window. But we know two things that help:
- each state gives r values in a typical range, and
- the state rarely changes from one window to the next, because introgressed segments are long compared with a 100-kb window (median about 10 Mb in our pools; see [Tract lengths](#265-how-long-are-the-tracts) below).

The HMM combines both to find the most likely series of states along each chromosome. It has the same three-state structure as the SNP HMM in step 14 (`14_hmm_pooled.py`); see Rabiner (1989) for HMMs in general. Only the observation is different: a k-mer score per window instead of SNP reads. The design follows the HMM that Zhong et al. (2025) built to call introgressions in 888 maize near-isogenic lines with B73 as the recurrent parent: three states (homozygous B73, heterozygous, homozygous donor), starting frequencies set by the expected genotype frequencies of the generation, transitions from the average recombination rate over all generations of line development, and, after a change, a move to the other two states in proportion to their expected frequencies. Like theirs, our model does not need the donor's genotype, and its key parameter is chosen by a grid search (theirs: the non-informative rate and error rates, checked against SNP-chip calls on 24 lines; ours: the recombination scale ρ, chosen by likelihood). Their observations are GBS SNP calls on single lines; ours are k-mer scores of pools of six plants, which needs the pool-aware emissions below.

**The four parts of the model:**

| Part | What it is here |
|---|---|
| **Positions** | The 100-kb windows along one chromosome, in order. Each chromosome is run separately. Windows with fewer than 2,000 single-copy B73 k-mers, or where the B73 controls themselves see very few k-mers (repeats, gaps), are skipped. |
| **Hidden states** | The pool's state in the window: **B73** (all 6 plants AA), **SEG** (the family segregates: each plant AA, AB or BB), **DONOR** (all 6 plants BB). |
| **Observation** | r, the k-mer score of the window (≈ 1 for B73, ≈ d for DONOR, in between for SEG; see the worked example above). |
| **Emission** | For each state: how likely the observed r is if the pool really is in that state. |
| **Transition** | The chance that the state changes between one window and the next. |
| **Start / prior** | The share of the genome expected in each state: B73 84.4%, SEG 6.25%, DONOR 9.4% (BC2S3 expectation). |

#### 2.6.1 Emissions: how well does r fit each state?

Each state expects r near a certain value, with some spread (noise). The spread σ_w is measured, not assumed:
- **Pool noise σ:** how much r scatters around 1 in the pool's B73 windows. It is a *robust* standard deviation:
  1. For every B73-like window, take how far r is from 1: |r − 1|.
  2. Take the median of those distances. This is the **MAD** (median absolute deviation).
  3. Multiply by 1.4826. For normally distributed noise, 1.4826 × MAD equals the ordinary SD.

  Unlike the ordinary SD, the MAD isn't inflated by the few introgressed or junk windows, because a median ignores extreme values. Good pools have σ ≈ 0.08–0.10. Pools with σ > 0.15 are flagged **NOISY**, which mostly means low depth.
- **Window factor:** windows where the 10 B73 controls disagree with each other are noisier in every pool, so σ_w = σ × (control spread in this window ÷ typical control spread), limited to 0.5–3.

| State | Expected r | Emission (probability of r) | In words |
|---|---|---|---|
| **B73** | 1 | Normal(1, σ_w) | r near 1, with the pool's noise |
| **DONOR** | a × d_w | Normal(a·d_w, √(σ_w² + (a·τ_w)²)) | r near the donor level, with the pool's noise plus the uncertainty of the donor level |
| **SEG** | between d_w and 1 | average of Normal(a·d_w + (1 − a·d_w)·g(k/12), √(σ_w² + ((1 − g)·a·τ_w)²)) over k = 1…11 donor copies, weighted by Binomial(12, ½) | for each possible number of donor copies k in the 6 plants, r near the matching in-between level; the k are weighted by how likely they are under 1:2:1 segregation |

What the terms mean:
- **d_w:** the share of the window's B73 k-mers that the donor genome keeps, measured from the teosinte assemblies (step 19). Each window has its own value.
- **τ_w (uncertainty of d_w):** our donor plants are not the assembly plants, so their real level can differ from d_w. Teosinte is outcrossing and very diverse, so two plants of the same taxon can differ a lot in a given window; even the two Zv assemblies differ by more than 0.2 in 28% of windows. τ_w is how much the taxon's two assemblies disagree in that window (their SD × √1.5, the expected spread of a new plant around the mean of two, averaged over ±5 windows). It adds to the noise of DONOR, and to SEG in proportion to how much SEG depends on d. Mean τ: Zv 0.18, Zx 0.16, Zd 0.06. Zh and Zl have only one assembly, so τ = 0 for them. The B73 controls get the widest τ (Zv) as a worst-case test. Why it was added: in a first run without τ_w, Zv and Zx tracts kept flipping between SEG and DONOR; τ_w halved that. Details: [`final_script/README.md`, section 5.7](final_script/README.md#57-method-history-family-donor-levels-then-donor-level-uncertainty-current).
- **a:** one scale factor per pool. It absorbs genome-wide differences between our donor plant and the assembly plants; for example, our plant may keep a few percent more B73 k-mers everywhere. **How it's fitted (maximum likelihood by grid search):**
  1. The whole HMM is run for the pool with each value of a: 0.5, 0.6, …, 1.4, 1.6, 1.8.
  2. Each run gives the **likelihood**: the probability of the pool's observed r values in all windows of all chromosomes under that model. The forward algorithm computes it.
  3. The a with the highest likelihood is kept.

  The same grid is repeated for every value of **ρ**, the switch rate of the transitions (how easily the state changes from one window to the next; see [Transitions](#262-transitions-how-often-can-the-state-change)). a and ρ both change the likelihood, so the best a is found separately for each ρ, and then the best ρ is chosen over all pools.

  The fitted a is mostly 0.9–1.0 in every taxon, so on average the assembly levels fit our lines well.
- **k/12:** the donor fraction of a SEG pool. 6 plants carry 12 chromosome copies. A SEG family segregates 1:2:1, so the number of donor copies k follows Binomial(12, ½). The model averages over all k instead of assuming exactly half (see [What SEG means](#21-what-a-sample-is)).
- **f and g(f):** f is the donor fraction of the pool's chromosomes in the window (f = k/12). g(f) is how often the pool still finds the B73 k-mers that the donor lacks, compared with a pure-B73 pool. Those k-mers sit only on the B73 chromosomes, a share (1 − f) of the pool, so they get only that share of the reads. From the coverage formula, g(f) = (1 − e^(−(1−f)λ)) / (1 − e^(−λ)), where λ is the pool's k-mer coverage (measured from the share of B73 k-mers it finds in B73 windows). In the worked example, f = ½ gives g = 0.54. g = 1 when the pool is all B73 (f = 0), and g = 0 when it's all donor (f = 1).
- **Outliers:** each state also has a 2% chance that the window is junk. A junk window's r is treated as equally likely anywhere from 0 to 2 (a uniform distribution). This is a standard robustness device, a "contamination" term in the emission. Without it, one odd window (e.g. a mis-assembled or duplicated region) could give near-certain evidence for a state just because it lies far out in a normal curve's tail.
  - **The 2% is a chosen value, not estimated from the data.** The B73 controls show it's on the safe side (44-pool test, 10 controls, about 21,000 windows each):
    - 0.77% of control windows are more than 3σ from 1, where a perfect normal would give 0.27%.
    - 0.10% are more than 4σ from 1, where a normal would give 0.006%.
  - So real outliers exist, but at well under 2%. A higher value makes the HMM more reluctant to switch on a few windows; a lower one makes it switch more readily.

**Example: one Zd window** (d = 0.27, σ = 0.09, a typical pool that finds 30% of the k-mers, as in section 2.5; Zd's τ is small, so it's left out here). For a given r, this is how well each state fits on its own, before neighbouring windows are considered:

| r | B73 | SEG | DONOR | Reading |
|---|---|---|---|---|
| 0.27 | 0% | 1% | **99%** | clearly DONOR |
| 0.45 | 1% | 59% | 41% | between SEG and DONOR: the neighbours decide |
| 0.66 | 1% | **99%** | 0% | clearly SEG (about half the plants' chromosomes donor) |
| 0.85 | 48% | 52% | 0% | B73 or SEG with few donor copies: the neighbours decide |
| 1.00 | **97%** | 3% | 0% | clearly B73 |

SEG covers a wide range (here 0.35–0.95) because a 6-plant pool can carry from 1 to 11 donor copies.

**τ_w and a are only needed because we don't have the donor parents.** Both correct for the same thing: d_w comes from assemblies of *other* plants of the same taxon, not from the plants our lines came from. τ_w allows for differences in single windows; a allows for a difference over the whole genome. Both can be switched off:
- **τ_w:** on with `SPREAD=1` (`20_kmer_hmm.py --d_spread`, the default), off with `SPREAD=0`.
- **a:** fitted over `--scale_grid`; `--scale_grid 1` fixes a = 1.

**If the donor parents are sequenced,** d_w can be measured on the parent itself, in the same way as for an assembly in step 19 (the share of the window's B73 k-mers the parent has). Then τ_w = 0 and a = 1, and the emissions are simply:
- **B73:** r ~ Normal(1, σ_w)
- **DONOR:** r ~ Normal(d_w, σ_w)
- **SEG:** r ~ Σ_k w_k × Normal(d_w + (1 − d_w)·g(k/12), σ_w), summed over k = 1…11 donor copies, with w_k = Binomial(12, ½) probabilities rescaled to sum to 1

Each state keeps the 2% outlier term: probability = 0.98 × the Normal above + 0.02 × Uniform(0, 2).

One caveat: a teosinte parent is an outcrossing, heterozygous plant, and the F1 got only one of its two chromosome sets (recombined). The parent's k-mers are the union of both sets, so where its two sets differ, the parent's d_w is a bit higher than that of the set our line carries. A small τ_w can still help there, so it stays an option rather than being removed. Reading d_w from a sequenced parent (a `--parent` input) is not in the code yet; for now d_w comes only from assemblies.

#### 2.6.2 Transitions: how often can the state change?

Between two neighbouring windows, the chance of a state change is **s = ρ × θ**:
- **θ** is the recombination fraction between the two window centres, from the NAM genetic map (Haldane: θ = ½(1 − e^(−2·Morgans))). On average a 100-kb step is 0.075 cM, but it's near 0 around centromeres and higher near chromosome ends.
- **ρ** scales it. Our lines went through several meioses (BC1, BC2, three selfings), so segment ends are more frequent than one meiosis would give. ρ is not guessed: every pool is run with ρ = 0.25, 0.5, 1, 1.5, 2, 3, 4 and 8, and the value with the highest total likelihood over all donor pools (low-depth pools left out) is used. That is **ρ = 1.5**.
- **Where to:** if the state changes, it moves to one of the other two states in proportion to their prior (e.g. from B73, more often to DONOR (9.4%) than to SEG (6.25%)), as in Zhong et al. (2025).

With ρ = 1.5, the chance of a change is about **0.1% per window** on average. The model stays in the same state unless several windows in a row clearly fit another state. I checked this with the model itself, on a Zd DONOR tract at average recombination. A SEG piece inside it is called only when about **4 windows in a row (400 kb)** clearly look SEG (r = 0.66); less clear windows (r = 0.55) need 5. One noisy window can't break a tract, because the 2% outlier term caps its evidence.

#### 2.6.3 Calling the states

- **Forward–backward** computes, for every window, the probability of each state given the whole chromosome. The window gets the state with the highest probability.
- **Tracts:** adjacent windows with the same state are joined. Skipped windows take the state before them.
- **Plant genotypes** follow from the pool states: AA = B73 + SEG/4, AB = SEG/2, BB = DONOR + SEG/4.

#### 2.6.4 What is fitted, and where

| Quantity | Fitted per | How | Code |
|---|---|---|---|
| r (k-mer score) | pool × window | pool count ÷ fit of the B73 controls | `20_kmer_hmm.py` |
| σ (noise) | pool, then × window factor | robust SD of r in B73 windows; control spread | `20_kmer_hmm.py` |
| λ (k-mer coverage) | pool | share of B73 k-mers seen | `20_kmer_hmm.py` |
| d_w (DONOR level) | donor taxon × window | learned from all pools of the taxon, shrunk towards the assemblies' level (2.6.6); in switched windows, from the pools called DONOR | `new_model/20c_pool_levels.py`, `new_model/fix1/fx1_3_choose.py` (start: `19_kmer_donor_reference.sh`) |
| τ_w (uncertainty of d_w) | donor taxon × window | disagreement of the taxon's two assemblies; in switched windows, the spread between our DONOR pools | `new_model/fix1/fx1_kmer_hmm.py` (`--d_spread`) |
| Counted k-mers | donor taxon × window | all B73 single-copy k-mers, or only the informative ones where they raise the contrast (2.6.6) | `new_model/fix1/fx1_3_choose.py` |
| a (scale on d) | pool | best likelihood over 0.5–1.8 | `20_kmer_hmm.py` |
| ρ (switch rate) | whole data set | best total likelihood | `21_kmer_collect.py` |
| Prior | fixed | BC2S3 expectation | `config.sh` |
| Plants per pool | fixed | 6 | `config.sh` (`KMER_PLANTS`) |
| Outlier share | fixed | 2% | `20_kmer_hmm.py` (`--eps`) |

"Donor taxon" is the pool's teosinte (Zd, Zh, Zl, Zv or Zx, from the sample sheet). Every pool of that taxon uses the same d_w and τ_w. In the current model d_w is learned from all the pools of the taxon at once (2.6.6), without using any calls, except in switched windows; nothing here comes from our BC1 or BC2 plants. The HMM itself (emissions, transitions, decoding) is the same as in step 20: `new_model/fix1/fx1_kmer_hmm.py` is a copy of `20_kmer_hmm.py` with the two inputs of 2.6.6 added.

**Why the prior is fixed.** It only sets the state each chromosome starts in, and which state a switch goes to. It doesn't set how much of the genome ends up in each state: with about 2,000 windows per chromosome, the data decide. DONOR came out at 4.9% although the prior expects 9.4%, so the prior isn't forcing the result.

#### 2.6.5 How long are the tracts?

Tract length depends on how many meioses broke the donor chromosome up: here the F1, BC1 and BC2 plants, then two rounds of selfing before the BC2S2 mother of each pool. Each crossover cuts a donor segment, so after these meioses most segments are many Mb long, but a few are short. Short ones come mostly from regions with high recombination (near chromosome ends, where 500 kb can be more than 1 cM) and from crossovers in different generations falling close together.

To check that short tracts are expected and not just a feature of our data, the crosses were simulated on the NAM genetic map (F1 → BC1 → BC2 → BC2S1 → BC2S2, crossovers placed at random along the map, 400 lines per chromosome, no selection), and the lengths of the introgressed tracts (SEG and DONOR joined) were compared with our calls (1,266 pools):

| Introgressed tracts | Median length | < 0.5 Mb | < 1 Mb | < 2 Mb | < 5 Mb |
|---|---|---|---|---|---|
| Simulated BC2S2 mothers (3,506 tracts) | 11.3 Mb | 4.4% | 8.4% | 15.3% | 30.8% |
| Our pools, k-mer HMM (12,668 tracts) | 9.6 Mb | 2.6% | 6.7% | 13.8% | 33.7% |

- **Short tracts are expected from the crossing scheme alone**, not only in our data: about 1 in 12 introgressed tracts should be under 1 Mb.
- We find slightly fewer under 1 Mb than the simulation. That fits the HMM's detection limit: a piece needs about 4 windows (400 kb) to be called ([Transitions](#262-transitions-how-often-can-the-state-change)), so the shortest tracts are missed.
- **Any single short tract can still be an error**, e.g. a window run where our donor plant differs from the assemblies, or a structural variant. They are a small part of the genome: tracts under 1 Mb hold 0.18% of the introgressed sequence (0.11% in the simulation). The current model has slightly more of them than the first one (1.6 instead of 1.3 tracts under 1 Mb per pool; section 2.6.6). See [Limitations](#213-limitations-and-future-work) for how they can be checked.

The simulation does not include crossover interference (which makes two crossovers close together in one meiosis rarer) or selection, so it is a rough guide, not an exact expectation. Script: `final_script/tract_length_sim.py` (seed 1); output `data/introgression_fix1/checks/tract_length_sim.txt`.

#### 2.6.6 The current model: informative k-mers and DONOR levels from our pools

The first model took each window's DONOR level d_w from teosinte assemblies of *other* plants and counted all B73 single-copy k-mers. The targeted sequencing (section 2.9.1) found one systematic error with it: in *parviglumis* and *mexicana*, a SEG tract between two DONOR tracts was really DONOR (BB) in about half of the cases, and these errors recurred at the same places in many pools (chr6 21–23 Mb in 48 pools, chr1 297 Mb in 28, chr7 54–65 Mb in 28, chr3 75–78 Mb in 19). There were two causes, and the current model changes one input for each. The HMM itself is unchanged.

**1. A wrong DONOR level (chr6, chr1).** Where our donor plants keep more of the B73 k-mers than the assembly plants, a DONOR window reads too high and looks like SEG. Fix: learn d_w from **all pools of the taxon at once**, without using any calls (`new_model/20c_pool_levels.py`). In each window the scores of a taxon's 60–570 pools are a mixture of three groups that share the window's level: B73 (r ≈ 1), SEG (r = d + (1 − d)·g(f), f from Binomial(12, ½)) and DONOR (r ≈ d). d is fitted by maximum likelihood (the group sizes by EM), then shrunk towards the assembly level by how many pools carry teosinte there, so windows with little evidence keep the assembly level. This is the idea of RTIGER's shared emission fit (one set of emission parameters learned from all samples; Campos-Martin et al. 2023), applied to the parameter that was wrong. On simulated pools the fitted level was within 0.02 of the truth in a hotspot where the assembly level was 0.18 off, and unchanged elsewhere.

**2. Too little contrast (chr7, chr3).** Where the teosinte keeps most of B73's k-mers (d_w ≈ 0.7–0.9), DONOR and SEG score almost the same: the shared k-mers are present in every pool whatever its state, so they add counting noise and no signal. *parviglumis* chr7 54–65 Mb had a contrast (gap between the SEG and DONOR levels ÷ noise) of 1.2, i.e. the two states were about one noise unit apart. Fix: count only the **informative** k-mers there, the B73 single-copy k-mers found in **none** of the taxon's assemblies (`new_model/fix1/`):

![How fix 1 works](new_model/fix1/fix1_method.svg)

- Every pool was recounted against its taxon's informative k-mers (the 10 B73 controls against all five taxa's), with the same KMC + KCFtools steps as step 17.
- Per taxon and window, the scores of all pools were computed both ways, and a window **switches** to the informative k-mers only if B73 has at least 200 of them there and they raise the contrast at least 1.2-fold (contrast measured from the pools called DONOR and the B73 controls). Elsewhere, the old counts and the learned levels of point 1 are kept.
- The assemblies have none of the informative k-mers by construction, so in switched windows the DONOR level, and its spread between our donors (τ_w), come from our pools called DONOR there.
- 33% (*diploperennis*) to 53% (*huehuetenangensis*) of the windows where both contrasts could be measured switched. At *parviglumis* chr7 54–65 Mb, 95% of windows switched and the contrast rose from 1.2 to 2.6. (Contrasts after switching are somewhat optimistic, because windows were chosen where the measured new contrast was higher; the targeted sequencing below is the unbiased check.)

**Result.** The two fixes were tested one by one and together, each time with the whole pipeline (HMM, calls) and the targeted sequencing (section 2.9.1, the three bad samples of its point 6 left out), which is independent of our calls:

| Check | First model | Learned levels only | Informative k-mers only | **Both (current)** |
|---|---|---|---|---|
| SEG tracts that are AA / AB / BB in the plant | 17.9 / 74.3 / 7.8% | 6.2 / 81.8 / 12.1% | 5.7 / 86.4 / 7.9% | **5.7 / 85.8 / 8.5%** |
| DONOR tracts that are BB | 90.7% | 93.7% | 95.5% | **95.2%** |
| B73 tracts that are AA | 97.1% | 96.2% | 97.0% | **96.6%** |
| SEG between two DONOR tracts that is BB | 38 tracts, 47% | 38 tracts, 71% | 28 tracts, 54% | **30 tracts, 57%** |
| Pools with a SEG tract between two DONOR tracts at chr6 21–23 / chr1 297 / chr7 54–65 / chr3 75–78 Mb | 48 / 28 / 28 / 19 | 11 / 14 / 23 / 16 | 48 / 28 / 7 / 8 | **9 / 11 / 14 / 8** |
| Genome share SEG / DONOR | 6.22 / 4.94% | 6.22 / 4.88% | 6.26 / 4.93% | **6.23 / 4.94%** |
| Switch rate ρ chosen | 1.5 | 1.5 | 1.5 | **1.5** |

- **Better:** SEG calls where the plant carries no teosinte fell from 18% to 6%, DONOR calls that are BB rose to 95%, and the hotspot errors fell from 123 to 42 pools. B73 calls are as good as before (the 0.5-point drop is within the noise of about 1,000 tracts).
- **Unchanged:** the genome shares (so the donor deficit of 2.8.1 does not depend on how the DONOR level is set) and the switch rate.
- **Costs:** about 0.3 more short tracts (< 1 Mb) per pool (1.31 → 1.62), spread over the genome rather than in switched windows; the targeted sequencing compared only tracts of 1 Mb or more, so whether they are right is not known. And 26 more pools pass the noise cut-off σ > 0.15 (119 instead of 93 left out): 25 *parviglumis* pools that sat just under it (e.g. 0.149 → 0.157, as switched windows have fewer k-mers) and PN17_SID1585 (*mexicana*, 0.077 → 0.226, to be checked).
- **Not fixed:** SEG between two DONOR tracts is less common (0.38% → 0.30% of the genome), but 57% of the remaining ones are BB; *parviglumis* chr3 75–78 Mb is unchanged (its islands read between DONOR and SEG in every count, so they are not a contrast problem).

Details, simulations and all intermediate runs: [`new_model/README.md`](new_model/README.md), [`new_model/pilot_fix1/README.md`](new_model/pilot_fix1/README.md), [`new_model/fix1/README.md`](new_model/fix1/README.md).

### 2.7 How k-mers solve mapping bias, low depth and pooling

Back to the three problems of the SNP route ([above](#22-the-problem-mapping-reads-and-calling-snps-at-08x-in-pools)):

| Problem (above) | How the k-mer approach avoids it |
|---|---|
| **1. Donor reads are lost** | No mapping. A pool's raw reads are counted as k-mers, and only the presence of B73 k-mers is checked, so donor reads can't be filtered away. |
| **2. Too few reads per site** | No per-site calls. Each 100-kb window adds up ~32,000 k-mers, so its score r is measured precisely from thousands of positions at once. |
| **3. Pools aren't diploid** | The model is built for pools: DONOR = d, SEG = d + (1 − d)·g(f), averaged over the donor fractions a 6-plant family can have (Binomial(12, ½)). |

**Telling DONOR from SEG.** It no longer depends on allele ratios at single sites, which reference bias distorts. It depends on the window's **level**: in the Zd example above, DONOR ≈ 0.27 and SEG ≈ 0.66, with B73 = 1. The DONOR level is measured from teosinte genomes and then from our own pools (2.6.6), not guessed.

#### 2.7.1 Why not a SNP-based AA/AB/BB HMM (e.g. RTIGER)?

A common alternative is to call SNPs and run a three-state HMM whose states are plant genotypes, AA, AB and BB: for example RTIGER (Campos-Martin et al. 2023) or a basic three-state HMM. These are sound tools for what they were built for: **single plants** from a cross between **two known, sequenced parents**, where **every marker separates the two parents** (RTIGER was developed on segregants of two *Arabidopsis* accessions, and its markers are "differences between the two genotypes"). Our data break each of these assumptions:

1. **Pools are not plants.** AA/AB/BB describe one diploid plant. A pool from an AB mother carries k of 12 donor copies. With the donor share of 0.46 measured here (2.8.1), 84% of SEG pools sit at 4–8 of 12 copies, close enough to an "AB" level, but 12% sit at ≤ 3/12 and look like AA, and 4% sit at ≥ 9/12 and look like BB. More importantly, k changes **along** a SEG region, at every crossover of any of the six plants. A plant-genotype HMM reads each step as a possible state change, so one SEG region breaks into AB, AA and BB pieces: extra short tracts and breakpoints, less SEG, and false B73 inside SEG. Even where it is right, its "AB" means "the mother was AB" (our SEG), not heterozygous plants, so it cannot be read as plant heterozygosity.
2. **The markers do not separate the parents.** These models treat each SNP as diagnostic: reference allele = B73, alternative allele = donor. Our donors are outbred teosinte plants that were never sequenced. At many SNPs the donor carries the B73 allele or is heterozygous, and plants of one accession carry different haplotypes. Inside a real DONOR region these SNPs show B73 reads, which the model takes as AA or AB evidence, so DONOR is pulled towards AB. The targeted sequencing shows how common this is: even inside our DONOR tracts and at high depth, only about 1 real SNP in 4 carries a donor allele (25.5 per 100 SNPs; [2.9.1](#291-independent-check-targeted-sequencing-of-83-plots)). The k-mer model expects this: a DONOR window keeps a share d_w of its B73 k-mers, measured from teosinte assemblies (2.6.1).
3. **One read per site, and donor reads lost.** At 0.8x most SNPs have 0 or 1 read, so a single site cannot tell AB from AA or BB. The calls then rest on how "sticky" the HMM is (in RTIGER, the rigidity: the minimum number of consecutive markers that must support a change of state) and change when it is tuned. On top of that, mapping and its filters remove donor reads first: in a known DONOR tract only ~11% of the expected donor bases survive ([2.2](#22-the-problem-mapping-reads-and-calling-snps-at-08x-in-pools)), which pushes DONOR towards AB and AB towards AA.
4. **Priors and transitions of another population.** Unless they are set to BC2S3 expectations, the starting frequencies and switch rates are those of the population type the tool assumes (three states for an F2-like outcross, or two states for a backcross), not those of BC2S3 families.

**What this predicts.** Problems 2 and 3 pull DONOR towards AB (too little BB, AB inside true DONOR), and problem 1 breaks SEG regions into pieces (too many short tracts and breakpoints); the results also move with the stickiness setting. A deficit of homozygous donor or an excess of heterozygosity found this way can therefore be partly made by the method. The k-mer calls were checked against independent high-depth data instead (2.9.1): our B73 tracts are clean (96.5%), our DONOR tracts are homozygous teosinte (86%) and our SEG tracts are intermediate (donor allele fraction 0.45), so the deficit we report (2.8.1) is not a calling artefact.

**Testing a SNP-based call set.** The same check works for any call set on these pools: rename AA/AB/BB to B73/SEG/DONOR (they are then mother states) and run step 34d against the targeted sequencing.

### 2.8 Results

1,266 donor pools (the 119 noisy pools are left out; see [Noisy pools](#211-noisy-pools-and-resequencing)).

| Taxon | Pools | Donor alleles | SEG (pool) | DONOR (pool) | AB (plant) | BB (plant) | AB/(AB+BB) |
|---|---|---|---|---|---|---|---|
| Zd | 227 | 5.6% | 5.1% | 3.1% | 2.5% | 4.3% | 37% |
| Zh | 61 | 8.1% | 5.9% | 5.2% | 2.9% | 6.6% | 31% |
| Zl | 98 | 7.1% | 6.8% | 3.7% | 3.4% | 5.4% | 39% |
| Zv | 308 | 8.5% | 6.5% | 5.2% | 3.2% | 6.9% | 32% |
| Zx | 572 | 8.9% | 6.5% | 5.7% | 3.2% | 7.3% | 31% |
| **All** | **1,266** | **8.0%** | **6.2%** | **4.9%** | **3.1%** | **6.5%** | **32%** |
| *Expected BC2S3* | | *12.5%* | *6.25%* | *9.4%* | *3.1%* | *10.9%* | *22%* |

> **What AB/(AB+BB) means.** Take only the plants that carry at least one teosinte allele at a locus: these are the AB plants and the BB plants. AB/(AB+BB) is the fraction of them that are **AB** (heterozygous). The rest, BB/(AB+BB), are homozygous teosinte. Example: out of 100 plants that carry teosinte at a locus, how many are AB.
> - **Expected for BC2S3:** AB = 3.13% and BB = 10.94% of the genome, so AB/(AB+BB) = 3.13 / (3.13 + 10.94) = **22%**: about 22 AB and 78 BB plants out of 100.
> - **Observed:** AB = 3.12% and BB = 6.49%, so AB/(AB+BB) = 3.12 / (3.12 + 6.49) = **32%**: about 32 AB and 68 BB out of 100.
> - **Why it's useful:** it doesn't depend on how much teosinte there is in total. If lines lost teosinte at random (drift, or segments lost before BC2), AB and BB would drop by the same factor and the ratio would stay at 22%. It goes up only if BB plants are lost more than AB plants.

![Observed vs expected, all taxa](data/introgression_fix1/figures/fig1_overall.png)

The figure has three groups, because these numbers are different kinds of measures:
- **Pool states (SEG, DONOR):** what the HMM calls, per pool and window.
- **Plant genotypes (AB, BB):** what the pool states mean for single plants. They are **not measured** directly: a pool can't say which plant is AB. They are computed from the pool states using Mendel: a SEG region holds ¼ AA, ½ AB and ¼ BB plants, so **AB = SEG/2** and **BB = DONOR + SEG/4**. So AB lives inside SEG.
- **Donor alleles:** the total share of teosinte alleles, = DONOR + SEG/2.

![Observed vs expected by taxon](data/introgression_fix1/figures/fig2_by_taxon.png)

![Per family](data/introgression_fix1/figures/fig3_by_family.png)

![Per chromosome](data/introgression_fix1/figures/fig4_by_chromosome.png)

![Donor share per pool](data/introgression_fix1/figures/fig5_pool_share.png)

Figures from step 24 (`24_kmer_figures.py`); the per-family values are in `data/introgression_fix1/figures/family_values.tsv`. The same figures for the first model are in `final_script/figures/results_spread/` (`old_model_README.md`).

**What this shows:**
- **SEG is as expected:** the HMM calls 6.24% of the genome SEG, against 6.25% expected for BC2S3; the F1-family bootstrap finds no difference (p = 0.94). The AB estimate (3.12%) is half of SEG by construction, so it is the same result, not a second one.
- **DONOR is about half of the expectation** (4.93% vs 9.38% of the genome; estimated BB 6.5% vs 10.9%; p < 0.0002 in every taxon). That deficit is the whole shortfall in donor alleles (8.0% vs 12.5%). The first model gave the same (SEG 6.23%, DONOR 4.95%), so this does not depend on how the DONOR level is set.
- **So a larger share of the remaining donor is estimated to be het:** 32% of donor-carrying plants are AB, against 22% expected. This is not extra het: SEG, and so AB, is as expected. The share goes up only because DONOR went down, so AB is a bigger part of a smaller total. This holds on every chromosome (29–35%).
- **Why is DONOR low?** See [Why is there less donor DNA than expected?](#281-why-is-there-less-donor-dna-than-expected) below.
- **F1 families:** only 2 of 98 F1 families reach the expected 12.5% donor alleles (Zx.0390_P2, Zx.0130_P1); 86% sit between 5% and 11%. AB/(AB+BB) is above 22% in 87 of the 98 F1 families, and in all F1 families of Zd, Zh and Zl.
- **Small taxa:** Zd (6 F1 families), Zl (5) and Zh (2) have few F1 families, so their intervals are less reliable. Zd's SEG (5.1%) is somewhat below the expectation (p = 0.0004), and Zh's (5.9%, p = 0.0002) and Zl's (6.8%, p = 0.03) differ from it too, Zh's because its two families give a very narrow interval.

**Where in the genome the introgressions are.** Each pool's calls were put on a 100-kb grid (B73 = 0, SEG = 0.5, DONOR = 1: the share of donor alleles among the pool's plants), then averaged over all 1,266 donor pools and over 1-Mb bins. Each row is a chromosome; the colour is the share of donor alleles in that bin (8.0% on average, 12.5% expected for BC2S3):

![Donor allele share along the genome, all pools](figures/fix1/results/introgression_landscape.png)

- **Most of the genome carries a similar amount of teosinte:** the median 1-Mb bin has 7.7% donor alleles, and 90% of bins lie between 5.4% and 11.5%. Only 11 of the 2,138 bins reach the BC2S3 expectation of 12.5% (9 on chr2, 2 on chr7).
- **More teosinte:** chr2 (11.2% on average over the chromosome) and chr9 (10.8%). The most is at the start of chr2 (0–8 Mb: 13.9%), and it is high there in every large taxon (*diploperennis*, *parviglumis* and *mexicana* 14.2–14.9%).
- **Less teosinte:** chr5 (6.1%), chr4 (6.5%) and chr6 (6.6%). The least is near the end of chr5 (185–193 Mb: 4.1%, down to 3.8% in single bins) and in the middle of chr4 (100–120 Mb: 5.0%).

The same per taxon, in 5-Mb bins along the whole genome:

![Donor allele share along the genome, per taxon](figures/fix1/results/introgression_landscape_by_taxon.png)

- The taxa differ in where they carry teosinte. For example, *luxurians* has much more than average on chr9 (10–40 Mb: 16.2%), *huehuetenangensis* on chr3 (14.5% over the chromosome, against 8.7% in all pools), while *diploperennis* has almost none in chr5 185–193 Mb (0.3%) and chr4 100–120 Mb (0.5%), two regions where *parviglumis* still has about 7.7%.
- The small taxa come from few F1 families (*huehuetenangensis* 2, *luxurians* 5, *diploperennis* 6), so a block of colour in their rows can be one family's tract carried by many pools, not a property of the taxon.
- A region with little teosinte in every taxon could mean selection or distortion against teosinte alleles there, or a region where DONOR is hard to call; we have not tested this.

**Number and length of tracts per taxon.** A tract is a run of 100-kb windows in the same state. Introgressed tracts join adjacent SEG and DONOR tracts into one piece of donor DNA (as in [Tract lengths](#265-how-long-are-the-tracts)).

| Taxon | Pools | SEG tracts (per pool) | SEG median length | DONOR tracts (per pool) | DONOR median length | Introgressed tracts per pool | Introgressed median length |
|---|---|---|---|---|---|---|---|
| Zd | 227 | 1,117 (4.9) | 7.5 Mb | 788 (3.5) | 8.0 Mb | 6.2 | 10.9 Mb |
| Zh | 61 | 358 (5.9) | 6.5 Mb | 300 (4.9) | 8.8 Mb | 6.9 | 12.1 Mb |
| Zl | 98 | 589 (6.0) | 8.5 Mb | 355 (3.6) | 9.8 Mb | 7.2 | 12.2 Mb |
| Zv | 308 | 2,905 (9.4) | 4.9 Mb | 2,379 (7.7) | 7.8 Mb | 11.5 | 9.3 Mb |
| Zx | 572 | 5,594 (9.8) | 4.2 Mb | 4,759 (8.3) | 6.8 Mb | 11.5 | 9.2 Mb |
| **All** | **1,266** | **10,563 (8.3)** | **4.9 Mb** | **8,581 (6.8)** | **7.4 Mb** | **10.0** | **9.6 Mb** |

- A pool carries about 10 pieces of donor DNA, with a median length of about 10 Mb.
- *parviglumis* and *mexicana* pools carry nearly twice as many tracts as *diploperennis*, *huehuetenangensis* and *luxurians* pools (about 11.5 vs 6–7 introgressed tracts), and their SEG tracts are shorter (4.2–4.9 vs 6.5–8.5 Mb). Their donor allele share is not twice as high (table above), so they carry more, shorter pieces. We have not yet looked into why.

Figures and tables from `final_script/24b_introgression_landscape.py` (`introgression_landscape_bins.tsv`, `tract_stats_by_taxon.tsv` in `figures/fix1/results/`). Single samples can be browsed in the [interactive viewer](#212-interactive-figures-step-23).

The full statistics, F1 family results and the method history are in [`final_script/README.md`](final_script/README.md).

#### 2.8.1 Why is there less donor DNA than expected?

The lines carry 8.0% donor alleles instead of 12.5%. All of the shortfall is in DONOR (4.93% of the genome vs 9.38% expected), while SEG is as expected (6.24% vs 6.25%). Three things can lower the donor share, and the data can tell them apart in part:
- **Selection.** At every selfing step one plant per row was advanced (ear to row), and the plants with the best-looking ears were kept. The 6 plants pooled for sequencing were also the best-looking plants of their plot.
- **Segregation distortion or lower survival** of seeds and seedlings carrying teosinte alleles.
- **Errors in the calls.** If DONOR windows were called B73 (e.g. where the donor keeps most of B73's k-mers, d_w close to 1), DONOR would fall without SEG changing.

**Check 1: how many donor copies do SEG pools carry?** A SEG pool comes from an AB mother. Under Mendel, half of its plants' 12 gene copies are teosinte (donor allele share p = ½). The SEG level of r depends on p: more donor copies means fewer B73 k-mers (r = d_w + (1 − d_w) g(f), section 2.6.1). For each pool, p was found so that the mean r expected over its SEG windows matches the mean r observed (method of moments, with the pool's own d_w as used by the HMM (2.6.6), scale a and k-mer coverage). Three checks come with it: a simulation of SEG windows at p = ½ with each pool's noise (the method should return ½), the mean r in DONOR windows against a·d_w, and the mean r in B73 windows against 1.

| Taxon | Pools | Donor allele share in SEG pools, median (95% CI) | DONOR level: r − a·d_w | B73 level: r − 1 | Simulation at p = ½ |
|---|---|---|---|---|---|
| Zd | 222 | 0.464 (0.455–0.472) | +0.007 | −0.001 | 0.500 |
| Zh | 61 | 0.484 (0.456–0.523) | +0.002 | +0.001 | 0.502 |
| Zl | 96 | 0.453 (0.443–0.472) | +0.002 | +0.001 | 0.501 |
| Zv | 308 | 0.451 (0.443–0.457) | +0.038 | +0.002 | 0.500 |
| Zx | 562 | 0.453 (0.446–0.458) | +0.033 | 0.000 | 0.500 |
| **All** | **1,249** | **0.455 (0.451–0.458)** | **+0.029** | **0.000** | **0.500** |

- **The method works:** the simulation returns ½, and the B73 level is where the model puts it. The DONOR level is close (offsets ≤ 0.04); in *parviglumis* and *mexicana*, DONOR windows read about 0.03–0.04 above a·d_w.
- **It is not an edge effect.** Near the ends of a SEG tract, stretches with 10–11 donor copies could be called DONOR and stretches with 1–2 copies B73, which would move the estimate. Rerun on the middle of long SEG tracts only (tracts ≥ 5 Mb, 1 Mb trimmed at each end; 1,207 pools, 1.36 million windows; `--core_mb 5`), it is unchanged: 0.454 (0.450–0.460) overall and 0.448–0.493 per taxon.
- **SEG pools carry 45–46% donor alleles, not 50%,** alike in all five taxa (the first model gave 0.459). Part of this can come from the DONOR offset (a DONOR r above a·d_w pulls the estimate down, by up to about 0.03–0.04 in *parviglumis* and *mexicana*), but *diploperennis* and *luxurians*, with offsets of 0.002–0.007, are still at 0.45–0.46.
- **This fits the choice of the best-looking plants.** Within a SEG plot the plants are AA, AB and BB; choosing the best-looking ones for the pool favours plants without the teosinte form of a trait, i.e. fewer BB, so fewer donor copies. Seeds or seedlings with teosinte alleles doing worse would give the same result; the two can't be told apart here.

**Check 2: what would a loss like this do over the whole crossing scheme?** If every generation passed on teosinte alleles at 0.46 instead of 0.5 (F1 → BC1 → BC2 → two selfings to the BC2S2 mother):

| | SEG | DONOR | Donor alleles |
|---|---|---|---|
| Mendel (0.5) | 6.25% | 9.38% | 12.5% |
| Transmission 0.46 in every generation | 5.22% | 6.70% | 9.10% |
| **Observed** | **6.24%** | **4.93%** | **8.0%** |

- A loss of this size explains most of the missing donor alleles (about 3.4 of the 4.5 points).
- But it gets the split wrong: it should lower SEG too (to about 5.2%), while SEG is on expectation and DONOR is even lower. So something acts **against BB more than against AB**.
- That is what choosing the best-looking plant to advance does, if teosinte effects are mostly recessive:
  - **Selecting against teosinte.** A BB plant shows the teosinte form of the traits its segment controls (e.g. smaller or harder ears), with no B73 copy to mask it, so it is less often chosen. An AB plant still has a B73 copy, which often hides the teosinte effect.
  - **Selecting for vigour.** An AB plant is a small "hybrid" between B73 and teosinte at that segment and can be more vigorous than either homozygote; picking the best-looking ears is in part picking for that vigour.
  - **Hidden harmful alleles.** Teosinte is outcrossing, so it carries recessive harmful alleles that do no harm in AB plants but show up in BB plants (the same reason maize inbreeding reduces vigour).

**Check 3: comparison with TeoNAM.** TeoNAM (Chen et al. 2019; 1,257 BC1S4 lines) is a maize × teosinte population made the same way (backcross, then selfing). It also kept less teosinte than Mendel predicts: 15% of SNP calls were homozygous teosinte (14.2–16.2% across its five subpopulations) against 23.4% expected, and 8.1% were heterozygous against 3.1% expected. The authors' explanation is the one we reach above: "unconscious selection for more heterozygous plants", and, at the region with the highest heterozygosity near *tga1*, selection against plants homozygous for the teosinte allele, which "have poor ear and kernel quality". Our lab's earlier lines also showed more heterozygosity than expected.

TeoNAM is independent support, not a yardstick: its numbers cannot be compared directly with ours.
- **Low-coverage method.** TeoNAM lines were genotyped by GBS (low-coverage sequencing at selected SNP sites), and missing calls were filled in by imputation (FSFHap in TASSEL5), which assumed a heterozygosity of 3.125%. We did not call SNPs or impute: each 100-kb window gets its state from its own k-mer counts, and the HMM links neighbouring windows.
- **Not a genome share.** Their 8.1% heterozygous figure is "across all SNP sites", and those sites are not spread evenly over the genome. Ours are shares of the genome (100-kb windows along every chromosome).
- **Uncertain heterozygous calls.** At low depth a heterozygous site often shows reads from only one allele (heterozygotes missed), while sequencing errors and repeats can create false ones, so their heterozygous share is uncertain in both directions.
- **Different population.** W22 instead of B73, inbred teosinte instead of outbred teosinte plants (which may carry more harmful recessive alleles), and BC1S4 single lines instead of BC2S3 pools.

To still ask what TeoNAM's numbers would mean for our pools, a calculation note from our lab (18 November 2026) fits a simple model to them. The model has three numbers, each 0.5 under Mendel:
- t: the share of the genome that is heterozygous after a backcross;
- k: the share of heterozygous genome that stays heterozygous after one selfing;
- q: of the heterozygous genome that becomes fixed, the share fixed for teosinte.

Two observed values cannot set three numbers, so one is kept at 0.5 in turn. This gives two bounding cases: **selection during selfing** (t = 0.5, k = 0.63, q = 0.36) and **selection at the backcross** (t = 0.38, k = 0.68, q = 0.5). For a pool, the expected donor share inside a heterozygous (SEG) segment is (1 − k)·q + k/2; it falls below ½ only when q < ½, i.e. when selfing fixes B73 more often than teosinte.

| Expected for our BC2S3 pools | DONOR (homozygous) | SEG (heterozygous) | Donor alleles | Donor share inside SEG pools |
|---|---|---|---|---|
| Mendel | 9.38% | 6.25% | 12.5% | 0.500 |
| From TeoNAM, selection during selfing | 5.4% | 10.0% | 9.8% | 0.447 |
| From TeoNAM, selection at the backcross | 3.9% | 6.6% | 7.2% | 0.500 |
| **Observed** | **4.93%** | **6.24%** | **8.0%** | **0.455** |

- DONOR and the donor allele share fall inside the range from TeoNAM.
- The SEG range rests almost entirely on TeoNAM's 8.1% heterozygous figure, the least certain of its numbers. If their true heterozygosity were lower, the range would move down to ours:

  | TeoNAM heterozygosity used | Expected SEG for our pools |
  |---|---|
  | 8% (as reported) | 6.6–10.0% |
  | 6% | 5.3–8.7% |
  | 5% | 4.6–7.9% |
  | **Observed** | **6.24%** |

- The donor share inside SEG pools (0.455) is closer to the "selection during selfing" case (0.447) than to the "selection at the backcross" case (0.500). The model uses one average for the whole genome; the truth is probably a mix of the two cases.

**Check 4: are the calls right?** Independent data say yes. 83 of our plots were also genotyped by high-depth targeted sequencing, and every tract of ours (≥ 1 Mb) was compared with it (details in [2.9.1](#291-independent-check-targeted-sequencing-of-83-plots)). **95.9% of our B73 tracts carry no teosinte (97.5% of their length), and only 0.8% are homozygous teosinte**; 89% of our DONOR tracts are homozygous teosinte. So the low DONOR share is not DONOR called as B73: the teosinte really is missing. Inside our SEG tracts the targeted sequencing gives a median donor allele fraction of 0.45, close to the 0.46 of Check 1.

**So:** selection for good-looking plants, at two levels, explains the pattern: choosing the plant to advance removes BB mothers (low DONOR), and choosing the best-looking plants for the pool lowers the donor share inside SEG pools (0.45–0.46). TeoNAM's authors explained their excess of heterozygous and maize genome the same way, and the targeted sequencing (Check 4) rules out calling errors as the cause. Lower survival of teosinte-carrying seeds and seedlings may add to it.

**Caveats.**
- The transmission calculation assumes the same loss in every generation and in both sexes; it is a rough guide.
- Call errors (DONOR called B73) are ruled out for tracts of 1 Mb or more (Check 4); shorter tracts were not checked.
- Because the pooled plants were chosen, a SEG pool's donor copies are not exactly Binomial(12, ½) as the HMM assumes. The shift (0.455 vs 0.5) is small compared with the spread the SEG emission already allows.

Script: `final_script/22b_seg_donor_copies.py` (Slurm: `22b_seg_donor_copies.sh`, which passes the model's DONOR levels with `MODEL=fix1`; per-pool and per-taxon tables: `data/introgression_fix1/seg_donor_copies*.tsv`).

### 2.9 Validation

- **B73 controls are clean:** all 10 have ≤ 0.56% non-B73 windows (donor share ≤ 0.28%). This is a fair test, because the k-mer route never used the controls to filter sites.
- **Known tracts:** they are found with the expected levels. For example, a Zd DONOR tract on chr2 has r ≈ 0.20 against an assembly value of 0.26 (first model).
- **Agreement with the SNP route:** donor shares agree where the SNP route works. Where they differ, the cause is on the SNP side: lost donor reads, κ collapsing, or noisy pools.
- **A second k-mer measure (section 2.10):** teosinte-specific k-mers in the reads give the pool's donor taxon (matches the pedigree in 93.9% of pools, 99.7% when the assignment is clear) and its donor share, which agrees with the HMM's across pools (Spearman ρ 0.91–0.98 within taxa, in both models). It uses other k-mers (teosinte k-mers present, not B73 k-mers missing) but the same reads, so it is complementary, not fully independent. It checks the total donor share, not the SEG/DONOR split.

- **Targeted sequencing of 83 plots (section 2.9.1):** an independent, high-depth genotype for the same plots agrees with our calls: 95.9% of our B73 tracts carry no teosinte, 89% of our DONOR tracts are homozygous teosinte, and 83% of our SEG tracts are heterozygous.

#### 2.9.1 Independent check: targeted sequencing of 83 plots

**Data.** 96 of our plots were also genotyped by targeted sequencing: a fixed panel of SNP regions, sequenced at high depth (paired-end reads; median 23 reads per SNP in our comparison). The samples are named after the plot (e.g. `Zd0010_P3_P3_P3211` = plot `Zd.0010_P3_P3_P3.2.1.1`), and all 96 match a plot of ours. 95 samples had both read files, 94 aligned, 86 belong to a pool with k-mer calls in the current model (88 in the first; two pools are now among the noisy pools), and **83** had enough SNPs to compare. The DNA is thought to come from one plant of each plot, but this is not yet confirmed (see "One plant or a bulk?" below).

**Methods.**
1. **Reads to SNPs** (steps 34, 34b). Adapters trimmed (fastp), reads aligned to B73 v5 (bwa mem; duplicates kept, since the panel may be amplicon-based), and variable sites found in all samples together (bcftools mpileup with mapping and base quality ≥ 20, multi-sample call, biallelic SNPs with QUAL ≥ 30): 656,483 SNPs. Kept: sites with ≥ 10 reads in ≥ 20 samples, and not carrying the alternative allele in more than half of the samples (at any one place most pools are B73, so such sites are paralogs or places where our B73 differs from the reference): **84,813 sites**.
2. **Donor allele fraction.** At each site, a sample's donor allele fraction is alt reads / all reads; the reference is B73, so the alternative allele is the non-B73 one. The donor genotypes are not needed. A site counts as **carrying donor alleles** if it has ≥ 3 alt reads and a fraction ≥ 0.2.
3. **Real SNPs only.** Most of the 84,813 sites are not real SNPs (see Results 1). A site is used for a sample only if it is homozygous for the alternative allele (fraction ≥ 0.9, ≥ 10 reads) in at least one *other* sample. Leaving the sample out means its own reads do not decide which of its sites are used. **7,269 real SNPs** remain (7,666 with the 85 samples of the first model).
4. **Our call at each site** (step 34c): the k-mer state (B73 / SEG / DONOR) of the sample's pool at that position.
5. **Tract comparison** (step 34d). Each tract of ours of 1 Mb or more is compared with the targeted-sequencing reads inside it. The sample's own noise (sites carrying donor alleles inside our B73 tracts; median 0.4 per 100 real SNPs) sets what counts as signal:
   - **carrier:** at least 2 sites carry donor alleles, and more than the noise allows (Poisson, p < 0.001); then **BB** if their median fraction is ≥ 0.8, **AB** if it is 0.25–0.75;
   - **AA:** not a carrier, with at least 15 sites covered;
   - **unclear:** otherwise (too few sites).

   A DONOR region whose teosinte carries the B73 allele at most panel SNPs can look AA; the AA share of our DONOR tracts measures this.
6. **Missed introgressions.** Inside our B73 tracts, 2-Mb windows (1-Mb steps) with more donor-carrying sites than the sample's noise allows (Poisson, p < 0.001).

**Results.**

![Targeted sequencing against the k-mer calls](figures/fix1/targseq_validation.png)

*a: sites carrying donor alleles per 100 SNPs, inside our B73, SEG and DONOR calls, for all sites and for real SNPs only. b: genotype of our tracts in the targeted sequencing. c: median donor allele fraction of each tract that carries donor alleles.*

| Our call | Tracts ≥ 1 Mb | Clearly called | AA (no teosinte) | AB (intermediate) | BB (homozygous teosinte) |
|---|---|---|---|---|---|
| **B73** | 1,471 | 1,013 | **95.9%** (97.5% of Mb) | 3.4% | 0.8% |
| **SEG** | 736 | 329 | 8.8% | **82.7%** (92.6% of Mb) | 8.5% |
| **DONOR** | 761 | 377 | 7.2% | 4.0% | **88.9%** (90.4% of Mb) |

All 83 samples. First model (85 samples): B73 96.5 / 2.9 / 0.6%, SEG 21.2 / **71.1** / 7.7%, DONOR 9.1 / 5.0 / **85.9%**. Without the three bad samples of point 6, the current model gives B73 96.6 / 3.0 / 0.4%, SEG 5.7 / 85.8 / 8.5%, DONOR 1.7 / 3.1 / 95.2% (the comparison of the models in 2.6.6 uses this set).

1. **91% of the targeted-sequencing "SNPs" were amplification errors.** Of the 84,813 sites, only 7,666 are homozygous for the alternative allele in any sample. The others carry the alternative allele at a fraction of about 0.3, at random sites, in single samples; they are not shared by relatives and do not come from another sample in the batch. This fits errors made early in the amplification, which are then copied. On all sites, about 1 in 100 sites inside our B73 tracts looked like donor alleles, and a first comparison called 74% of our B73 tracts "AB". On real SNPs the contrast is clean: donor alleles at 0.6 per 100 SNPs in our B73 tracts, 16.5 in SEG and 25.3 in DONOR (figure, a).
2. **Our B73 calls are right.** 95.9% of our B73 tracts (97.5% of their length) show no teosinte, and 0.8% are homozygous teosinte. So DONOR is not being called B73, and the low DONOR share (section 2.8.1) is real.
3. **Our DONOR calls are right.** 89% of our DONOR tracts (90% of their length) are homozygous teosinte (86% in the first model), with a median donor allele fraction of **1.00** in the carrier tracts (figure, c); the targeted sequencing has no bias towards the B73 allele. Most of the 7% that look AA (21 of 27) come from the three bad samples in point 6; without them, 95% of DONOR tracts are BB.
4. **Our SEG calls are heterozygous.** 83% of our SEG tracts are AB (71% in the first model), with a median donor allele fraction of **0.45** (300 carrier tracts; figure, c). SEG tracts where the plant carries no teosinte fell from 21% to 9%. This is close to the 0.455 measured from the k-mers inside SEG pools (Check 1 in 2.8.1), so two independent methods give the same number.
5. **Missed introgressions are rare.** Inside our B73 tracts, 176 of 20,730 windows of 2 Mb (in 55 samples) carry more donor alleles than noise allows, against about 21 expected by chance. That is under 1% of the B73-called genome. Their median donor allele fraction is 0.37, so they are short intermediate (SEG-like) pieces at most, not missed DONOR.
6. **Three samples are bad and explain most mismatches:**
   - `Zv0480_P3_P2_P1211` is **contaminated**: 53% of its stray donor alleles (inside our B73 tracts) are homozygous-donor sites of `Zx0550_P4_P3_P5311`, a *mexicana* sample. This explains why its DONOR regions show intermediate fractions (mean 0.44).
   - `Zx0090_P2_P5_P1111` and `Zx0160_P2_P4_P1111` look like **sample swaps or wrong plants**: 8 and 10 of their DONOR tracts show no teosinte at all. Their targeted-sequencing DNA does not seem to come from the pool we sequenced.
7. **The remaining error: some SEG tracts between two DONOR tracts are DONOR.** In the targeted sequencing (three bad samples left out), SEG tracts that lie between two DONOR tracts are BB in **57%** of cases (17 of 30 clearly called), against 6.5% for SEG tracts next to one DONOR tract and 2.1% for SEG tracts with no DONOR neighbour (4 of 195). This happens only in *parviglumis* (4 of 8) and *mexicana* (13 of 19); for *diploperennis*, *huehuetenangensis* and *luxurians* none of these tracts is BB. In the first model the same check gave 47% (18 of 38): the current model removed 8 of these tracts, but mostly correct ones, so the share of errors among the rest went up. These tracts are now 0.30% of the genome (4.8% of all SEG; median 1.6 Mb), against 0.38% before. If about half of them are DONOR, DONOR rises from 4.93% to about 5.1% and SEG falls from 6.24% to about 6.1%, which does not change the DONOR deficit in 2.8.1.

   **The errors sit at the same places in many pools.** Counting, for every 1-Mb bin, the pools whose SEG tract between two DONOR tracts covers it, the first model had four hotspots: chr6 21–23 Mb (48 pools), chr1 297 Mb (28), chr7 54–65 Mb (28) and chr3 75–78 Mb (19). This is what led to the current model (2.6.6), which brings them to 9, 11, 14 and 8 pools:
   - **chr6 and chr1: the DONOR level was wrong.** There, our donors keep more of the B73 k-mers than the assembly plants, so a DONOR window read like SEG. The DONOR level learned from all our pools removes most of these errors.
   - **chr7 and chr3: too little contrast.** There, the teosinte keeps 61–78% of the B73 k-mers, so with all B73 k-mers counted a DONOR window scores almost like a SEG window (about one noise unit apart in *parviglumis* chr7). Counting only the informative k-mers raised the contrast (*parviglumis* chr7: 1.2 → 2.6) and turned most chr7 islands into DONOR (e.g. `Zv0210_P1_P5_P1111` chr7 at 54–66 Mb, in the gallery below). ***parviglumis* chr3 75–78 Mb is not fixed:** its islands read halfway between DONOR and SEG in every count, all B73 k-mers or informative only, so the problem there is not contrast; it may be family-specific teosinte haplotypes that share these k-mers with B73. The targeted sequencing at chr3 70–81 Mb is mixed (BB in some plots, AA or AB in others; most of them *mexicana*).
   - None of the QTL that hold in the joint model lies in these four regions.

   **Contrast: how far apart DONOR and SEG are in one window.** The HMM separates DONOR from SEG by the gap between their expected scores, about (1 − d_w) × g(½), where d_w is the share of B73 k-mers the teosinte keeps in the window. Divided by the window's noise, this gives the contrast: how many noise units apart the two states are. Measured from our pools (windows with at least 5 pools called DONOR; median, and share of windows below 3):

   | Taxon | All B73 k-mers | Current model (switched windows: informative k-mers) |
   |---|---|---|
   | *diploperennis* | 4.2 (9% below 3) | 4.7 (3%) |
   | *luxurians* | 4.3 (15%) | 4.8 (7%) |
   | *mexicana* | 4.1 (19%) | 4.7 (5%) |
   | *huehuetenangensis* | 3.3 (41%) | 4.0 (15%) |
   | *parviglumis* | 3.3 (41%) | 3.8 (21%) |

   *parviglumis*, the closest relative of maize, shares the most k-mers with B73, so its DONOR and SEG are still the hardest to tell apart (the "current" column is somewhat optimistic, since windows switched where the measured new contrast was higher). The contrast is per window: the HMM pools the evidence of a whole stretch, so low contrast makes DONOR vs SEG harder, not impossible. **B73 vs introgressed is not affected**, and with it the ancestry dosage used in the GWAS except for the ½-vs-1 step.

**Examples: our calls against the targeted sequencing.** "Ours" is our k-mer call; "theirs" is the targeted sequencing. Each panel is one sample and one chromosome, in three rows:
- **OUR call:** our HMM call of the 6-plant pool, from the 0.8x skim-seq: B73 (grey), SEG (blue) or DONOR (orange);
- **THEIR call:** the targeted sequencing's call for each targeted region (SNPs within 200 kb of each other): **AB** (blue) if its donor-carrying SNPs sit below 0.8, **BB** (orange) if at 0.8 or above. A short light tick means no donor allele was seen there: either the plant is AA, or the teosinte carries the B73 allele at those SNPs, so the region says nothing;
- **THEIR SNPs:** the donor allele fraction at each real SNP, the raw data behind the second row. Dots come in clumps: one clump is usually one amplicon (10–20 SNPs within about 150 bp, read from the same reads), so one odd amplicon gives a clump of odd dots and its own small box in the second row.

**The text above our row is the comparison**, one per tract of ours that the targeted sequencing calls clearly: "ours SEG → theirs AB" means we call the pool SEG there and the plant is AB over that tract (step 34d, from all its SNPs, not from single dots or boxes). Short tracts get the short form "SEG→AB". What to expect: our B73 → theirs AA and our DONOR → theirs BB, because every plant of such a plot must carry that genotype. Our SEG means the plot segregates, so one plant of it can be AA, AB or BB (1 : 2 : 1); a single "ours SEG → theirs AA" or "→ theirs BB" is therefore not an error. Only the proportions over all SEG tracts can be checked (point 4, and "One plant or a bulk?" below).

The outlined tract is the one in the title. **Four good examples:**

![Our k-mer calls against the targeted sequencing, four good examples](figures/fix1/targseq_examples.png)

1. **Ours SEG, theirs AB** (`Zl0050_P1_P5_P1211`, chr5 0–193 Mb). Every amplicon with donor alleles in the tract sits at 0.3–0.6 (104 donor-carrying SNPs, median 0.42).
2. **Ours DONOR, theirs BB** (`Zv0120_P2_P3_P2111`, chr7 0–92 Mb). Every amplicon with donor alleles sits at about 1 (24 SNPs, median 1.0), and the B73 rest of the chromosome has none.
3. **Ours SEG then DONOR, theirs AB then BB** (`Zh0030_P4_P3_P1111`, chr6). SEG at 8–75 Mb is AB (median 0.50), DONOR at 75–142 Mb is BB (median 1.0).
4. **Ours SEG; theirs between AB and BB along the whole tract** (`Zv0480_P3_P1_P1211`, chr5 11–150 Mb). The 23 amplicons with donor alleles sit anywhere from 0.3 to 0.9, most at 0.6–0.9 (99 SNPs, median 0.78, so step 34d calls the tract "unclear"). One plant cannot switch between AB and BB this often within one tract, and this sample is clean (no stray donor alleles in its B73 tracts, unlike the contaminated `Zv0480_P3_P2_P1211` of point 6). A fraction of about 0.7 all along fits DNA of several plants of a segregating plot, e.g. 6 plants carrying 8–9 of their 12 copies as teosinte. Our SEG call fits that plot. This example is evidence that at least some targeted-sequencing DNA is not one plant (see "One plant or a bulk?" below).

**Four bad examples** (the three bad samples of point 6 left out):

![Our k-mer calls against the targeted sequencing, four bad examples](figures/fix1/targseq_examples_bad.png)

1. **Ours DONOR, theirs AB** (`Zd0010_P3_P4_P2311`, chr3 28–157 Mb). 16 amplicons with donor alleles, most at 0.3–0.7 (38 SNPs, median 0.60; only 7 at 0.8 or above). The plot most likely segregates here and our pool read too high, e.g. a pool that by chance held few AA plants (2.9.1, point 4). This is the largest DONOR → AB tract of the 80 good samples.
2. **Ours SEG between two DONOR tracts, theirs BB** (`Zv0030_P1_P3_P1111`, chr1, zoomed to 280–309 Mb). The systematic error of point 7, at the chr1 297 Mb hotspot: our 2.2-Mb SEG tract has SNPs at about 1, like the DONOR beside it.
3. **Our tract ends too early** (`Zv0370_P3_P4_P1111`, chr1, zoomed to 270–309 Mb). Our SEG tract ends at 297.0 Mb; the targeted sequencing finds an AB amplicon at 298.6 Mb (11 SNPs, median 0.41), inside what we call B73. The donor tract probably goes on for at least 1.5 Mb. This is a boundary error near the same chr1 hotspot.
4. **A disagreement made by the check itself** (`Zv0420_P3_P1_P1111`, chr9, zoomed to 0–60 Mb). We call DONOR at 1.4–21 Mb. Two amplicons say BB (6.3 Mb: 18 SNPs at 0.98; 10.8 Mb: 2 SNPs at 0.90) and one says AB (18.7 Mb: 54 SNPs at 0.47). Because step 34d takes the median over SNPs, the one dense amplicon outvotes the other two, and the tract counts as DONOR → AB. Counted per amplicon it would be BB at two of three places. So part of the "ours DONOR → theirs AB" share (3%) can come from single amplicons with many SNPs, not from our calls.

<details><summary>All 83 samples: every chromosome with an introgression (one figure per sample)</summary>

Same three rows and the same "ours → theirs" labels as above, for every chromosome where we call SEG or DONOR or the targeted sequencing finds donor alleles. The three bad samples of point 6 are included and marked BAD SAMPLE in their titles: `Zv0480_P3_P2_P1211` (contaminated), `Zx0090_P2_P5_P1111` and `Zx0160_P2_P4_P1111` (swapped). Made with `final_script/34f_targseq_examples.py --gallery_dir`.

![Zd0010_P3_P3_P3211](figures/fix1/targseq_all/Zd0010_P3_P3_P3211.png)
![Zd0010_P3_P3_P3411](figures/fix1/targseq_all/Zd0010_P3_P3_P3411.png)
![Zd0010_P3_P4_P2311](figures/fix1/targseq_all/Zd0010_P3_P4_P2311.png)
![Zd0020_P3_P4_P2111](figures/fix1/targseq_all/Zd0020_P3_P4_P2111.png)
![Zd0020_P3_P4_P3211](figures/fix1/targseq_all/Zd0020_P3_P4_P3211.png)
![Zd0030_P2_P5_P5211](figures/fix1/targseq_all/Zd0030_P2_P5_P5211.png)
![Zd0040_P1_P4_P3111](figures/fix1/targseq_all/Zd0040_P1_P4_P3111.png)
![Zh0020_P3_P5_P1211](figures/fix1/targseq_all/Zh0020_P3_P5_P1211.png)
![Zh0030_P4_P1_P2211](figures/fix1/targseq_all/Zh0030_P4_P1_P2211.png)
![Zh0030_P4_P3_P1111](figures/fix1/targseq_all/Zh0030_P4_P3_P1111.png)
![Zl0010_P1_P2_P1411](figures/fix1/targseq_all/Zl0010_P1_P2_P1411.png)
![Zl0010_P1_P4_P1411](figures/fix1/targseq_all/Zl0010_P1_P4_P1411.png)
![Zl0020_P3_P1_P1411](figures/fix1/targseq_all/Zl0020_P3_P1_P1411.png)
![Zl0020_P3_P4_P3511](figures/fix1/targseq_all/Zl0020_P3_P4_P3511.png)
![Zl0040_P4_P4_P2411](figures/fix1/targseq_all/Zl0040_P4_P4_P2411.png)
![Zl0050_P1_P5_P1211](figures/fix1/targseq_all/Zl0050_P1_P5_P1211.png)
![Zv0030_P1_P3_P1111](figures/fix1/targseq_all/Zv0030_P1_P3_P1111.png)
![Zv0080_P2_P1_P1211](figures/fix1/targseq_all/Zv0080_P2_P1_P1211.png)
![Zv0120_P2_P3_P2111](figures/fix1/targseq_all/Zv0120_P2_P3_P2111.png)
![Zv0150_P4_P1_P1111](figures/fix1/targseq_all/Zv0150_P4_P1_P1111.png)
![Zv0190_P3_P4_P1111](figures/fix1/targseq_all/Zv0190_P3_P4_P1111.png)
![Zv0200_P3_P2_P1111](figures/fix1/targseq_all/Zv0200_P3_P2_P1111.png)
![Zv0210_P1_P4_P2211](figures/fix1/targseq_all/Zv0210_P1_P4_P2211.png)
![Zv0210_P1_P5_P1111](figures/fix1/targseq_all/Zv0210_P1_P5_P1111.png)
![Zv0220_P3_P1_P2211](figures/fix1/targseq_all/Zv0220_P3_P1_P2211.png)
![Zv0220_P3_P2_P1111](figures/fix1/targseq_all/Zv0220_P3_P2_P1111.png)
![Zv0250_P3_P2_P1311](figures/fix1/targseq_all/Zv0250_P3_P2_P1311.png)
![Zv0260_P2_P2_P1111](figures/fix1/targseq_all/Zv0260_P2_P2_P1111.png)
![Zv0310_P3_P2_P1311](figures/fix1/targseq_all/Zv0310_P3_P2_P1311.png)
![Zv0310_P3_P3_P1211](figures/fix1/targseq_all/Zv0310_P3_P3_P1211.png)
![Zv0370_P3_P4_P1111](figures/fix1/targseq_all/Zv0370_P3_P4_P1111.png)
![Zv0410_P2_P3_P2111](figures/fix1/targseq_all/Zv0410_P2_P3_P2111.png)
![Zv0420_P3_P1_P1111](figures/fix1/targseq_all/Zv0420_P3_P1_P1111.png)
![Zv0420_P3_P3_P1211](figures/fix1/targseq_all/Zv0420_P3_P3_P1211.png)
![Zv0440_P4_P4_P1111](figures/fix1/targseq_all/Zv0440_P4_P4_P1111.png)
![Zv0470_P2_P5_P1211](figures/fix1/targseq_all/Zv0470_P2_P5_P1211.png)
![Zv0480_P3_P1_P1211](figures/fix1/targseq_all/Zv0480_P3_P1_P1211.png)
![Zv0480_P3_P2_P1211](figures/fix1/targseq_all/Zv0480_P3_P2_P1211.png)
![Zx0030_P2_P1_P2111](figures/fix1/targseq_all/Zx0030_P2_P1_P2111.png)
![Zx0030_P2_P2_P1211](figures/fix1/targseq_all/Zx0030_P2_P2_P1211.png)
![Zx0030_P2_P2_P2111](figures/fix1/targseq_all/Zx0030_P2_P2_P2111.png)
![Zx0090_P2_P5_P1111](figures/fix1/targseq_all/Zx0090_P2_P5_P1111.png)
![Zx0140_P2_P4_P1211](figures/fix1/targseq_all/Zx0140_P2_P4_P1211.png)
![Zx0160_P2_P4_P1111](figures/fix1/targseq_all/Zx0160_P2_P4_P1111.png)
![Zx0160_P2_P5_P1211](figures/fix1/targseq_all/Zx0160_P2_P5_P1211.png)
![Zx0190_P2_P4_P1111](figures/fix1/targseq_all/Zx0190_P2_P4_P1111.png)
![Zx0190_P2_P5_P1211](figures/fix1/targseq_all/Zx0190_P2_P5_P1211.png)
![Zx0220_P2_P4_P1111](figures/fix1/targseq_all/Zx0220_P2_P4_P1111.png)
![Zx0340_P2_P4_P2111](figures/fix1/targseq_all/Zx0340_P2_P4_P2111.png)
![Zx0350_P1_P4_P2211](figures/fix1/targseq_all/Zx0350_P1_P4_P2211.png)
![Zx0360_P2_P2_P1311](figures/fix1/targseq_all/Zx0360_P2_P2_P1311.png)
![Zx0360_P2_P5_P1211](figures/fix1/targseq_all/Zx0360_P2_P5_P1211.png)
![Zx0370_P2_P4_P1311](figures/fix1/targseq_all/Zx0370_P2_P4_P1311.png)
![Zx0380_P2_P4_P1111](figures/fix1/targseq_all/Zx0380_P2_P4_P1111.png)
![Zx0390_P2_P2_P1111](figures/fix1/targseq_all/Zx0390_P2_P2_P1111.png)
![Zx0390_P2_P3_P2211](figures/fix1/targseq_all/Zx0390_P2_P3_P2211.png)
![Zx0400_P1_P4_P2211](figures/fix1/targseq_all/Zx0400_P1_P4_P2211.png)
![Zx0430_P1_P5_P1211](figures/fix1/targseq_all/Zx0430_P1_P5_P1211.png)
![Zx0440_P1_P5_P1211](figures/fix1/targseq_all/Zx0440_P1_P5_P1211.png)
![Zx0440_P2_P5_P1211](figures/fix1/targseq_all/Zx0440_P2_P5_P1211.png)
![Zx0450_P1_P5_P2111](figures/fix1/targseq_all/Zx0450_P1_P5_P2111.png)
![Zx0470_P2_P3_P1211](figures/fix1/targseq_all/Zx0470_P2_P3_P1211.png)
![Zx0490_P4_P2_P1211](figures/fix1/targseq_all/Zx0490_P4_P2_P1211.png)
![Zx0490_P4_P3_P1211](figures/fix1/targseq_all/Zx0490_P4_P3_P1211.png)
![Zx0500_P2_P2_P1111](figures/fix1/targseq_all/Zx0500_P2_P2_P1111.png)
![Zx0500_P2_P3_P1111](figures/fix1/targseq_all/Zx0500_P2_P3_P1111.png)
![Zx0500_P2_P5_P1211](figures/fix1/targseq_all/Zx0500_P2_P5_P1211.png)
![Zx0510_P1_P3_P1211](figures/fix1/targseq_all/Zx0510_P1_P3_P1211.png)
![Zx0510_P1_P3_P2111](figures/fix1/targseq_all/Zx0510_P1_P3_P2111.png)
![Zx0520_P2_P4_P1111](figures/fix1/targseq_all/Zx0520_P2_P4_P1111.png)
![Zx0540_P3_P5_P1111](figures/fix1/targseq_all/Zx0540_P3_P5_P1111.png)
![Zx0540_P3_P5_P3211](figures/fix1/targseq_all/Zx0540_P3_P5_P3211.png)
![Zx0540_P3_P5_P6111](figures/fix1/targseq_all/Zx0540_P3_P5_P6111.png)
![Zx0550_P4_P2_P2511](figures/fix1/targseq_all/Zx0550_P4_P2_P2511.png)
![Zx0550_P4_P2_P4211](figures/fix1/targseq_all/Zx0550_P4_P2_P4211.png)
![Zx0550_P4_P3_P5311](figures/fix1/targseq_all/Zx0550_P4_P3_P5311.png)
![Zx0560_P2_P3_P2111](figures/fix1/targseq_all/Zx0560_P2_P3_P2111.png)
![Zx0560_P4_P3_P5111](figures/fix1/targseq_all/Zx0560_P4_P3_P5111.png)
![Zx0570_P2_P2_P1311](figures/fix1/targseq_all/Zx0570_P2_P2_P1311.png)
![Zx0570_P2_P3_P1311](figures/fix1/targseq_all/Zx0570_P2_P3_P1311.png)
![Zx0580_P2_P5_P2411](figures/fix1/targseq_all/Zx0580_P2_P5_P2411.png)
![Zx0580_P2_P5_P5111](figures/fix1/targseq_all/Zx0580_P2_P5_P5111.png)
![Zx0580_P3_P5_P1211](figures/fix1/targseq_all/Zx0580_P3_P5_P1211.png)

</details>

**How to read the SNPs:** they come in clusters (one targeted region each, often 5–30 SNPs), so points stack at one position; and SNPs of one region can differ (for example 0.3 and 0.95 at chr4 186.9 Mb), probably from amplification bias. So a region is read from all its SNPs together, and a tract from its median (as in the table).

**One plant or a bulk?** This matters only for the SEG row, and it waits for confirmation of what DNA was sequenced:
- **If it is one plant of the plot,** a SEG plot's plant should be AA : AB : BB = 25 : 50 : 25 (child of an AB mother). We see 9 : 83 : 9 among the clearly called tracts (21 : 71 : 8 in the first model), which would mean both homozygotes are strongly depleted. An AB plant's donor allele fraction should also centre on 0.5; we see 0.45 (middle half of the well-covered tracts 0.36–0.52).
- **If it is a bulk of several plants,** the SEG fraction should be about 0.46, which is what we see, and nearly all SEG tracts would look intermediate (AB), which is also what we see.
- Single samples point the same way: in good example 4 (`Zv0480_P3_P1_P1211`, chr5) the fraction sits between AB and BB at amplicon after amplicon, which one plant cannot do.
- The data lean towards a bulk, but the B73 and DONOR rows hold either way: they are the validation.

**Limits.** A tract's call is the median over its SNPs, so one amplicon with many SNPs can outvote the rest (bad example 4); counting per amplicon would be fairer and is not done yet. Only tracts of 1 Mb or more were compared, so short tracts are still unchecked; many tracts had too few panel SNPs for a clear call (31–55% unclear); and two samples were not aligned (one had only an R2 file, one job ran out of time) and were not rerun.

Scripts: `final_script/34_targseq_align.sh`, `34b_targseq_call.sh`, `34c_targseq_compare.py` (site-level comparison), `34d_targseq_tracts.py` (tracts), `34e_targseq_figure.py` (summary figure), `34f_targseq_examples.py` (example panels: `--set good`, `--set bad`; `--gallery_dir` for all samples). Results: `data/introgression_fix1/targseq/` (`tracts_targseq.tsv`, `hidden_b73_windows.tsv`, `per_sample.tsv`, `sites_used.tsv.gz`); the comparison of the models: `data/introgression_fix1/model/compare_old_new.txt`.

### 2.10 Population structure from k-mers (step 32)

The HMM calls say *whether* a window carries teosinte, but not *which* teosinte, because they are built only from B73 k-mers that are present or missing. A PCA of the ancestry calls therefore does not separate the five taxa (taxon explains ≤ 2.5% of the variance of each of the first 10 PCs; the PCs follow tracts shared within families).

To see which teosinte each pool carries, we used the k-mers that the donor DNA brings in and B73 does not have, which differ between teosinte taxa. Following the idea of reading population structure from k-mer content without alignment (Hrytsenko et al. 2025), step 32 uses **sourmash** (Brown & Irber 2016; Irber et al. 2022) to compare every pool with the teosinte assemblies:
1. **Sketch** (32a): k = 31, FracMinHash with scaled = 1,000 (every k-mer whose hash is in the lowest 1/1,000 of the hash space; the same k-mers are kept in every sample, so sketches compare directly). One sketch for B73, one per assembly (Zd Gigi, Zd Momo, Zh RIMHU001, Zl RIL003, Zv TIL01, Zv TIL11, Zx TIL18, Zx TIL25), one per pool from its trimmed reads.
2. **Taxon-specific k-mers** (32b): for each assembly, its k-mers minus B73's; for each taxon, the union over its assemblies minus the k-mers of the other four taxa.
3. **Profile per pool** (32b): the share of each taxon's specific k-mers found in the pool's reads, f_taxon. Its largest value is the pool's **k-mer taxon**, which is checked against the pedigree taxon. The share of all teosinte k-mers found, f_teo, measures how much donor DNA the pool has, independently of the HMM.
4. **PCA** (32c): each pool's f values divided by their sum (depth and donor amount cancel; only *which* teosinte is left), centred log-ratio, `prcomp`. Points coloured by pedigree taxon; arrows = taxa.

The results below are for the 1,385 teosinte donor pools; the 29 landrace pools are left out.

![PCA of the teosinte k-mer profiles](data/introgression_fix1/sourmash/pca/kmer_pca.png)

- **The k-mer content separates the taxa.** PC1 and PC2 (79% of the variance) pull *diploperennis* and *luxurians* apart from the rest; PC3 separates *huehuetenangensis*; PC4 separates *parviglumis* from *mexicana*, which overlap most, as expected for the two closest subspecies.
- **The k-mer taxon matches the pedigree taxon in 1,300 of 1,385 pools (93.9%).**

  | Pedigree taxon | Pools | Same k-mer taxon | Other k-mer taxon |
  |---|---|---|---|
  | Zd *diploperennis* | 249 | 216 (86.7%) | Zx 19, Zv 14 |
  | Zh *huehuetenangensis* | 61 | 61 (100%) | |
  | Zl *luxurians* | 118 | 114 (96.6%) | Zv 2, Zx 2 |
  | Zv *parviglumis* | 380 | 340 (89.5%) | Zx 40 |
  | Zx *mexicana* | 577 | 569 (98.6%) | Zv 7, Zd 1 |

- **The mismatches are mostly weak assignments.** `ratio` = the pool's largest f_taxon ÷ its second largest. Pools that disagree have a median ratio of 1.08, i.e. two taxa are nearly tied, against 1.70 for pools that agree. Among pools with a clear assignment (ratio ≥ 1.5), 859 of 862 agree (99.7%). Pools that disagree also carry less donor DNA (HMM donor share, median 5.4% vs 7.8%).
- **The *parviglumis* pools that look like *mexicana* come from high-elevation accessions.** Of the 40, 34 are from accessions above 1,400 m, and the share rises with elevation (Spearman ρ = 0.60 across the 36 Zv accessions, p = 10⁻⁴):

  | Zv accessions at | Pools with k-mer taxon Zx |
  |---|---|
  | < 1,000 m (12 accessions) | 0 of 120 |
  | 1,000–1,400 m (13) | 6 of 130 (5%) |
  | 1,400–1,600 m (7) | 10 of 86 (12%) |
  | 1,600–2,000 m (4) | 24 of 44 (55%) |

  > **Worth noting: three *parviglumis* accessions look like *mexicana* in most of their pools.**
  > - Zv.0200 (Cerro Nananchi, Morelos, 1,400 m): 6 of 6 pools
  > - Zv.0250 (Ruinas de Malinalco, México, 1,901 m): 10 of 11 pools
  > - Zv.0270 (Amatlán de Quetzalcóatl, Morelos, 1,654 m): 9 of 12 pools
  >
  > This could reflect the known gene flow between the two subspecies (Aguirre-Liguori et al. 2019), or how well the two Zv and two Zx assemblies represent each taxon; it is a hypothesis, not evidence of *mexicana* introgression into these accessions: higher-elevation *parviglumis* is closer to the elevation range of *mexicana*. They are also close to *mexicana* on the map: 23, 25 and 45 km from the nearest Chalco collection site, among the five closest of the 36 *parviglumis* accessions (median 97 km). A labelling error would more likely affect single pools than most pools of an accession. These three are worth checking against the seed records. In the k-mer HMM, these pools also fit a lower DONOR level than the Zv assemblies give (scale factor a: median 0.8, against 0.9 for the other Zv pools), as expected if their donor DNA is partly *mexicana*-like (Zx keeps fewer B73 k-mers than Zv: d = 0.37 vs 0.43).
  >
  > One other explanation should be ruled out first. The taxon-specific k-mers come from only two assemblies per taxon. If the two *parviglumis* assemblies (TIL01, TIL11) were collected far from these accessions, these accessions would share fewer of their specific k-mers and look relatively more *mexicana*-like even without gene flow. Where TIL01, TIL11, TIL18 and TIL25 were collected should be checked.
- **In *diploperennis*,** the mismatches come mostly from Zd.0030 and Zd.0040 (both from Milpillas, 1,390 m, listed as race "Zea-diploperennis-Na"): their assignments are weaker (median ratio 1.8 and 1.5, against 3.0 and 3.5 for Zd.0010 and Zd.0020), so they are less like the two Zd assemblies.

![Per accession: pools whose k-mer taxon differs from the pedigree taxon, by elevation](data/introgression_fix1/sourmash/pca/kmer_mismatch_elevation.png)

Next, we asked whether this *mexicana*-like signal follows elevation across all *parviglumis* accessions, not only in the pools that were assigned to the wrong taxon. For this, each pool got a ***mexicana* index**: of the taxon-specific k-mers it shares with either *mexicana* or *parviglumis*, the fraction that are *mexicana*, f_Zx / (f_Zx + f_Zv). An index near 0 means the pool looks like *parviglumis*, near 1 like *mexicana*. We took the median over the pools of each accession and plotted it against the accession's elevation (right panel below; the left and middle panels are the PCA coloured by elevation).
- **Within *parviglumis*, the index rises with elevation** (Spearman ρ = 0.49, 36 accessions, p = 0.003), and the highest accessions reach the values of *mexicana*. So the two subspecies form a continuous gradient, not two separate groups.
- **The index also rises from west to east** (with longitude, ρ = 0.73). Elevation still adds to it once longitude is accounted for (partial ρ = 0.45, p = 0.006), but with 36 accessions, geography and elevation cannot be fully separated.
- **Within *mexicana*, the index does not change with elevation** (ρ = 0.12, 37 accessions).
- **How much teosinte DNA a pool carries does not depend on its donor's elevation** (HMM donor share, ρ = 0.00), so the gradient is about *which* teosinte, not how much.

![k-mer PCA coloured by elevation, and the mexicana index along elevation](data/introgression_fix1/sourmash/pca/kmer_pca_elevation.png)

We then looked within *mexicana*, at its four races: Chalco, Mesa Central, Durango and Nobogame (section 2.1). The taxon-level profile cannot separate them, since all four are *mexicana*. Instead, we used each pool's share of the k-mers of each of the 8 assemblies separately, and compared the two *mexicana* assemblies, TIL18 and TIL25 (figure below; collection sites from `data/metadata_BZea.csv`, where 75 of the 126 teosinte longitudes lack the minus sign; the scripts correct this, since all sites are in Mexico and Guatemala):
- **A, collection sites:** Chalco comes from the Valley of Mexico in the southeast of the range, Mesa Central from the central plateau to its west, and Durango and Nobogame from the north.
- **B, PCA of the 577 *mexicana* pools:** Chalco separates from the other races along PC2, and Mesa Central lies between Chalco and the northern races.
- **C, TIL18 vs TIL25:** Chalco pools share relatively more k-mers with TIL18 (median log ratio +0.20), while Mesa Central (−0.07), Durango (−0.18) and Nobogame (−0.16) are closer to TIL25.
- **D, *mexicana* index along elevation:** the *mexicana* races and the *parviglumis* accessions on one axis.

So Chalco differs from the other three races in its k-mer content, as expected from geography. It is also the race that carries the late-flowering allele at chr5 221.8 Mb (section 3.6.14). Which collections TIL18 and TIL25 come from should be checked before reading more into the direction of panel C.

![The four mexicana races](data/introgression_spread/sourmash/pca/mexicana_races.png)

Finally, the k-mer profile gives a second measure of the HMM's donor share. The HMM measures the donor share from *B73* k-mers that are missing; sourmash measures it from *teosinte* k-mers that are present (f_teo). They use different k-mers from the same reads, so the check is complementary rather than fully independent. After correcting f_teo for sequencing depth (residual of log f_teo on log sketch size, within taxon), the two agree closely across pools: Spearman ρ = 0.91 (Zd), 0.97 (Zh), 0.97 (Zl), 0.93 (Zv) and 0.98 (Zx) (1,266 pools, NOISY pools left out; the first model gave 0.91–0.98).

![Teosinte k-mers vs the HMM donor share](data/introgression_fix1/sourmash/pca/kmer_vs_hmm.png)

Per-pool shares of each taxon's k-mers: `data/introgression_fix1/sourmash/pca/kmer_taxon_shares.png`. Tables: `kmer_taxon_profile.tsv` (per pool), `pca/kmer_taxon_agreement.tsv`, `pca/kmer_accession_agreement.tsv` (per accession, with elevation), `pca/kmer_vs_hmm.tsv`, `pca/kmer_pca_scores.tsv`.

```bash
N=$(awk -F'\t' 'NR>1 && ($2=="control"||$2=="donor")' ${D11}/pools.tsv | wc -l)
J=$(sbatch --parsable --array=1-$(( 9 + (N + 9) / 10 ))%60 32a_sourmash_sketch.sh)   # needs: pip install sourmash
sbatch --dependency=afterok:$J 32b_kmer_taxon_profile.sh
```
Output in `${D32}` (copied to `data/introgression_fix1/sourmash/`): `kmer_taxon_profile.tsv` (per pool: f per taxon, f_teo, containment per assembly, k-mer taxon, agreement) and `pca/`. The PCA step alone, on a laptop: `POOL_SUMMARY=<pool_summary.tsv> Rscript final_script/32c_kmer_pca.R kmer_taxon_profile.tsv pca data/metadata_BZea.csv` (needs R with data.table, ggplot2, patchwork).

### 2.11 Noisy pools and resequencing

119 of the 1,385 donor pools are flagged **NOISY** (σ > 0.15; 93 in the first model). They are left out of the results above. Almost all of them are on three plates with low sequencing depth:

| Plate | Pools | NOISY | Median σ | Mean single-copy k-mers per pool | NOISY by taxon |
|---|---|---|---|---|---|
| PN2 | 3 | 0 | 0.123 | 189 M | — (rest of PN2 is Zm, skipped) |
| PN3 | 85 | 3 | 0.092 | 192 M | Zx 3 |
| PN4 | 93 | 0 | 0.084 | 211 M | — |
| PN5 | 92 | 0 | 0.100 | 170 M | — |
| PN6 | 93 | 2 | 0.112 | 145 M | Zx 1, Zd 1 |
| PN7 | 91 | 0 | 0.083 | 213 M | — |
| PN8 | 92 | 0 | 0.090 | 199 M | — |
| PN9 | 94 | 2 | 0.111 | 162 M | Zl 2 |
| PN10 | 93 | 1 | 0.086 | 245 M | Zv 1 |
| **PN11** | **93** | **36** | **0.133** | **150 M** | **Zl 18, Zv 18** |
| PN12 | 95 | 5 | 0.098 | 211 M | Zv 5 |
| **PN13** | **88** | **32** | **0.143** | **115 M** | **Zd 18, Zv 14** |
| PN14 | 94 | 1 | 0.091 | 207 M | Zd 1 |
| **PN15** | **93** | **36** | **0.145** | **109 M** | **Zv 34, Zd 2** |
| PN16 | 95 | 0 | 0.088 | 214 M | — |
| PN17 | 96 | 1 | 0.085 | 217 M | Zx 1 (PN17_SID1585; see below) |

**Why 26 more than in the first model.** In the windows that switch to informative k-mers (2.6.6) there are fewer k-mers, so a pool's σ rises a little. 25 *parviglumis* pools that sat just under the cut-off (σ 0.140–0.150) moved just over it (0.150–0.160). One *mexicana* pool, PN17_SID1585, rose from 0.077 to 0.226, which points to a problem with its recount rather than its reads; it is to be checked.

**Recommendation: resequence PN11, PN13 and PN15.**
- **Depth:** PN13 and PN15 have about half the k-mers of a typical plate (109–115 M vs ~200 M), and PN11 about 75%. On PN11, 39% of pools are NOISY, on PN15 39% and on PN13 36%.
- **Taxa:** the noisy pools are concentrated in a few taxa, so leaving them out takes a large bite out of those taxa:
  - Zl: 18 of its 20 noisy pools are on PN11.
  - Zd: 18 of 22 are on PN13.
  - Zv: 66 of 72 are on PN11, PN13 and PN15; *parviglumis* loses 72 of its 380 pools.

**Other samples to check.** Besides depth, a few samples look wrong for other reasons. They are kept in the results but flagged; together they change the main numbers by less than 0.1 percentage point.
- **Donor pools with almost no donor.** 5 pools are as clean as the B73 controls (donor share ≤ 0.11%): PN4_SID334 and PN4_SID335 (both from BC2 plant `Zx.0020_P3_P1_P2`), PN14_SID1313 and PN14_SID1314 (both from `Zd.0010_P2_P5_P4`), and PN7_SID649. Losing every donor segment by chance is very unlikely, and the two pairs each come from one BC2 plant. So these BC2 plants may have been B73 (an outcross or a seed mix-up), or the samples may be mislabelled B73 or check plants. Worth checking against the field and seed records.
- **Sample sheet.** PN13_SID1225 is listed twice, as `Zd.0020_P3_P1_P2.5.1.1` and as PurpleCheck; it is left out. The sheet also has 21 PurpleCheck and 13 B73 samples. The checks are not donor lines and are left out. 10 of the B73 samples are the controls; the other 3 (PN1_SID85, PN3_SID236, PN5_SID468) are not in the k-mer run and should be checked.

The full statistics and the per-plate table are in [`final_script/README.md`](final_script/README.md).

### 2.12 Interactive figures (step 23)

[`data/introgression_fix1/introgression_viewer.html`](data/introgression_fix1/introgression_viewer.html) is one self-contained page with every pool's chromosome plot from the current model ([download it](https://raw.githubusercontent.com/nirwan1265/BZea_genotyping/main/data/introgression_fix1/introgression_viewer.html) and open it in a browser). The first model's page is `final_script/figures/introgression_viewer_spread.html`.
- **Drop-down menus:** taxon → F1 family → "All pools" or one pool. There is also a search box for a pool ID.
- **Family view:** each chromosome shows one track per pool of the F1 family, so the segments the pools share line up.
- **Pool view:** the 10 chromosomes of one pool.
- **Hover** a segment for its position, length and state. A table under the figure lists the F1 family's pools (σ, flag, donor share, SEG, DONOR, tracts).
- **Download SVG / PNG** saves the figure on screen, legend included, for papers and slides.
- B73 is blue, SEG yellow and DONOR red, the same colors as in `blocks_<GROUP>.pdf`.

To view it, download the file and open it in a browser. GitHub shows HTML files as code, not as a page.

### 2.13 Limitations and future work

**Short tracts aren't proven yet.** The k-mer HMM calls tracts down to about 400 kb. Short tracts are expected from the crossing scheme ([Tract lengths](#265-how-long-are-the-tracts)), and they are shared by relatives as often as long tracts (below), but single short calls haven't been checked with independent data.

**SEG vs DONOR is checked with targeted sequencing** ([2.9.1](#291-independent-check-targeted-sequencing-of-83-plots)): our DONOR tracts are homozygous teosinte (89%; 95% without the three bad samples), our SEG tracts heterozygous (83% AB, donor allele fraction 0.45) and our B73 tracts clean (96%). What is still open is whether the targeted-sequencing DNA is one plant or a bulk, which decides how the SEG row is read, and tracts shorter than 1 Mb, which were not compared. The check found one systematic error of the first model, SEG tracts between two DONOR tracts that are really DONOR, recurring at four places; the current model (2.6.6) removes most of it, but 57% of the remaining SEG tracts between two DONOR tracts are still BB, and *parviglumis* chr3 75–78 Mb is not fixed. The current model also has slightly more short tracts (1.6 instead of 1.3 tracts under 1 Mb per pool), which are not checked.

**The DONOR level still starts from the assemblies.** The model needs to know how much of each window's B73 k-mers a teosinte genome keeps (d_w). The best source would have been **the actual donor parents** of our lines, sequenced. We don't have them, so d_w starts from teosinte assemblies from MaizeGDB (PanAnd), other plants of the same taxa, and is then learned from our pools (2.6.6). Limits that remain:
- One level per taxon and window: within a taxon, families carry different teosinte haplotypes (the per-accession levels tried in `new_model/` did not do better, and *parviglumis* chr3 75–78 Mb may be such a case).
- The informative k-mers are defined by the assemblies (missing from all of them); k-mers that our donors lack but an assembly carries are not used.
- In switched windows, the DONOR level is learned from pools the current model's first steps called DONOR, so it is partly circular; the targeted sequencing, which is independent, says the result is better.
- Zv and Zx still switch between SEG and DONOR about twice as often as the other taxa.
- Zm (landraces) can't be done yet, because it needs a maize-based donor level (e.g. from NAM founder assemblies).

A first check of whether the short tracts are real uses **relatives**. A real introgression is inherited: pools that descend from the same plant received their donor segments from that plant, so a real tract in one pool should often be found at the same place in its relatives. A false tract (noise) appears at random places, so it should be found in relatives no more often than in unrelated pools.

The relatives are read from the pedigree name. For example, in `Zv.0200_P3_P2_P1.1.1.1`:
- `Zv.0200` is the accession (the teosinte donor);
- `P3` is the F1 plant, `P2` the BC1 plant and `P1` the BC2 plant;
- `.1.1.1` are the three rounds of selfing that led to this plot.

So `Zv.0200_P3_P2_P1.1.1.1`, `Zv.0200_P3_P2_P1.2.1.1` and `Zv.0200_P3_P2_P1.3.1.1` all descend from the same BC2 plant, `Zv.0200_P3_P2_P1`, through different selfed lines.

The check works like this (`final_script/24c_relatives_check.py`; current model, 1,266 pools, 19,144 tracts; SEG and DONOR counted together as "introgression"):
1. Take one tract in one pool.
2. Take a group of other pools with a given relationship to it, e.g. all pools from the same BC2 plant.
3. Count how many of those pools also have an introgression somewhere inside the tract.
4. Repeat for every tract in every pool, separately for short (< 10 Mb) and long (≥ 10 Mb) tracts.

The table gives the result as a percentage: for example, 31.8% means that for a short tract in one pool, on average 31.8% of the other pools from the same BC2 plant have an introgression in the same place.

| The other pools | Short tracts (< 10 Mb) | Long tracts (≥ 10 Mb) | First model: short | First model: long |
|---|---|---|---|---|
| Same taxon, other accession (unrelated: the baseline) | 15.7% | 20.7% | 16.0% | 20.8% |
| Same accession, other F1 plant | 17.0% | 21.0% | 17.0% | 21.3% |
| Same F1 plant, other BC1 plant | 16.3% | 20.8% | 16.5% | 21.0% |
| Same BC1 plant, other BC2 plant | 17.4% | 23.8% | 18.1% | 24.2% |
| **Same BC2 plant** | **31.8%** | **40.6%** | **31.9%** | **40.4%** |

(Tables: `data/introgression_fix1/checks/relatives_check.tsv`, `data/introgression_spread/checks/relatives_check.tsv`. The numbers in `old_model_README.md` (31.4% / 43.8% for the same BC2 plant) came from an earlier version of this count that was not kept; the first-model columns here are recounted with step 24c, so the two models are compared with the same rule.)

- **The baseline is not zero.** Introgressions are common and some genome regions carry them in many pools, so even unrelated pools share some by chance. That's why every row is compared with the baseline, not with zero.
- **Pools from the same BC2 plant share tracts about twice as often as unrelated pools** (31.8% vs 15.7% for short tracts, 40.6% vs 20.7% for long tracts).
- **Short tracts rise by the same factor as long tracts.** If the short tracts were mostly noise, they would not be inherited, and their sharing would stay at the baseline even among close relatives. Since they rise like the long tracts, the short tracts are inherited like the long ones. This is pedigree-consistency evidence, not proof that each call is right: it doesn't tell how many calls are false, or how genome position and donor–assembly similarity affect sharing.
- **The current model's tracts are inherited as well as the first model's.** Every row is within 0.7 points of the first model, although the current model has slightly more short tracts (1.6 instead of 1.3 under 1 Mb per pool; 2.6.5).
- **Pools from the same accession but a different BC2 plant share hardly more than the baseline** (17.0% vs 15.7%, 21.0% vs 20.7%). If the short tracts were an artifact of the accession's donor plant being different from the teosinte assemblies, they would appear at the same places in every pool of that accession. They don't, so that artifact is not what makes them.
- **Open question:** pools from the same BC1 plant but different BC2 plants should also share more than the baseline, about half the increase seen for the same BC2 plant, because each BC2 plant inherits about half of its BC1 parent's donor segments. They share only a little more (23.8% vs 20.7% for long tracts: about 15% of the increase seen for the same BC2 plant, not half). Why is not known yet.

**Future work:**
1. **Read check, with data we already have.** The SNP route (steps 01–12) counted each pool's B73 and teosinte reads at SNP sites. Those reads come from mapping, which loses many donor reads, but they still show which alleles are present:
   - **B73 regions:** almost only B73 reads.
   - **DONOR regions:** mostly teosinte reads.
   - **SEG regions:** a mix of both, because the pool holds both B73 and teosinte chromosomes.
   - **Short tracts:** if real, they should show teosinte reads like longer tracts of the same state; if noise, they look like the B73 around them.
2. **BC1 ancestors.** Every real tract in a pool must lie inside a teosinte region of the BC1 plant it descends from. This needs the link between each pool and its sequenced BC1 plant.
3. **Gold standard: KASP markers on single plants.** A few markers inside short tracts and inside SEG tracts, genotyped on single plants from the same plots:
   - **Short tract:** real if plants carry the teosinte allele there.
   - **SEG tract:** plants should be AA, AB and BB in about 1:2:1.
   - **DONOR tract:** all plants BB.
   - **B73 region:** all plants AA.
4. **Donor plants.** If seed or tissue of the donor parents (or of the same accessions) can be found, 10–20x short-read sequencing gives their own d_w directly, and their own informative k-mers. No assembly is needed.
5. **Teosinte-side k-mers.** Counting k-mers present in the teosinte but missing from B73 (placed on B73 by alignment) would add signal on the DONOR side, where the current model only sees B73 k-mers going missing; this could help where contrast is still low (e.g. *parviglumis* chr3).
6. **Noisy plates.** Resequence PN11, PN13 and PN15.

---

## Section 3 — k-mer ancestry GWAS and joint QTL mapping

We ran **two analyses** on the same genotypes (k-mer teosinte ancestry per 100-kb window, Section 2) and the same phenotypes (two-year SpATS BLUEs):

| | **Analysis 1: k-mer ancestry GWAS** (step 26, GridLMM) | **Analysis 2: joint QTL model** (step 29.2, R/qtl2) |
|---|---|---|
| Question | Where in the genome are the QTL? | How many separate QTL, how large is each, and does the effect differ between teosinte taxa? |
| What is tested | Every window, one at a time | The QTL found by Analysis 1, plus any second QTL that a conditional scan finds inside one of its regions; all fitted together |
| Effect | One effect per teosinte copy, shared by all taxa | The effect of each QTL adjusted for the others, and its effect within each taxon |
| Significance | 5% genome-wide threshold from permutations | Bonferroni over the QTL of that trait; the QTL were chosen from the same data, so this characterizes the scan's signals rather than confirming them |
| Output | Manhattan plot, regions, genes in each region | QTL table, per-taxon effects, test that taxa differ |

Both use taxon as a fixed effect and a leave-one-chromosome-out (LOCO) kinship built from the same ancestry calls. At the same windows the two packages give the same answer (DTS: −log10 p 5.12 / 9.68 / 9.08 in GridLMM and 5.09 / 9.67 / 9.08 in R/qtl2), so Analysis 2 describes exactly the QTL of Analysis 1.

### 3.1 Combining a k-mer ancestry GWAS with joint QTL mapping

The marker is not a SNP. It is the **teosinte ancestry** of each plot in each 100-kb window, from the k-mer HMM (Section 2): 0 = B73, 1 = SEG, 2 = DONOR copies per plant. Testing local ancestry against a phenotype is what human genetics calls **admixture mapping** (Smith & O'Brien 2005).

We combine two steps, each doing what it is good at:
- **The k-mer ancestry GWAS finds the peaks.** It is GWAS-style: one mixed-model test per window across the whole genome, with a kinship matrix and a permutation threshold.
- **The joint QTL model describes them.** It fits the peaks together, as in QTL mapping, to tell which peaks are separate QTL, how large each is, and whether the teosinte taxa differ.

Neither analysis can point to a single gene, and the joint model does not narrow the regions. The limit is the length of the donor tracts. A tract is several Mb long (median about 10 Mb, section 2.6.5), so within a plot, all the windows inside one tract have the same ancestry. When the scan tests those windows, it sees almost the same genotypes in each, so they all get similar p-values, and the peak is a broad hill rather than a sharp spike. The scan reports the hill as a region with a support interval (the windows within 1.5 −log10 p of the top window), which here spans about 0.3 to 15 Mb. What the joint model adds is not a narrower region, but the answer to whether two hills are two separate QTL or one QTL seen twice through linked tracts. The genes inside a support interval are therefore candidates, not the identified cause.

### 3.2 Phenotypes

Field trials at Clayton, NC in **2023** (CLY23_D4: 3 reps, 4,432 plots, 1,329 lines) and **2025** (CLY25_B5: 2 reps, 2,898 plots, 1,410 lines); 1,316 lines in both. Each BC2S3 plot is grown from one BC2S2 mother.

BLUEs are two-stage (`final_script/28_spats_blues.R`; details in `final_script/README.md`, section 11):
1. **Per year and trait, SpATS** (Rodríguez-Álvarez et al. 2018): `y = Rep + Border + line (fixed) + 2-D P-spline surface(column, row) + random row + random column + e`. Border = plot next to a border (`B`) cell in the field map. Plots noted as no germination, hybrid, wrong line or sick are removed, plus one pass of outlier removal (|residual| > 3.5 robust SD). Heritability (H², Oakey et al. 2006) from the same model with line random.
2. **Across years**: the year effect from a weighted `lm(BLUE ~ line + year)` (weights 1/SE²), then the weighted mean of the year-adjusted BLUEs.

| Trait | H² 2023 | H² 2025 | r between years |
|---|---|---|---|
| DTS (days to silking) | 0.88 | 0.84 | 0.69 |
| DTA (days to anthesis) | 0.89 | 0.82 | 0.68 |
| ASI (DTS − DTA) | 0.61 | 0.57 | 0.43 |
| PH (plant height) | 0.86 | 0.86 | 0.74 |
| EH (ear height) | 0.78 | 0.74 | 0.63 |
| EN (ear number: nodes with ears) | 0.68 | 0.42 | 0.33 |
| Prolif (prolificacy: most ears on one node) | 0.61 | 0.55 | 0.21 |
| NBR (nodes with brace roots) | 0.52 | 0.47 | 0.34 |
| SPAD (leaf greenness, SPAD meter, 36 days after silking); 2023 only | 0.40 | | |
| LAE (leaves above the ear); 2025 only | | 0.40 | |
| Kinki (how much the stem bends, score 0–2); 2023 only | 0.29 | | |
| StPu (stem pubescence), StPi (stem pigmentation); present/absent, line means | | | 0.16, 0.09 |
| BW (blade width), BL (blade length), SL (sheath length); 2023 only | 0.51, 0.56, 0.73 | | |

### 3.3 Analysis 1: k-mer ancestry GWAS

#### 3.3.1 Genotypes and kinship

- **Plots:** donor pools with a phenotype; low-depth (NOISY) pools are dropped.
- **Windows:** 21,309 of 21,330 100-kb windows with at least 20 plots carrying donor.
- **Additive predictor:** dosage 0 / 1 / 2 (B73 / SEG / DONOR). A SEG plot holds ¼ AA, ½ AB, ¼ BB plants, so its average is 1 teosinte copy per plant; the phenotype is a plot mean, so this is the right scale.
- **Dominance predictor:** SEG = 1 if the plot is SEG, else 0. The additive predictor assumes that a heterozygous plant (AB) is exactly halfway between AA and BB. Dominance means it is not: the teosinte allele can act more strongly (dominant) or more weakly (recessive) when there is only one copy. The SEG plots are the only plots with heterozygous plants, so they carry this information:
  - On the plot-mean scale, with a = the additive effect per teosinte copy and d = the dominance deviation of an AB plant, a B73 plot averages 0, a DONOR plot 2a, and a SEG plot (¼ AA, ½ AB, ¼ BB plants) ¼·0 + ½(a + d) + ¼·2a = a + d/2.
  - So adding SEG to the model with coefficient c gives c = d/2: c > 0 means the teosinte allele is partly dominant, c < 0 partly recessive, c = 0 purely additive.
  - It is tested only where a window has ≥ 5 SEG and ≥ 5 DONOR plots.
  - **Result:** at none of the 120 QTL of the 11 traits does the dominance term stand out (largest −log10 p = 1.95; Bonferroni over 120 QTL would need 3.38). 15 of the 120 have p < 0.05, where 6 are expected by chance. So the QTL behave additively as far as we can tell.
  - **Why the power is low:** only about 75 plots per QTL window are SEG, and in a SEG plot only half the plants are AB, so the heterozygote is diluted. A small dominance effect would not be detected. The permutations were run for the additive test only, so there is no genome-wide threshold for dominance, and the dominance scan is not used to find QTL.
- **Kinship, leave-one-chromosome-out (LOCO):** for the chromosome being tested, K = Z Z′ / p from the centred dosage of all windows on the **other nine** chromosomes, scaled to mean diagonal 1. In BC2S3 lines on a B73 background, relatedness is mostly shared donor tracts, which the dosage captures directly; no SNP calls or imputation are needed.

#### 3.3.2 Why LOCO, and not one K for the whole genome

A kinship matrix that includes the tested region already "contains" the QTL: lines that share the donor tract at the QTL look related, and the random effect absorbs the QTL's effect. This is **proximal contamination** (Listgarten et al. 2012; Yang et al. 2014). LOCO removes it by leaving the tested chromosome out of K (also the default in R/qtl2, Broman et al. 2019).

In a backcross population this matters a lot, because tracts are long and few loci carry most of the donor relatedness.

The table compares the two choices of K on the same data (current model: DTS, 1,260 plots, taxon fixed; first model in brackets), using **λ_GC** (the genomic inflation factor; Devlin & Roeder 1999). λ_GC is the median test statistic over all windows divided by the median expected if no window had any effect. λ_GC = 1 means the test statistics look as expected under no effect; above 1, they are larger than expected (inflated); below 1, smaller (deflated, the model is too conservative). In most GWAS a λ_GC well above 1 is read as a warning of uncorrected population structure.

| K | λ_GC | Best window chr10 (ZmCCT10 region, 94–98 Mb) | Best window chr5 (~221.7 Mb) |
|---|---|---|---|
| **LOCO** | 4.10 (3.97) | −log10 p **9.1**, +1.0 day per teosinte copy (9.1) | **9.4** (9.7) |
| Whole genome | 0.59 (0.65) | 2.7 (2.9) | 6.0 (6.1) |

With the whole-genome K, the chr10 flowering signal almost disappears and λ_GC falls **below 1**: the model over-corrects. LOCO is the fix.

**λ_GC under LOCO (about 4) needs checking.** λ_GC assumes most markers are null and independent. In this population each QTL is tagged by hundreds of linked windows (tracts of several Mb, and with LOCO the tested chromosome's own background is not in K), and the traits are polygenic, so a high median test statistic is expected even without false positives (Yang et al. 2011). But a high λ_GC can also mean the test is miscalibrated, so this explanation is not proven. Significance is set by **permutation** (next section), and whether that threshold holds when lines are related is still to be checked.

#### 3.3.3 Model

GridLMM (Runcie & Crawford 2019), maximum likelihood, one chromosome at a time with its LOCO K:

`y = taxon (fixed) + b · dosage_w + u + e`,  `u ~ N(0, σ²_g K_LOCO)`,  `e ~ N(0, σ²_e I)`

- b = effect per teosinte copy (trait units, e.g. days); likelihood-ratio test for b (1 df).
- Dominance: the same model plus c · SEG_w; likelihood-ratio test for c.
- GridLMM fits the variance proportion on a grid instead of re-optimising it for every window, which makes ~21,000 tests with a 1,258 × 1,258 K fast (~3 min per trait on 4 cores). The code started from `scripts/Analyze_QTL_Founder2.R`.

We used **taxon** (5 levels: Zd, Zh, Zl, Zv, Zx) as the fixed effect, not the **F1 family**. An F1 family is all the plots that descend from one F1 plant (one teosinte plant × B73); the plots come from about 95 F1 families of 82 accessions. Plots of the same F1 family share many of their tracts (relatives check, Section 2, Limitations). With F1 family (95 levels) as a fixed effect, a tract carried by a whole family is confounded with the family mean and absorbed, and 94 extra parameters raise the threshold. Family relatedness is already in K. For DTS (829 plots, one year): the peaks on chr7, chr8/9 and chr10 were similar either way, but the 5% threshold rose from 3.8 to 5.0 and the chr3 peak (4.4) was lost. F1 family as a fixed effect is kept as a check (`FIXED_COL=F1FAMILY`).

#### 3.3.4 Significance: permutations

For each trait, phenotypes are shuffled **within taxon** (keeping the taxon means), the whole genome scan is rerun, and the largest −log10 p in the genome is kept. The 95th percentile of these 1,000 genome-wide maxima is the threshold. It controls the **family-wise error rate (FWER) at 5% across the whole genome**: if the trait had no QTL, a scan would put any window above it in only 5% of cases (Churchill & Doerge 1994). Because each permutation keeps the real correlation between neighbouring windows, the threshold accounts for it; Bonferroni over 21,309 windows (−log10 p = 5.63) would treat them as independent and be too strict. The thresholds come from **1,000 permutations** per trait. With *mexicana* split by race (section 3.6.14), phenotypes are shuffled within race.

**Caveat: the plots are related.** Permutation assumes the phenotypes can be swapped between plots. Plots from the same F1 family share tracts and phenotypes, so shuffling within taxon also breaks that family structure, and the threshold may not give exactly 5% FWER. A check that keeps the relatedness is planned: simulate phenotypes with no QTL from the fitted model (taxon effects + u ~ N(0, σ²_g K) + e), run the whole scan on each, and count how often any window passes the threshold (it should be about 5%). The strongest QTL (e.g. ZmCCT10, −log10 p 9–14) are far above any plausible threshold; QTL close to the threshold are the ones the check could change.

#### 3.3.5 Regions and genes

Step 27 joins windows above the threshold into regions (gap < 5 Mb), keeps one region per locus, and gives each a support interval: the windows within 1.5 −log10 p of the top window. This is a drop in −log10 p, not in LOD, so it is not the usual 1.5-LOD interval, and how often it contains the causal gene has not been calibrated. Every B73 v5 gene overlapping the interval is listed from the GFF3. Output: `genes/regions_genes_<TRAIT>_ALL.tsv` (one row per gene, with the region's position, top −log10 p and effect). Step 30 (`final_script/30_genes_go.R`) adds each gene's window −log10 p and effect and its GO terms: `data/introgression_fix1/GO/genes_GO_all_<TRAIT>.tsv`.

### 3.4 Analysis 2: joint QTL model with per-taxon effects

The scan in Analysis 1 tests each window on its own. A peak can therefore be a real QTL, or a "shadow" of a stronger QTL elsewhere on the chromosome that it shares tracts with. Analysis 2 takes the QTL from the scan and tests them **together**, to answer three questions for each one: is it still significant when the other QTL are in the model, how large is its effect, and is the effect the same in all teosinte taxa? It uses the package R/qtl2 (Broman et al. 2019), with the same plots, phenotypes, taxon and LOCO kinship as Analysis 1 (`final_script/29_joint_qtl.R`).

**Input QTL.** For each region of the scan, the top window and its support interval (`peaks_<TRAIT>_ALL.tsv`, step 27).

**Genotypes.** At each window, a plot's teosinte allele frequency p = dosage / 2 (B73 0, SEG 0.5, DONOR 1). R/qtl2 takes this as the probability of the teosinte allele (B73 allele: 1 − p). The kinship is LOCO, from the same values (`calc_kinship(probs, "loco")`).

**Step 1: more than one QTL in a region? (conditional scans)** A region of the scan is a run of windows above the threshold. With long tracts, one region can span tens of Mb and could hold two QTL. To check, the whole genome is scanned again with taxon and all QTL found so far as covariates:

`y = taxon + Σ_i γ_i · p_i + β · p_w + u + e` (for every window w; i = the QTL already in the model)

If the best window is still above the genome-wide threshold of Analysis 1, it is a QTL that the others do not explain: it is added to the list, and the scan is repeated, until no window passes. Two rules stop it from adding the wrong things: a window must be significant in the ordinary scan too (so this can split a region, but never create a new one), and it must not be in strong LD with a QTL already in the model (r² ≤ 0.5), since two nearly identical predictors give unstable, opposite effects. This is forward selection by conditional scans (Broman & Sen 2009), which uses the main scan's threshold, as is standard; the conditional scan has its own null distribution, so this threshold is approximate. Close peaks with no independent signal stay one QTL. Added QTL are listed in `qtl2_<TRAIT>_added.tsv` and drawn as open triangles in the figures.

**Step 2: joint test and effect of each QTL.** Each QTL j is then fitted with taxon and **all the other QTL** as covariates (`fit1`):

`y = taxon + Σ_{i≠j} γ_i · p_i + μ + a_j · (2 p_j − 1) + u + e`,  `u ~ N(0, σ²_g K_LOCO)`,  `e ~ N(0, σ²_e I)`

- a_j is the **effect per teosinte copy**: going from B73 (p = 0) to DONOR (p = 1) changes the fitted value by 2a_j, i.e. a_j per copy (set with `contrasts = cbind(mu = c(1, 1), a = c(-1, 1))`), with its standard error.
- The **joint test** compares this model with the same model without QTL j: LOD = log10 of the likelihood ratio, converted to −log10 p with 1 df (χ² = 2 ln(10) · LOD).
- If QTL j's signal disappears once the others are in the model, it was a shadow of another QTL, not a separate one.

**Step 3: effect within each taxon.** The same model as step 2 (without the taxon term), fitted only to the plots of one taxon, gives that taxon's effect per copy. It is estimated only where the taxon has ≥ 5 carrier plots (plots with SEG or DONOR at the QTL).

**Step 4: do the taxa differ?** At the QTL window, the model with a separate effect for each taxon (`scan1(..., intcovar = taxon)`)

`y = taxon + μ + Σ_t a_t · 1[taxon = t] · (2 p − 1) + u + e`

is compared with the model with one shared effect a. The LOD difference is converted to −log10 p with (number of taxa − 1) df. A significant result means the teosinte allele at that QTL has a different effect depending on the taxon it came from.

**Significance.** The joint model does not search the genome: for a trait with k QTL it makes k joint tests and k "taxa differ" tests. Each is significant if −log10 p > −log10(0.05 / k) (Bonferroni over the k QTL of that trait; e.g. 2.08 for 6 QTL). Bonferroni controls the family-wise error rate over a fixed set of k tests. Here the k QTL were picked from the scan of the same phenotypes, so the joint p-values are not an independent confirmation, and the effects of QTL found this way tend to be overestimated (winner's curse). The joint model is best read as **characterizing** the scan's signals: which of them remain once the others are fitted, how large each is, and whether the taxa differ. "Hold" in the tables below means the QTL passes this test; independent data would be needed to confirm a QTL.

Every statistical step is an R/qtl2 function; the script only converts the k-mer calls into R/qtl2's format.

#### 3.4.1 Why the joint model starts from the peaks, and not from all windows at once

A joint model needs a short list of QTL to fit together. It cannot take all 21,309 windows: that is far more predictors than the 1,258 plots, and neighbouring windows are nearly identical (one tract covers dozens of them), so their effects cannot be separated. The genome scan is the step that picks the list: one window per peak (the top window), with its support interval. The joint model then asks, for those few QTL together, which ones still matter once the others are in the model, how large each is, and in which taxa. This two-step design (scan, then fit the selected QTL jointly) is how QTL mapping is usually done; TeoNAM's joint linkage mapping chooses its QTL by stepwise selection instead (Chen et al. 2019). Using the scan's peaks keeps the two analyses on exactly the same QTL.

#### 3.4.2 Why some QTL show fewer taxa

A per-taxon effect needs carriers: plots of that taxon with the teosinte allele at the QTL. Taxa with fewer than 5 carrier plots at a QTL get no estimate. For example, DTS Q2 (chr10) has carriers from Zd (23), Zv (24) and Zx (66) but only 4 from *huehuetenangensis* and 2 from *luxurians*, so those two are not shown. Which taxa carry a tract depends on where the backcrosses happened to keep donor segments.

This is the mixed-model counterpart of TeoNAM's joint linkage mapping (Chen et al. 2019). Effects are estimated per taxon, not per accession: a typical window has ~90 carrier plots across 82 accessions, about one per accession. Details: `final_script/README.md`, section 12.

#### 3.4.3 Fine-mapping inside the support interval: how likely is each window?

The support interval (the windows within 1.5 −log10 p of the peak) says where a QTL could be, but treats every window in it the same. To rank the windows, each QTL's support interval was fine-mapped with **SuSiE** ("sum of single effects"; Wang et al. 2020; R package `susieR`). SuSiE fits all windows of the interval together and returns:
- **PIP** (posterior inclusion probability) of each window: the probability that the window holds the causal variant. With one QTL in the interval, the PIPs add up to about 1.
- **95% credible set:** the fewest windows that hold the QTL with at least 95% probability. SuSiE gives one set per independent effect, so two QTL in one interval give two sets.
- **Top windows** (our summary of the PIPs): the fewest windows that hold half of the QTL's probability.

**Model.** The same as in the scan and the joint model: y = taxon + the trait's other QTL (their top windows) + polygenic background (LOCO kinship) + the windows of the interval + error. The polygenic heritability is estimated by REML without the interval; y, the windows and the covariates are then transformed so that the residuals are independent, the covariates are removed, and SuSiE (up to 10 effects; residual and effect variances estimated by `susieR`) is fitted to the result.

**What limits the resolution.** Windows that lie inside the same tract in every pool have the same dosage in every pool, so no method can split them: they share the probability equally. The resolution is set by where the tracts of different pools end inside the interval, not by the 100-kb window size. Calling SNPs instead of counting k-mers would not change this, because every SNP inside a tract carries the same information as its window.

**Result for the 81 QTL that hold in the joint model** (medians; first model, 78 QTL: 4.65 Mb / 120 genes, 3.65 Mb / 101 genes, 1.2 Mb / 30 genes):

| | Width | Genes |
|---|---|---|
| Support interval | 4.5 Mb | 119 |
| 95% credible set | 3.3 Mb | 95 |
| Top windows (half of the probability) | 0.9 Mb | 27 |

- **The credible sets are only a little narrower than the support intervals.** The intervals were not too wide: for most QTL, neighbouring windows really cannot be told apart in this population.
- **One QTL per interval.** 77 of the 81 QTL have one credible set. Four have none, i.e. their probability is spread over the whole interval: BL Q1 (chr2 9.3 Mb; 9.1-Mb interval), PH Q14 (chr10 6.5 Mb; 10.2 Mb), Prolif Q5 (chr2 34.4 Mb; 29.7 Mb) and SL Q11 (chr9 157.3 Mb; 4.6 Mb).
- **For some QTL the probability falls on a few windows** (held QTL whose top windows hold 10 genes or fewer):

| QTL | Support interval (genes) | Top windows (genes) | Best window PIP | Candidate named in 3.6, and where it falls |
|---|---|---|---|---|
| BW Q10 | chr9 0.7–1.8 Mb (5) | 0.7–0.9 Mb (0) | 0.26 | — |
| EN Q8 | chr4 119.2–121.3 Mb (12) | 119.2–119.5 Mb (1) | 0.19 | — |
| Prolif Q1 | chr1 22.6–23.3 Mb (24) | 23.2–23.3 Mb (2) | 0.57 | — |
| ASI Q9 | chr5 217.4–218.8 Mb (56) | 217.8–218.0 Mb (4) | 0.27 | — |
| ASI Q10 | chr7 165.2–166.1 Mb (26) | 165.5–165.7 Mb (4) | 0.30 | Zm00001eb323880 (*bhlh90*): outside the top windows (PIP 0.08) |
| Prolif Q7 | chr3 11.0–11.7 Mb (11) | 11.2–11.5 Mb (4) | 0.26 | Zm00001eb122680 (CYP734A6): in the top windows (PIP 0.16) |
| PH Q13 | chr9 139.0–141.0 Mb (61) | 139.6–139.9 Mb (5) | 0.21 | — |
| EN Q18 | chr10 148.2–149.3 Mb (60) | 148.5–148.7 Mb (7) | 0.34 | — |
| EN Q12 | chr7 165.1–166.6 Mb (56) | 165.2–165.5 Mb (7) | 0.19 | — |
| ASI Q7 | chr4 65.0–71.3 Mb (91) | 66.9–68.0 Mb (8) | 0.05 | — |
| DTS Q1 | chr5 221.7–222.0 Mb (16) | 221.7–221.8 Mb (8) | 0.64 | Zm00001eb257680 (HDR1-like): in the top windows (PIP 0.64) |
| EN Q1 | chr1 22.6–23.7 Mb (34) | 22.9–23.1 Mb (8) | 0.30 | Zm00001eb007670 (*ZmD14*): outside the top windows (PIP 0.01) |
| BL Q2 | chr2 34.7–36.0 Mb (26) | 35.0–35.3 Mb (9) | 0.20 | — |
| DTA Q2 | chr5 221.7–222.3 Mb (26) | 221.7–221.9 Mb (10) | 0.43 | Zm00001eb257680 (HDR1-like): in the top windows (PIP 0.43) |
| Prolif Q15 | chr9 25.4–32.7 Mb (156) | 31.6–32.5 Mb (10) | 0.07 | — |
| SL Q1 | chr1 67.6–71.1 Mb (74) | 69.6–70.4 Mb (10) | 0.12 | — |

- **How to use it.** A gene in the top windows is better supported by our data than a gene elsewhere in the interval, but PIP ranks windows, it does not identify the gene: most top windows still hold several genes.
- **A low PIP does not make a gene a poor candidate.** It means our data do not single out its window, for several reasons:
  - PIP compares windows only by how well their dosage pattern across the pools fits the phenotype. Within one interval the patterns differ in a few pools only (where tracts end), so a small amount of noise moves probability from one window to another, especially when the effect is small.
  - The tract ends come from the HMM and are placed to within a few windows; a causal window can lose probability to its neighbour because a few tract ends are slightly off.
  - 5% of the probability lies outside the 95% credible set by definition, and when the probability is spread thinly (highest window ≤ 0.03), fine-mapping says nothing about which gene it is.
  - Evidence on what a gene does (mutants, expression, other mapping studies) is independent of PIP. A gene with a known role and a low PIP can still be the better candidate than an unknown gene with a high one.

  So in the write-ups below, PIP is given for each named candidate as one more piece of evidence, next to the biology, not as a verdict. Each gene's PIP is in the gene tables (`data/introgression_fix1/qtl_genes/<TRAIT>_qtl_genes.tsv`: `window_pip`, `in_credible_set`, `in_top_windows`), and the trait tables in 3.6 give the top windows of every QTL.

Scripts: `final_script/33a_qtl_gene_tables.py` (genes per QTL), `33_susie_finemap.R` (per trait, ~15 s on a laptop) and `33b_pip_to_gene_tables.py`; results in `data/introgression_fix1/finemap_susie/` (`susie_<TRAIT>_{qtl,cs,pip,genes}.tsv`).

### 3.5 How to run

The current model's calls are made by `new_model/fix1/` (steps fx1_1 to fx1_6; see its README). Everything downstream reads them when `MODEL=fix1` is set (`final_script/config.sh`); without it, the first model's calls are used. No script activates an environment: jobs use the env that is active when you submit.

```bash
# Cluster, from new_model/fix1/ (login node):
bash fx1_7a_stats.sh     # python env active (numpy, scipy, pandas, matplotlib): shares (22), SEG donor copies (22b), viewer (23), figures (24)
bash fx1_7b_gwas.sh      # GWAS env active (GridLMM, qtl2): inputs, scan, 1,000 permutations, peaks + genes, joint QTL,
                         # for taxon and for mexicana split by race (submit_assoc_all.sh), per year and for leaf N 2025
# per year for every trait (from final_script/, GWAS env):
MODEL=fix1 TRAITS="DTS DTA ASI PH EH EN BW BL SL SPAD SPAD2 Prolif Kinki NBR LAE StPu StPi" bash submit_year_check.sh
# Laptop, from the repo root, after copying the outputs to data/introgression_fix1 (rsync commands in new_model/fix1/README.md):
bash new_model/fix1/fx1_8_laptop.sh
#   genes + GO (30), genes per QTL (33a), SuSiE fine-mapping (33) and PIPs into the gene tables (33b),
#   side-by-side figures (31), DTS by year (31b), landscape (24b), targeted sequencing (34e), k-mer PCA (32c), checks (incl. relatives, 24c)
T="--sites data/introgression_fix1/targseq/sites_used.tsv.gz --tracts data/introgression_fix1/tracts_all.tsv.gz --targseq_tracts data/introgression_fix1/targseq/tracts_targseq.tsv"
python final_script/34f_targseq_examples.py $T --set good --out figures/fix1/targseq_examples.png       # the four good examples
python final_script/34f_targseq_examples.py $T --set bad --out figures/fix1/targseq_examples_bad.png    # the four bad examples
python final_script/34f_targseq_examples.py --sites data/introgression_fix1/targseq/sites_used.tsv.gz \
  --tracts data/introgression_fix1/tracts_all.tsv.gz --gallery_dir figures/fix1/targseq_all \
  --targseq_tracts data/introgression_fix1/targseq/tracts_targseq.tsv          # all targeted-sequencing samples
python final_script/35_trait_tables.py --base data/introgression_fix1 --out <dir>  # the trait tables of 3.6
```
All results are in `data/introgression_fix1/` (figures in `figures/fix1/`):
- `assoc_gridlmm_spats/`: scans, summaries, 1,000 permutations (`perm_<TRAIT>_seed<N>.tsv`), Manhattan plots, `genes/` (peaks and every gene in each region)
- `joint_qtl/noperm/`: joint QTL tables, per-taxon effects, genes per QTL
- `finemap_susie/`: SuSiE fine-mapping per QTL (credible sets, PIP per window, top windows and their genes)
- `qtl_genes/`: every gene in each QTL interval, with GO terms, its window's −log10 p and its window's PIP
- `GO/`: every gene in the regions with GO terms and its window's −log10 p (`genes_GO_all_<TRAIT>.tsv`)
- `assoc_gridlmm_spats_subtaxon/`, `joint_qtl/noperm_subtaxon/`, `GO_subtaxon/`: the same with *mexicana* split by race
- `*_y2023`, `*_y2025`: each year on its own; `*_N2025`: predicted leaf nitrogen (2025)

### 3.6 Results

**Overview (1,098–1,260 plots per trait; thresholds from 1,000 permutations):**

| Trait | Plots | 5% threshold (−log10 p) | QTL | Hold in joint model | Strongest QTL |
|---|---|---|---|---|---|
| DTS | 1258 | 5.50 | 2 | 2 | chr5 221.8 Mb (−log10 p 9.4, +1.00 d per copy) |
| DTA | 1258 | 5.00 | 6 | 3 | chr10 101.2 Mb (−log10 p 9.4, +0.93 d per copy) |
| ASI | 1258 | 4.42 | 10 | 6 | chr7 165.6 Mb (−log10 p 13.1, -0.37 d per copy) |
| NBR | 1248 | 4.04 | 8 | 3 | chr10 98.2 Mb (−log10 p 9.3, +0.18 per copy) |
| SPAD | 1098 | 4.22 | 8 | 3 | chr1 57.9 Mb (−log10 p 6.2, -3.68 per copy) |
| PH | 1257 | 4.17 | 15 | 10 | chr10 99.5 Mb (−log10 p 11.0, +6.18 cm per copy) |
| EH | 1256 | 4.11 | 16 | 9 | chr4 196.6 Mb (−log10 p 7.2, +3.03 cm per copy) |
| EN | 1256 | 4.49 | 18 | 12 | chr2 18.9 Mb (−log10 p 18.3, +0.19 per copy) |
| Prolif | 1260 | 4.72 | 15 | 9 | chr1 23.2 Mb (−log10 p 19.4, +0.21 per copy) |
| SL | 1170 | 4.20 | 13 | 11 | chr7 141.7 Mb (−log10 p 9.2, -0.33 per copy) |
| BL | 1171 | 4.15 | 7 | 6 | chr2 35.1 Mb (−log10 p 10.7, -1.82 per copy) |
| BW | 1171 | 4.21 | 12 | 7 | chr5 41.2 Mb (−log10 p 9.4, -0.32 per copy) |

"Hold in joint model" = QTL whose joint −log10 p passes Bonferroni over that trait's QTL. The other five traits are in [3.6.15](#3615-other-traits).

**Reading the figures.** Left: the k-mer ancestry GWAS; −log10 p per 100-kb window, dotted line = 5% permutation threshold, QTL numbered Q1, Q2 … by position; an open triangle marks a QTL added by the conditional scans (a second QTL inside a significant region, section 3.4). Middle: the same QTL in the joint model; effect per teosinte copy (trait units) with 95% CI, black = joint model across all taxa, colours = within each taxon (taxa with ≥ 5 carrier plots). Right: two columns of **−log10 p** (not effects): "Joint test" = is the QTL still significant with all other QTL in the model; "Taxa differ" = does its effect differ between taxa. For example, DTS Q1 "9.5 | 13.9" means p = 10⁻⁹·⁵ for the QTL and p = 10⁻¹³·⁹ for the difference between taxa. Positive effect = the teosinte allele increases the trait (for DTS/DTA: later flowering). **Fine-mapping PIP** (section 3.4.3): the probability of the window holding each known gene, in the order the genes are listed; **top** = in the top windows (the fewest windows holding half of the QTL's probability), CS = in the 95% credible set but outside the top windows, outside CS = outside the 95% credible set, flat = the probability is spread over the whole interval (highest window ≤ 0.03), so fine-mapping does not separate its genes. A low PIP does not rule a gene out (3.4.3).

#### 3.6.1 Days to silking (DTS)

![DTS: k-mer ancestry GWAS and joint QTL model](figures/fix1/joint_qtl_DTS.png)

Two QTL pass the 1,000-permutation threshold (−log10 p > 5.50), and both hold in the joint model (Bonferroni over 2 QTL: −log10 p > 1.60). They are the same two QTL as in the first model.

| QTL | Known genes | Fine-mapping PIP | Effect (joint) | Joint −log10 p | Taxa differ | What it says |
|---|---|---|---|---|---|---|
| **Q1 chr5 221.8 Mb** | Zm00001eb257680 (HDR1-like) | 0.644 (**top**) | +1.00 d | 9.5 | 13.9 (yes) | Differs between taxa (Zd +1.95, Zl +4.77, Zv +0.41, Zx +0.50). Strong, but mostly *luxurians* and *diploperennis*; Zv and Zx only about +0.5 d. |
| **Q2 chr10 97.8 Mb** | Zm00001eb418700 (ZmCCT10) | 0.019 (**top**) | +1.01 d | 9.0 | 1.0 (no) | No significant difference between taxa (Zd +1.18, Zv +0.78, Zx +0.94). The ZmCCT10 (Zm00001eb418700) region: teosinte alleles delay flowering everywhere. |

**Suggestive peaks** (above −log10 p 4.5 but below the 5.50 threshold, outside the two QTL): chr5 217.2 Mb (5.5), chr1 84.5 Mb (5.3; the DTA Q1 peak), chr10 81.8 Mb (5.3), chr9 10.9 Mb and chr9 1.4 Mb (5.0 each; both DTA QTL, teosinte earlier), chr10 69.2 Mb (4.8), chr10 125.0 Mb (4.6).

**Candidate genes:**
- **Q2, chr10:** *ZmCCT10* (**Zm00001eb418700**, annotated "GHD7_2"; 96.18 Mb), the main teosinte photoperiod gene; teosinte alleles delay flowering under long days (Hung et al. 2012; Yang et al. 2013). *Fine-mapping (3.4.3): PIP 0.019; in the top windows, but these are wide (93.3–97.8 Mb, 43 genes), so fine-mapping does not single it out.*
- **Q1, chr5:** **Zm00001eb257680**, annotated as **Protein HEADING DATE REPRESSOR 1** (SnRK1-interacting protein; GO "regulation of flower development"), in the top window (221.7–221.8 Mb, 8 genes). Its functional ortholog in rice, *HDR1* (*Heading Date Repressor 1*), is a nuclear protein that suppresses flowering under long days: it upregulates *Hd1* and downregulates *Ehd1*, which lowers the florigens *Hd3a* and *RFT1*; it interacts with the kinase OsK4, another long-day repressor (Sun et al. 2016). The *hdr1* mutant flowers early in long days only, i.e. HDR1 acts through the photoperiod pathway, which fits a teosinte allele that delays silking. Not characterised for flowering in maize. Orthology still to be confirmed with a gene tree. The same window is the DTA Q2 peak. *Fine-mapping (3.4.3): PIP 0.644 (the highest in the interval); in the top windows (8 genes).*

<details><summary>All DTS QTL with per-taxon effects and all genes in the interval</summary>

| QTL | Position | Interval (Mb) | Scan −log10 p | Effect per copy (joint) | Joint −log10 p | Taxa differ −log10 p | Per-taxon effects (carrier plots) | Genes in interval | Top windows, Mb (genes) | Best window PIP |
|---|---|---|---|---|---|---|---|---|---|---|
| Q1 | chr5 221.8 Mb | 221.7–222.0 | 9.4 | +1.00 d ± 0.16 | **9.5** | **13.9** | Zd +1.95 (10), Zl +4.77 (13), Zv +0.41 (46), Zx +0.50 (51) | [16 genes](data/introgression_fix1/qtl_genes/DTS_qtl_genes.tsv) | 221.7–221.8 (8) | 0.64 |
| Q2 | chr10 97.8 Mb | 92.9–105.9 | 9.1 | +1.01 d ± 0.17 | **9.0** | 1.0 | Zd +1.18 (23), Zv +0.78 (24), Zx +0.94 (66) | [197 genes](data/introgression_fix1/qtl_genes/DTS_qtl_genes.tsv) | 93.3–97.8 (43) | 0.02 |

Top windows: the fewest windows holding half of the QTL's fine-mapping probability (section 3.4.3), shown from the first to the last (they need not be next to each other). Every gene in each interval, with its GO annotation, the k-mer GWAS −log10 p of its window and its window's PIP: [`data/introgression_fix1/qtl_genes/DTS_qtl_genes.tsv`](data/introgression_fix1/qtl_genes/DTS_qtl_genes.tsv).

</details>

<details><summary>Check: DTS in each year separately (2023 and 2025)</summary>

To see whether the DTS QTL hold up in each year, the whole analysis (inputs, LOCO GridLMM scan, 1,000 permutations, joint QTL model) was rerun on each year's SpATS BLUEs alone (`submit_year_check.sh`). The DTS BLUEs of the two years correlate at r = 0.69.

![DTS: k-mer ancestry GWAS on the combined BLUEs and on each year](figures/fix1/DTS_by_year.png)

GridLMM −log10 p of the window (effect per teosinte copy, days); \* = above that run's threshold within 2 Mb (thresholds 5.50 both years, 5.12 in 2023, 4.33 in 2025). The joint-model column says where the locus holds (Bonferroni over that run's QTL).

| Locus | Both years | 2023 only | 2025 only | Holds in the joint model |
|---|---|---|---|---|
| **chr10 94.6–97.8 Mb (*ZmCCT10*)** | 9.1 (+1.03)\* | 12.4 (+1.41)\* | 12.5 (+0.83)\* | both years, 2023, 2025 (at 95.6 Mb) |
| **chr5 221.8 Mb** | 9.4 (+1.00)\* | 9.2 (+1.15)\* | 4.8 (+0.47)\* | both years, 2023; 2025 at 225.4 Mb |
| chr5 225.4 Mb | 6.6 (+0.82)\* | 5.7 (+0.89)\* | 5.5 (+0.50)\* | 2025 |
| chr9 7.3 Mb | 3.7 (−0.49) | 2.9 (−0.50) | 9.4 (−0.55)\* | 2025 |
| chr8 133.8 Mb | 2.2 (+0.44) | 0.7 (+0.23) | 6.8 (+0.57)\* | 2025 |
| chr7 184.7 Mb (*dlf1*) | 2.4 (+0.37) | 2.4 (+0.45) | 4.8 (+0.38)\* | 2025 (taxa differ, −log10 p 4.8) |
| chr1 3.8 Mb | 2.2 (+0.40) | 2.1 (+0.46) | 5.0 (+0.43)\* | 2025 |
| chr7 163.1 Mb | 3.0 (−0.49) | 3.4 (−0.61) | 4.3 (−0.40)\* | 2025 |
| chr4 184.6 Mb | 2.5 (−0.44) | 1.2 (−0.34) | 4.5 (−0.42)\* | 2025 |
| chr2 226.8 Mb | 3.3 (+0.47) | 5.9 (+0.76)\* | 0.0 (−0.01) | 2023 |

- **Both QTL of the combined analysis are found in each year.** *ZmCCT10* is significant in both years, with the teosinte allele flowering later in both (+1.41 d in 2023, +0.83 d in 2025). The chr5 221.8 Mb QTL is strong in 2023 and above the threshold in 2025, but its effect is smaller in 2025 (+0.47 vs +1.15 d), and the 2025 joint model places the chr5 QTL 3.6 Mb away, at 225.4 Mb, where the taxa do not differ.
- **Other QTL are found in one year only.** The strongest, chr9 7.3 Mb (2025), has the same direction in 2023 (teosinte earlier by about half a day) but is not significant there. The 2025 QTL at chr7 184.7 Mb holds *dlf1* (*delayed flowering1*, Zm00001eb331880, 184.53 Mb). The 2025 QTL at chr8 133.8 Mb (interval 127.9–136.1 Mb) is 7 Mb from the florigen gene *ZCN8* (Zm00001eb353250, 126.68 Mb), outside its interval.
- **Why single years find more QTL than the combined BLUEs:** each run has its own threshold, and 2025's is lower (4.33, against 5.12 and 5.50). The 2025 scan is also more inflated (λ_GC 7.1, against 3.8 in 2023 and 4.1 combined). Effects that differ between years (genotype × year) partly cancel in the combined BLUEs. So QTL seen in one year only are weaker evidence than the two seen in every run.

Results: `data/introgression_fix1/assoc_gridlmm_spats_y2023`, `_y2025` and `joint_qtl/noperm_y2023`, `_y2025` (all traits). Figure: `final_script/31b_year_compare.py`.

</details>

#### 3.6.2 Days to anthesis (DTA)

![DTA: k-mer ancestry GWAS and joint QTL model](figures/fix1/joint_qtl_DTA.png)

Six QTL pass the threshold (−log10 p > 5.00); three hold in the joint model (Bonferroni over 6: −log10 p > 2.08). The first model had three, all holding; two of them are the same here (chr5 221.8 Mb and *ZmCCT10*), and the third (chr9 10.9 Mb, *MADS74*) is still above the scan threshold but no longer holds once the new chr9 1.4 Mb QTL is in the model.

| QTL | Known genes | Fine-mapping PIP | Effect (joint) | Joint −log10 p | Taxa differ | What it says |
|---|---|---|---|---|---|---|
| **Q6 chr10 101.2 Mb** | Zm00001eb418700 (ZmCCT10) | 0.009 (**top**) | +0.77 d | 4.4 | 1.9 (no) | No significant difference between taxa (Zd +1.07, Zv +0.44, Zx +0.54). Same region as DTS Q2 (ZmCCT10 at 96.2 Mb is inside the interval, 93.0–106.0 Mb); the strongest DTA scan peak (9.4). |
| **Q1 chr1 84.5 Mb** | — | 0.03 (best window; top windows 75.3–85.8 Mb, 55 genes) | +0.71 d | 5.3 | 17.2 (yes) | New in this model. Differs strongly between taxa: almost all *luxurians* (Zl +3.47; Zh +0.96, Zv +0.59, Zx +0.16, Zd −0.09). Suggestive for DTS (5.3). |
| **Q2 chr5 221.8 Mb** | Zm00001eb257680 (HDR1-like) | 0.429 (**top**) | +0.62 d | 5.0 | 11.2 (yes) | Differs between taxa (Zd +1.48, Zl +2.95, Zv +0.14, Zx +0.20). Same QTL as DTS Q1 (Zm00001eb257680), mostly *luxurians* and *diploperennis*. It delays pollen shedding less than silking (DTS +1.00 d), so it also lengthens ASI. |
| Q4 chr9 10.9 Mb | Zm00001eb372960 (MADS74) | 0.011 (CS) | −0.32 d | 1.3 | 0.5 (no) | Teosinte flowers earlier in every taxon (Zd −0.67, Zh −0.27, Zl −1.40, Zv −0.38, Zx −0.13). Above the scan threshold (5.5) but does not hold in the joint model; it shares tracts with Q3 (chr9 1.4 Mb, also earlier, also not holding). |

**Candidate genes:**
- **Q2, chr5 221.8 Mb, shared with DTS Q1:** Zm00001eb257680, the ortholog of rice *HEADING DATE REPRESSOR 1* (see DTS). *Fine-mapping (3.4.3): PIP 0.429 (the highest in the interval); in the top windows (10 genes).*
- **Q6, chr10 101.2 Mb, the DTS Q2 region:** Zm00001eb418700, *ZmCCT10* (see DTS). *Fine-mapping (3.4.3): PIP 0.009; in the top windows, but they are wide (94.1–104.6 Mb, 78 genes), so fine-mapping does not separate its genes.*
- **Q4, chr9 10.9 Mb (above the scan threshold for DTA, suggestive for DTS; does not hold in the joint model):** Zm00001eb372960 (*MADS74*; 8.36 Mb, inside the support interval 6.1–13.7 Mb). The teosinte allele makes plants shed pollen earlier, alike in every taxon. *MADS74* is an *SVP*-like MADS-box transcription factor: its Arabidopsis relative in the Gramene gene tree is *SHORT VEGETATIVE PHASE* (*SVP*, AT2G22540), and its annotation names the rice SVP-like MADS-box transcription factors 22 and 55. *SVP* is a flowering repressor, so an earlier-flowering teosinte allele fits a weaker repressor. In the first model this QTL held (joint −log10 p 5.6); here it does not (1.3), so it is a weaker candidate than before. *Fine-mapping (3.4.3): PIP 0.011; in the 95% credible set, but outside the top windows (9.5–13.7 Mb, 102 genes).*
  - SVP delays flowering by repressing the florigen gene *FT*, binding directly to the CArG motifs in its sequence. It also links flowering to temperature: plants flower later at 16 °C than at 23 °C, but *svp* mutants flower at almost the same time at both (Lee et al. 2007). At 16 °C, SVP switches on a second repressor, *TEMPRANILLO2*, which also represses *FT* and *TWIN SISTER OF FT* (Marín-González et al. 2015).
  - Natural alleles of *SVP* change flowering time: the Arabidopsis flowering QTL *FAQ1* is *SVP*, where one amino-acid change makes a loss-of-function, early-flowering allele, and how large its effect is depends on the genetic background (Méndez-Vigo et al. 2013).
  - *mads73* (ASI Q9, chr5 217.8 Mb) is a second *SVP*-like gene among our QTL.

<details><summary>All DTA QTL with per-taxon effects and all genes in the interval</summary>

| QTL | Position | Interval (Mb) | Scan −log10 p | Effect per copy (joint) | Joint −log10 p | Taxa differ −log10 p | Per-taxon effects (carrier plots) | Genes in interval | Top windows, Mb (genes) | Best window PIP |
|---|---|---|---|---|---|---|---|---|---|---|
| Q1 | chr1 84.5 Mb | 73.8–85.9 | 5.3 | +0.71 d ± 0.15 | **5.3** | **17.2** | Zd -0.09 (14), Zh +0.96 (6), Zl +3.47 (10), Zv +0.59 (19), Zx +0.16 (59) | [233 genes](data/introgression_fix1/qtl_genes/DTA_qtl_genes.tsv) | 75.3–85.8 (55) | 0.03 |
| Q2 | chr5 221.8 Mb | 221.7–222.3 | 5.4 | +0.62 d ± 0.14 | **5.0** | **11.2** | Zd +1.48 (10), Zl +2.95 (13), Zv +0.14 (47), Zx +0.20 (51) | [26 genes](data/introgression_fix1/qtl_genes/DTA_qtl_genes.tsv) | 221.7–221.9 (10) | 0.43 |
| Q3 | chr9 1.4 Mb | 0.9–2.2 | 5.6 | -0.37 d ± 0.17 | 1.5 | 0.5 | Zd -0.69 (21), Zh -0.27 (10), Zl -1.40 (22), Zv -0.09 (35), Zx -0.47 (84) | [12 genes](data/introgression_fix1/qtl_genes/DTA_qtl_genes.tsv) | 0.9–1.5 (1) | 0.11 |
| Q4 | chr9 10.9 Mb | 6.1–13.7 | 5.5 | -0.32 d ± 0.16 | 1.3 | 0.5 | Zd -0.67 (24), Zh -0.27 (8), Zl -1.40 (22), Zv -0.38 (56), Zx -0.13 (94) | [203 genes](data/introgression_fix1/qtl_genes/DTA_qtl_genes.tsv) | 9.5–13.7 (102) | 0.02 |
| Q5 | chr10 69.2 Mb | 67.3–75.6 | 6.0 | +0.24 d ± 0.20 | 0.7 | 0.5 | Zd +0.47 (25), Zv +0.17 (20), Zx +0.32 (52) | [115 genes](data/introgression_fix1/qtl_genes/DTA_qtl_genes.tsv) | – | 0.00 |
| Q6 | chr10 101.2 Mb | 93.0–106.0 | 9.4 | +0.77 d ± 0.19 | **4.4** | 1.9 | Zd +1.07 (24), Zv +0.44 (24), Zx +0.54 (64) | [194 genes](data/introgression_fix1/qtl_genes/DTA_qtl_genes.tsv) | 94.1–104.6 (78) | 0.02 |

Every gene in each interval, with its GO annotation, the k-mer GWAS −log10 p of its window and its window's PIP: [`data/introgression_fix1/qtl_genes/DTA_qtl_genes.tsv`](data/introgression_fix1/qtl_genes/DTA_qtl_genes.tsv).

</details>

#### 3.6.3 Anthesis–silking interval (ASI)

![ASI: k-mer ancestry GWAS and joint QTL model](figures/fix1/joint_qtl_ASI.png)

ASI = anthesis–silking interval, DTS − DTA. 10 regions, 6 hold in the joint model (Bonferroni over 10: −log10 p > 2.30). Many other chr2 and chr3 peaks drop in the joint model: they share one signal through linked tracts. Five of the six are the QTL of the first model; the sixth, chr4 67.0 Mb, replaces the first model's chr4 52.0 Mb QTL (*tga1*), as the lead of the same long chr4 region (below).

| QTL | Known genes | Fine-mapping PIP | Effect (joint) | Joint −log10 p | Taxa differ | What it says |
|---|---|---|---|---|---|---|
| **Q10 chr7 165.6 Mb** | Zm00001eb323880 (bhlh90) | 0.083 (CS) | -0.36 d | 13.5 | 0.3 (no) | No significant difference between taxa (Zd -0.47, Zl -0.57, Zv -0.36, Zx -0.30). The strongest ASI QTL: teosinte shortens the interval. |
| **Q9 chr5 217.8 Mb** | Zm00001eb255670 (mads73) | 0.273 (**top**) | +0.34 d | 9.1 | 3.4 (yes) | Differs between taxa (Zl +1.03, Zv +0.24, Zx +0.26). Lengthens ASI, mostly *luxurians*; next to the DTS/DTA chr5 221.8 Mb QTL. |
| **Q1 chr1 4.8 Mb** | Zm00001eb001510 (remorin) | 0.050 (CS) | +0.25 d | 6.5 | 0.0 (no) | No significant difference between taxa (Zd +0.18, Zh +0.22, Zl +0.14, Zv +0.32, Zx +0.23). Lengthens ASI. |
| **Q7 chr4 67.0 Mb** | — | 0.049 (best window; top windows 66.9–68.0 Mb, 8 genes) | -0.29 d | 5.2 | 0.2 (no) | No significant difference between taxa (Zh -0.52, Zl -0.96, Zv -0.37, Zx -0.23; Zd +0.03). Shortens ASI. The lead of the long chr4 region that held *tga1* in the first model. |
| **Q8 chr4 206.2 Mb** | Zm00001eb201580 (ZmFKF1b) | 0.004 (CS) | -0.21 d | 4.8 | 0.3 (no) | No significant difference between taxa (Zd -0.36, Zl +0.52, Zv -0.28, Zx -0.15). Shortens ASI. |
| **Q2 chr2 26.1 Mb** | Zm00001eb076990 (ckx12); Zm00001eb077080 (GA2ox2) | 0.012 (CS); 0.015 (CS) | -0.16 d | 3.7 | 0.1 (no) | No significant difference between taxa (Zd -0.17, Zh -0.51, Zl -0.16, Zv -0.20, Zx -0.13). Shortens ASI. |

**Candidate genes:**
- **Q10, chr7 165.6 Mb (the strongest ASI QTL):** Zm00001eb323880 (*bhlh90*, bHLH transcription factor 90; 165.74 Mb, in a window with −log10 p 12.6 next to the top window). MaizeGDB annotates it as the ortholog of Arabidopsis *FLOWERING BHLH 4* (Neupane et al. 2025); Gramene lists Arabidopsis *FBH3* (AT1G51140) as its closest homolog and rice *OsFBH1* / *POH1* (Os08g0506700) as its closest rice gene. In Arabidopsis the four FBH proteins bind E-boxes in the *CONSTANS* promoter and activate *CO*: their overexpression causes early flowering regardless of photoperiod, and *CO* levels drop in the *fbh* quadruple mutant (Ito et al. 2012). In rice, POH1 forms a dimer with HBP1, both bind the *Hd1* promoter and activate *Hd1*; their overexpression delays flowering by suppressing *Ehd1*, *Hd3a* and *RFT1* (Yin et al. 2023). *Fine-mapping (3.4.3): PIP 0.083 (highest window 0.30); in the 95% credible set, but outside the top windows (165.5–165.7 Mb, 4 genes).*
- **Q9, chr5 217.8 Mb:** Zm00001eb255670 (*mads73*, MADS-box transcription factor; 217.89 Mb, in a window with −log10 p 9.5). In its Gramene gene tree, the Arabidopsis gene is *SHORT VEGETATIVE PHASE* (*SVP*, AT2G22540; also *AGL22*), a MADS-box repressor of flowering: *Fine-mapping (3.4.3): PIP 0.273 (the highest in the interval); in the top windows (217.8–218.0 Mb, 4 genes), so the fine-mapping now points to it more clearly than in the first model (0.136, 16 genes).*
  - SVP represses the florigen gene *FT* by binding the CArG motifs in its sequence. It mediates the response of flowering to ambient temperature: plants flower later at 16 °C than at 23 °C, but *svp* mutants flower at almost the same time at both temperatures (Lee et al. 2007).
  - At low temperature (16 °C), SVP also switches on another floral repressor, *TEMPRANILLO2* (*TEM2*), which represses *FT* and *TWIN SISTER OF FT* (Marín-González et al. 2015).
  - Natural variation in *SVP* changes flowering time: the Arabidopsis flowering QTL *FAQ1* is *SVP*, where a single amino-acid change makes an early-flowering allele, and its effect depends on the genetic background (Méndez-Vigo et al. 2013).
- **Q1, chr1 4.8 Mb:** Zm00001eb001510 (remorin family protein; Zm00001d027403 in B73 v4; 4.45 Mb, in a window with −log10 p 5.8). In a maize panel, Wang et al. (2026) found a chr1 segment (v4 4,418,001–4,420,000 bp) both by GWAS for days to silking and days to tasseling and by a genetic-differentiation scan between subpopulations (Fst 0.53), and listed its seven genes, including this one, among 41 flowering-time candidates (their Table 2). In the same segment they point to a NAC transcription factor (*NAC78*, Zm00001d027395) and a histone H4 gene (Zm00001d027394). *Fine-mapping (3.4.3): PIP 0.050 (highest window 0.09); in the 95% credible set, but outside the top windows (3.9–4.8 Mb, 29 genes; in the first model it was in them).*
- **Q8, chr4 206.2 Mb:** Zm00001eb201580 (*ZmFKF1b*, *FLAVIN-BINDING, KELCH REPEAT, F-BOX 1b*; Zm00001d053091 in B73 v4; 216.2 Mb, in a window with −log10 p 4.2). The gene is 10 Mb from the top window but inside the support interval (202.3–226.0 Mb). FKF1 is a blue-light receptor whose expression follows the circadian clock. It has a light-sensing LOV domain, an F-box (which lets it tag proteins for breakdown) and Kelch repeats. In *Arabidopsis* it stabilizes CONSTANS (CO) and raises its expression, and CO switches on the florigen gene *FT*. Chen et al. (2024) showed that the maize copy works in the same pathway: *Fine-mapping (3.4.3): PIP 0.004 (highest window 0.02); in the 95% credible set, but outside the top windows (206.2–218.4 Mb, 82 genes).*
  - CRISPR knockouts of *ZmFKF1b*, and of its paralog *ZmFKF1a* (Zm00001d007445 in v4), each flowered significantly later than the wild type under short days, both for anthesis and for silking.
  - Both proteins bind ZmCONZ1 (the maize CO homolog) and ZmGI1 (GIGANTEA1), and both raise the expression of *ZmCONZ1* and of the maize florigen gene *ZCN8*.
  - *ZmFKF1b* lies in a known flowering-time QTL and in regions earlier found to be under selection, and its allele frequencies shifted during Chinese maize breeding from the 1960s–70s to the 1980s–90s.
  - **The link to geography:** the coding sequence has three main haplotypes. Hap_Z58 (as in the inbred Zheng58) is found mostly in temperate maize, at up to 92%. Hap_C7 (as in Chang7-2) is spread evenly over temperate and tropical maize. In an F2 population, plants with Hap_C7 flowered earlier, and Hap_C7 raised *ZmCONZ1* expression more. The authors propose that Hap_C7 suits a wider range of climate zones and could help adapt tropical maize to temperate regions.

  In our lines the teosinte allele shortens ASI by 0.21 d per copy, and the taxa do not differ significantly (from −0.15 d in *mexicana* to −0.36 d in *diploperennis*; *luxurians* goes the other way, +0.52 d, but with only 6 carrier lines). At the *ZmFKF1b* window itself the signal is for ASI (−log10 p 4.2), not for DTS (0.5) or DTA (0.0). Our donors were collected from 80 to 2,713 m, in very different climates, and we asked (first model) whether this QTL's effect changes with the donor's environment: it does not change with elevation, annual mean temperature, minimum temperature of the coldest month, rainfall, rainfall seasonality, aridity or soil pH (all p ≥ 0.20; `old_scripts/34_env_tests.R`, results in `data/env/tests/B_qtl_effect_x_env.tsv`). So in our population the teosinte allele acts the same way wherever the donor came from. The k-mer data can't tell which *ZmFKF1b* haplotype each donor carries, so the haplotype–climate link of Chen et al. (2024) can't be tested here.
- **Q2, chr2 26.1 Mb:** Zm00001eb076990 (*ckx12*, cytokinin oxidase/dehydrogenase 12; Zm00001d002989 in B73 v4; 29.06 Mb, in a window with −log10 p 6.3). CKX enzymes break down cytokinin, the hormone that keeps meristem cells dividing; in rice, lower expression of the homolog *OsCKX2* (*Gn1a*) lets cytokinin build up in the inflorescence meristem and increases grain number (Ashikari et al. 2005). Ni et al. (2024) mapped kernel row number in an advanced backcross of a two-rowed *mexicana* teosinte (MT1, donor) into the inbred Mo17: *Fine-mapping (3.4.3): PIP 0.012 (highest window 0.05); in the 95% credible set, but outside the top windows (25.3–28.2 Mb, 33 genes).*
  - The chr2 short-arm QTL *qKRN2* split into two closely linked QTL, *qKRN2-1* (fine-mapped to 272 kb) and *qKRN2-2* (775 kb), which explain 33.5% and 35.3% of the variance in kernel row number in their populations. The Mo17 alleles add about 1.4 and 1.2 rows.
  - Of the genes in the *qKRN2-1* interval, Zm00001d002989 (this gene) was expressed 3.7 times more in young ears (3–4 mm) of the line carrying the teosinte allele than in the line carrying the Mo17 allele. RNA in situ hybridization placed its expression in the primordia of the inflorescence meristem and of the spikelet-pair meristems, which set the number of kernel rows. They propose it as the causal gene: more CKX means less cytokinin, fewer spikelet-pair meristems and fewer rows.

  In our lines the teosinte allele at Q2 shortens ASI by 0.16 d per copy, and the taxa do not differ significantly. Why a gene acting in the developing ear would change ASI isn't known; we did not score kernel row number. The same interval also holds EH Q1 (+2.06 cm), and Prolif Q4–Q5 lie on either side of it (see [Regions shared between traits](#3613-regions-shared-between-traits)).
- **No longer a QTL of its own: chr4 52.0 Mb, Zm00001eb175150 (*tga1*, *teosinte glume architecture1*; 46.65 Mb).** In the first model this was ASI Q9 (joint −log10 p 3.0). The scan signal is still there (−log10 p 6.1 at *tga1*'s window and at 52.0 Mb), but chr4 is significant over a long stretch, and its lead is now at 67.0 Mb (Q7, interval 65.0–71.3 Mb), which does not include *tga1*; the conditional scan does not add a second QTL at 52 Mb. *tga1* is a classic maize domestication gene: it is the reason maize kernels are naked, while teosinte kernels are sealed inside a hard fruitcase made of a hardened glume and the rachis.
  - *tga1* encodes an SBP-box transcription factor. Maize and teosinte differ at a single fixed site in the gene, which changes one amino acid (Lys in teosinte, Asn in maize). This change turns maize TGA1 into a transcriptional repressor. Both forms bind the same GTAC motif in DNA, but the maize form makes more stable dimers. Because it is the only fixed difference, it is most likely the change that made naked kernels (Wang H et al. 2015).
  - TGA1 sits in a gene network with the other big domestication gene, *teosinte branched1* (*tb1*): TB1 binds the *tga1* promoter and activates it. TGA1 in turn regulates MADS-box genes, a family that sets the identity of floral and inflorescence organs: of eight tested, five were expressed more in immature ears with the teosinte allele of *tga1*, and TGA1 bound the promoters of five. One of them is *zmm19*, the *Tunicate* locus, whose overexpression makes the glumes grow around each kernel. Several genes in the network show signatures of selection during domestication (Studer et al. 2017).

  Since the Lys/Asn difference is fixed between maize and teosinte, every donor tract at *tga1* should carry the teosinte allele. In the first model, the teosinte segment shortened ASI by 0.24 d per copy. How *tga1* would change ASI isn't known: it acts in the developing ear, so the effect could come through ear development. Whether the chr4 signal is one QTL (at 67 Mb) or two (also near *tga1*) is not resolved by our data.
- **Not in a QTL interval in this model:** Zm00001eb132420 (*LNK2*, a circadian clock gene) was in a chr3 QTL of the first model that did not hold.

The chr3 95.8 Mb peak (−log10 p 4.9) lies inside the chr3 108.3 Mb region (93.8–121.6 Mb), and the conditional scan does not add it.

<details><summary>All ASI QTL with per-taxon effects and all genes in the interval</summary>

| QTL | Position | Interval (Mb) | Scan −log10 p | Effect per copy (joint) | Joint −log10 p | Taxa differ −log10 p | Per-taxon effects (carrier plots) | Genes in interval | Top windows, Mb (genes) | Best window PIP |
|---|---|---|---|---|---|---|---|---|---|---|
| Q1 | chr1 4.8 Mb | 3.6–5.6 | 6.1 | +0.25 d ± 0.05 | **6.5** | 0.0 | Zd +0.18 (15), Zh +0.22 (10), Zl +0.14 (12), Zv +0.32 (36), Zx +0.23 (64) | [85 genes](data/introgression_fix1/qtl_genes/ASI_qtl_genes.tsv) | 3.9–4.8 (29) | 0.09 |
| Q2 | chr2 26.1 Mb | 25.3–29.5 | 7.1 | -0.16 d ± 0.04 | **3.7** | 0.1 | Zd -0.17 (35), Zh -0.51 (5), Zl -0.16 (21), Zv -0.20 (48), Zx -0.13 (107) | [126 genes](data/introgression_fix1/qtl_genes/ASI_qtl_genes.tsv) | 25.3–28.2 (33) | 0.05 |
| Q3 | chr2 68.2 Mb | 60.7–78.5 | 5.0 | -0.09 d ± 0.05 | 1.2 | 0.0 | Zd -0.21 (33), Zh -0.12 (7), Zl -0.10 (15), Zv -0.06 (43), Zx -0.07 (88) | [178 genes](data/introgression_fix1/qtl_genes/ASI_qtl_genes.tsv) | 62.1–75.7 (72) | 0.01 |
| Q4 | chr2 191.2 Mb | 184.6–193.2 | 4.7 | -0.12 d ± 0.05 | 1.9 | 0.2 | Zd +0.01 (17), Zh -0.14 (6), Zl -0.40 (16), Zv -0.13 (47), Zx -0.12 (112) | [200 genes](data/introgression_fix1/qtl_genes/ASI_qtl_genes.tsv) | 190.8–193.2 (53) | 0.03 |
| Q5 | chr3 70.7 Mb | 61.9–75.4 | 4.9 | -0.02 d ± 0.15 | 0.1 | 0.5 | Zd -0.41 (25), Zh -0.17 (13), Zl -0.41 (8), Zv -0.33 (25), Zx -0.14 (67) | [88 genes](data/introgression_fix1/qtl_genes/ASI_qtl_genes.tsv) | – | 0.00 |
| Q6 | chr3 108.3 Mb | 93.8–121.6 | 5.2 | -0.20 d ± 0.15 | 0.8 | 0.6 | Zd -0.37 (25), Zh -0.17 (13), Zl -0.41 (8), Zv -0.38 (28), Zx -0.14 (72) | [285 genes](data/introgression_fix1/qtl_genes/ASI_qtl_genes.tsv) | – | 0.00 |
| Q7 | chr4 67.0 Mb | 65.0–71.3 | 6.6 | -0.29 d ± 0.07 | **5.2** | 0.2 | Zd +0.03 (5), Zh -0.52 (8), Zl -0.96 (5), Zv -0.37 (29), Zx -0.23 (42) | [91 genes](data/introgression_fix1/qtl_genes/ASI_qtl_genes.tsv) | 66.9–68.0 (8) | 0.05 |
| Q8 | chr4 206.2 Mb | 202.3–226.0 | 5.2 | -0.21 d ± 0.05 | **4.8** | 0.3 | Zd -0.36 (17), Zl +0.52 (6), Zv -0.28 (32), Zx -0.15 (77) | [420 genes](data/introgression_fix1/qtl_genes/ASI_qtl_genes.tsv) | 206.2–218.4 (82) | 0.02 |
| Q9 | chr5 217.8 Mb | 217.4–218.8 | 9.4 | +0.34 d ± 0.06 | **9.1** | **3.4** | Zl +1.03 (13), Zv +0.24 (43), Zx +0.26 (42) | [56 genes](data/introgression_fix1/qtl_genes/ASI_qtl_genes.tsv) | 217.8–218.0 (4) | 0.27 |
| Q10 | chr7 165.6 Mb | 165.2–166.1 | 13.1 | -0.36 d ± 0.05 | **13.5** | 0.3 | Zd -0.47 (17), Zl -0.57 (14), Zv -0.36 (45), Zx -0.30 (78) | [26 genes](data/introgression_fix1/qtl_genes/ASI_qtl_genes.tsv) | 165.5–165.7 (4) | 0.30 |

Bonferroni over the 10 QTL: −log10 p > 2.30 (bold). Top windows: the fewest windows holding half of the QTL's fine-mapping probability (section 3.4.3). Every gene in each interval, with its GO annotation, the k-mer GWAS −log10 p of its window and its window's PIP: [`data/introgression_fix1/qtl_genes/ASI_qtl_genes.tsv`](data/introgression_fix1/qtl_genes/ASI_qtl_genes.tsv).

</details>

#### 3.6.4 Nodes with brace roots (NBR)

![NBR: k-mer ancestry GWAS and joint QTL model](figures/fix1/joint_qtl_NBR.png)

Number of stem nodes that carry brace roots. 8 regions, 3 hold (Bonferroni over 8: −log10 p > 2.20). The first model had 4 holding; chr10 98.2 Mb and chr8 126.2 Mb are the same, chr9 4.2 Mb is new, and chr1 260.8 Mb (2.1) and chr9 140 Mb no longer hold.

| QTL | Known genes | Fine-mapping PIP | Effect (joint) | Joint −log10 p | Taxa differ | What it says |
|---|---|---|---|---|---|---|
| **Q8 chr10 98.2 Mb** | Zm00001eb418930 | 0.003 (outside CS) | +0.18 | 9.4 | 2.7 (yes) | Differs between taxa (Zd +0.29, Zv +0.23, Zx +0.08). The strongest NBR QTL; smallest in *mexicana*. |
| **Q6 chr8 126.2 Mb** | — | 0.15 (best window; top windows 126.2–126.7 Mb, 14 genes) | +0.12 | 4.4 | 1.1 (no) | No significant difference between taxa (Zd +0.10, Zh +0.33, Zl +0.15, Zv +0.09, Zx +0.09). |
| **Q7 chr9 4.2 Mb** | — | 0.07 (best window; top windows 4.0–9.5 Mb, 49 genes) | −0.09 | 4.4 | 3.2 (yes) | Fewer brace-root nodes; differs between taxa, mostly *diploperennis* (Zd −0.39; Zh −0.09, Zl −0.07, Zv −0.06, Zx −0.07). |

**Candidate gene:**
- **Q8, chr10 98.2 Mb (the strongest NBR QTL):** Zm00001eb418930 (GRMZM2G174736 in B73 v3; 97.41 Mb, in a window with −log10 p 8.0; the top window has 9.3). Sun et al. (2022) ran a GWAS of brace root angle and brace root diameter in 508 maize inbreds grown in three environments (Liaoning 2016 and 2017, Jilin 2016). Their significant chr10 SNP for brace root diameter (chr10.S_95437751, in older B73 coordinates) explains 13.95% of the variance, and GRMZM2G174736 is one of their candidate genes for brace root diameter:
  - They annotate it as a hydroxyproline-rich glycoprotein (HRGP, also called extensin), a structural protein of the cell wall. HRGPs help assemble the cell wall, are expressed in specific tissues and stages, and some act in root development.
  - It is one of three cell-wall candidates they highlight for brace root diameter, with an expansin (GRMZM2G445169, which loosens the cell wall) and a homolog of the *Arabidopsis* receptor kinase FEI (GRMZM2G479243). They propose that these genes shape brace roots through the cell wall.
  - In public RNA-seq of 13 tissues, including brace roots, it is expressed at moderate to high levels, but not only in brace roots. Its function has not been tested.

  The current B73 v5 annotation lists the gene as uncharacterized (protein family: SICKLE).

  In our lines the teosinte allele adds 0.18 brace-root nodes per copy (joint −log10 p 9.4), and the effect differs between taxa: +0.29 in *diploperennis*, +0.23 in *parviglumis* and +0.08 in *mexicana* (in the first model *diploperennis* was the smallest, +0.05).

  **Is this just flowering time?** The QTL is 2 Mb from *ZmCCT10* (Zm00001eb418700, 96.18 Mb), whose teosinte allele delays flowering (DTS Q2, DTA Q6), and a later-flowering plant can make more nodes. Two results argue against that:
  - NBR is hardly correlated with flowering time in our lines (r = 0.10 with DTS, 0.11 with DTA). In a linear model with taxon (1,243 pools), adding DTS as a covariate barely changes the Q8 effect (+0.18 → +0.16 nodes per copy, p = 5 × 10⁻⁸). Adding DTS, DTA and plant height together lowers it to +0.11 (p = 2 × 10⁻⁴), mostly because of plant height (r = 0.34 with NBR), which has its own QTL in this region (PH Q15).
  - The flowering effect of this region is the same in all taxa (DTS Q2: taxa differ −log10 p 1.0), but the NBR effect differs between taxa (2.7).
  - The position no longer helps: in this model the DTS peak is at 97.8 Mb, 0.5 Mb from the NBR peak (3.7 Mb in the first model).
  - **Fine-mapping does not favour this gene.** Its window has PIP 0.003 and lies outside the 95% credible set; the QTL's probability is spread over 98.0–99.4 Mb (highest window 0.05; 27 genes in the top windows, listed with `in_top_windows = yes` in `data/introgression_fix1/qtl_genes/NBR_qtl_genes.tsv`). This lowers the support from our own data, but does not remove the gene: the case for it comes from the brace-root GWAS of Sun et al. (2022), not from the position of our peak, and a low PIP does not rule a gene out (section 3.4.3). The genes in the top windows are equally worth a look.

  So the brace-root effect is not explained by flowering time. Our tracts are too long to separate Zm00001eb418930 from *ZmCCT10*, though: only 1 of the 1,243 pools has a different state at the two genes (117 carry a donor segment there). Script: `final_script/nbr_q8_flowering_check.py`; output `data/introgression_fix1/checks/nbr_q8.txt`.

The five chr1 peaks (all with negative effects) share their signal and do not hold individually. **Q4 chr1 219.8 Mb** has no overall effect but differs between taxa (joint −log10 p 0.2, taxa differ 2.5: *luxurians* carriers have fewer brace-root nodes, −0.86). *ZmD14* (Zm00001eb007670) is in the interval of Q2 (chr1 23.9 Mb; joint 1.8, taxa differ 1.9), which does not hold. The first model's chr10 25.4 Mb QTL (only *diploperennis* carriers, taxa differ 6.1) is not among this model's peaks.

<details><summary>All NBR QTL with per-taxon effects and all genes in the interval</summary>

| QTL | Position | Interval (Mb) | Scan −log10 p | Effect per copy (joint) | Joint −log10 p | Taxa differ −log10 p | Per-taxon effects (carrier plots) | Genes in interval | Top windows, Mb (genes) | Best window PIP |
|---|---|---|---|---|---|---|---|---|---|---|
| Q1 | chr1 12.2 Mb | 8.7–12.4 | 4.3 | -0.07 ± 0.03 | 1.9 | 1.2 | Zd -0.04 (16), Zh -0.14 (10), Zl +0.06 (6), Zv -0.03 (36), Zx -0.07 (59) | [127 genes](data/introgression_fix1/qtl_genes/NBR_qtl_genes.tsv) | 8.9–12.3 (63) | 0.04 |
| Q2 | chr1 23.9 Mb | 21.9–25.7 | 4.0 | -0.07 ± 0.03 | 1.8 | 1.9 | Zd +0.10 (14), Zh -0.18 (9), Zl -0.10 (12), Zv -0.15 (28), Zx -0.06 (56) | [102 genes](data/introgression_fix1/qtl_genes/NBR_qtl_genes.tsv) | 23.3–24.9 (30) | 0.05 |
| Q3 | chr1 210.3 Mb | 198.3–216.0 | 4.3 | -0.06 ± 0.04 | 1.1 | 1.3 | Zd -0.22 (8), Zh -0.26 (5), Zl +0.12 (8), Zv -0.01 (30), Zx -0.08 (78) | [444 genes](data/introgression_fix1/qtl_genes/NBR_qtl_genes.tsv) | – | 0.01 |
| Q4 | chr1 219.8 Mb | 216.2–222.5 | 4.5 | -0.02 ± 0.04 | 0.2 | **2.5** | Zd -0.18 (9), Zh -0.14 (5), Zl -0.86 (7), Zv -0.00 (33), Zx +0.04 (72) | [116 genes](data/introgression_fix1/qtl_genes/NBR_qtl_genes.tsv) | – | 0.00 |
| Q5 | chr1 260.8 Mb | 259.3–264.4 | 5.3 | -0.09 ± 0.04 | 2.1 | 1.5 | Zd -0.11 (13), Zh +0.17 (6), Zl +0.46 (9), Zv -0.16 (34), Zx -0.07 (77) | [134 genes](data/introgression_fix1/qtl_genes/NBR_qtl_genes.tsv) | 260.4–262.7 (42) | 0.03 |
| Q6 | chr8 126.2 Mb | 125.7–127.6 | 4.2 | +0.12 ± 0.03 | **4.4** | 1.1 | Zd +0.10 (23), Zh +0.33 (10), Zl +0.15 (9), Zv +0.09 (32), Zx +0.09 (61) | [54 genes](data/introgression_fix1/qtl_genes/NBR_qtl_genes.tsv) | 126.2–126.7 (14) | 0.15 |
| Q7 | chr9 4.2 Mb | 2.1–9.7 | 4.1 | -0.09 ± 0.02 | **4.4** | **3.2** | Zd -0.39 (26), Zh -0.09 (10), Zl -0.07 (21), Zv -0.06 (46), Zx -0.07 (88) | [185 genes](data/introgression_fix1/qtl_genes/NBR_qtl_genes.tsv) | 4.0–9.5 (49) | 0.07 |
| Q8 | chr10 98.2 Mb | 93.0–99.6 | 9.3 | +0.18 ± 0.03 | **9.4** | **2.7** | Zd +0.29 (23), Zv +0.23 (23), Zx +0.08 (64) | [107 genes](data/introgression_fix1/qtl_genes/NBR_qtl_genes.tsv) | 98.0–99.4 (27) | 0.05 |

Bonferroni over the 8 QTL: −log10 p > 2.20 (bold). Top windows: the fewest windows holding half of the QTL's fine-mapping probability (section 3.4.3). Every gene in each interval, with its GO annotation, the k-mer GWAS −log10 p of its window and its window's PIP: [`data/introgression_fix1/qtl_genes/NBR_qtl_genes.tsv`](data/introgression_fix1/qtl_genes/NBR_qtl_genes.tsv).

</details>

#### 3.6.5 Leaf greenness after silking (SPAD; 2023 only)

![SPAD: k-mer ancestry GWAS and joint QTL model](figures/fix1/joint_qtl_SPAD2.png)

SPAD is the SPAD-meter reading 36 days after silking (2023 only; H2 0.40). The meter reads how green the leaf is, a proxy for its chlorophyll; at this stage the leaves begin to senesce, so a high value means the leaf keeps its chlorophyll longer (stay-green), and nitrogen moving out of the leaf to the grain lowers it. 8 regions, 3 hold (Bonferroni over 8: −log10 p > 2.20). The first model had 4: chr10 148 Mb, chr7 123 Mb and chr9 (now at 96.8 Mb) are the same; chr1 148.6 Mb no longer holds (below).

| QTL | Known genes | Fine-mapping PIP | Effect (joint) | Joint −log10 p | Taxa differ | What it says |
|---|---|---|---|---|---|---|
| **Q6 chr7 123.3 Mb** | Zm00001eb313460 (senescence-associated protein) | 0.039 (**top**) | −2.98 | 4.7 | 0.2 (no) | No significant difference between taxa (Zd −1.08, Zl −1.00, Zv −1.70, Zx −3.53). The gene is next to the top window. |
| **Q7 chr9 96.8 Mb** | — | 0.01 (best window) | +2.59 | 4.7 | 1.8 (no) | Greener leaves (Zd +2.35, Zl +3.16, Zv +0.58, Zx +3.54). The interval spans 30.7 Mb and 505 genes, too wide to point to one gene. |
| **Q8 chr10 148.2 Mb** | Zm00001eb432440 (FRS6-like); Zm00001eb432850 (chlorophyll a/b-binding protein) | 0.088 (**top**); 0.040 (**top**) | +2.15 | 4.1 | 0.3 (no) | No significant difference between taxa (Zd +1.73, Zh +2.43, Zl +3.37, Zv +0.82, Zx +2.99). Same region as the best PredN region. |

**Candidate genes:**
- **Q8, chr10 148.2 Mb:** Zm00001eb432440 (FAR1-RELATED SEQUENCE 6-like; 148.33 Mb). It lies in the top windows of the QTL (window −log10 p 4.3; the windows of the interval range from 2.9 to 4.4), next to a second FAR1-related gene (Zm00001eb432450). FAR1 and FHY3 are light-signalling transcription factors; in *Arabidopsis* they switch on *HEMB1*, an enzyme of chlorophyll synthesis, when seedlings first turn green (Tang et al. 2012). This maize gene is not characterized; the ATP-synthase GO terms it carries in MaizeGDB are automatic predictions and do not fit a FAR1 protein. The teosinte allele makes leaves greener (+2.15 SPAD units per copy, alike across taxa). *Fine-mapping (3.4.3): PIP 0.088; in the top windows (148.2–149.5 Mb, 52 genes).*
  - Second candidate: Zm00001eb432850 (chlorophyll a/b-binding protein, chloroplastic; 148.97 Mb, window 4.0). Chlorophyll a/b-binding (light-harvesting) proteins hold most of the leaf's chlorophyll in the antenna of photosystem II, so more of them means more chlorophyll per leaf area, which is what the SPAD meter reads. *Fine-mapping (3.4.3): PIP 0.040; in the top windows (148.2–149.5 Mb, 52 genes).*
- **The same region in predicted leaf nitrogen.** The end of chr10 is also the best region for leaf nitrogen predicted from leaf spectra (PredN, 2025; top window 146.85 Mb, inside the Q8 interval; −log10 p 2.8, below its threshold of 4.23; see [Other traits](#3615-other-traits)). As for SPAD, the teosinte allele goes up there.
- **Q6, chr7 123.3 Mb:** Zm00001eb313460 (senescence-associated protein; 123.47 Mb, −log10 p 4.5, next to the top window). The teosinte allele makes leaves less green (−2.98 per copy). *Fine-mapping (3.4.3): PIP 0.039; in the top windows (122.2–126.2 Mb, 56 genes).*
- **chr1: strong in the scan, but no QTL holds.** Much of chr1 from 35 to 164 Mb is above the scan threshold, all with the teosinte allele making leaves less green (−1.5 to −3.7 per copy). The joint model splits it into five QTL (Q1–Q5) that share tracts, and none holds on its own (joint −log10 p 0.0–1.4). In the first model, one of them (chr1 148.6 Mb, interval 144.5–163.4 Mb) held just above its cut (2.4 against 2.08); here that stretch is outside every QTL interval (its windows are at −log10 p 4.4–5.7, median 4.8). Four genes there act on chlorophyll, photosynthesis or amino-acid recycling; they are candidates for the chr1 signal as a whole, not for a QTL that holds:
  - Zm00001eb027950 (*rccr1*, *red chlorophyll catabolite reductase1*; 144.52 Mb, window −log10 p 4.5), an enzyme of the pathway that breaks chlorophyll down during senescence. In the maize NAM population it is the candidate gene at a locus for the first principal component of leaf metabolites, a trait that correlates most with glutamate, chlorophyll a and malate (r = 0.74, 0.62 and 0.61); in a predicted protein-interaction network its cluster is enriched for chlorophyll breakdown (Poretsky et al. 2024).
  - Zm00001eb029690 (*bsd2*, *bundle sheath defective2*; 162.00 Mb, window 5.0), a chloroplast protein with a DnaJ-like zinc-finger domain that is needed for Rubisco. *bsd2* mutants have no detectable Rubisco although the *rbcL* mRNA is abundant and on polysomes, so BSD2 acts after translation; the mutant plants are pale green and their bundle-sheath chloroplasts are abnormal (Brutnell et al. 1999). Without Rubisco the mutant fixes no carbon, yet PEPC kinase and PEPC in the mesophyll are still activated by light, so that signal does not need a working Calvin cycle (Smith et al. 1998).
  - Zm00001eb028550 (*ppd1*, PsbP domain-containing protein1; 152.01 Mb, window 4.8). In *Arabidopsis*, PPD1 is a protein of the thylakoid lumen that is essential for photosynthetic growth: it helps assemble photosystem I, and RNAi lines lose assembled photosystem I in step with the loss of PPD1 (Roose et al. 2014). It is one of the known photosystem I assembly factors, with Ycf3, Ycf4, PYG7 and Y3IP1 (Li M et al. 2026).
  - Zm00001eb028540 (cytosol aminopeptidase domain-containing protein, similar to *Arabidopsis* LAP2, AT4G30920; 151.94 Mb, window 4.8). In *Arabidopsis*, loss of the leucine aminopeptidase LAP2 slows growth, makes leaves senesce early and makes plants more sensitive to stress; the mutants keep free leucine by making more of it, which shifts glutamate metabolism and lowers GABA, both nitrogen-rich compounds (Waditee-Sirisattha et al. 2011a, b). LAP2 is also one of the chloroplast aminopeptidases of the peptide-degradation cascade that turns the targeting peptides cut off imported proteins back into free amino acids (Teixeira et al. 2016). Early senescence when it is lost fits a trait measured during senescence.

These candidates come from the gene annotation only; none has been tested for SPAD in maize, and SPAD was measured once, in one year. Every gene in each interval, with its GO terms: [`data/introgression_fix1/qtl_genes/SPAD2_qtl_genes.tsv`](data/introgression_fix1/qtl_genes/SPAD2_qtl_genes.tsv).

<details><summary>All SPAD QTL with per-taxon effects and all genes in the interval</summary>

| QTL | Position | Interval (Mb) | Scan −log10 p | Effect per copy (joint) | Joint −log10 p | Taxa differ −log10 p | Per-taxon effects (carrier plots) | Genes in interval | Top windows, Mb (genes) | Best window PIP |
|---|---|---|---|---|---|---|---|---|---|---|
| Q1 | chr1 40.0 Mb | 34.8–40.5 | 5.1 | -1.53 ± 0.76 | 1.4 | 0.2 | Zd +0.91 (23), Zh +2.74 (8), Zl -4.55 (9), Zv -1.70 (24), Zx -2.55 (53) | [149 genes](data/introgression_fix1/qtl_genes/SPAD2_qtl_genes.tsv) | 34.8–40.0 (59) | 0.03 |
| Q2 | chr1 57.9 Mb | 56.8–60.1 | 6.2 | -1.66 ± 1.02 | 1.0 | 0.1 | Zd -2.94 (19), Zh -0.04 (9), Zl -0.78 (7), Zv -0.27 (18), Zx -0.21 (45) | [62 genes](data/introgression_fix1/qtl_genes/SPAD2_qtl_genes.tsv) | – | 0.00 |
| Q3 | chr1 78.0 Mb | 76.3–105.0 | 4.7 | +0.27 ± 1.19 | 0.1 | 0.1 | Zd -0.30 (11), Zh -1.29 (5), Zl +0.76 (8), Zv -3.00 (14), Zx +0.26 (52) | [448 genes](data/introgression_fix1/qtl_genes/SPAD2_qtl_genes.tsv) | – | 0.02 |
| Q4 | chr1 106.2 Mb | 105.1–108.9 | 4.9 | -0.03 ± 1.50 | 0.0 | 0.1 | Zd +0.53 (12), Zh +3.27 (6), Zl +0.34 (6), Zv -1.24 (16), Zx -1.82 (45) | [48 genes](data/introgression_fix1/qtl_genes/SPAD2_qtl_genes.tsv) | 105.1–108.9 (24) | 0.03 |
| Q5 | chr1 128.8 Mb | 114.5–143.8 | 6.1 | -2.95 ± 1.44 | 1.4 | 0.2 | Zd -3.83 (13), Zh -7.69 (5), Zl -4.19 (7), Zv -2.55 (11), Zx -0.73 (44) | [197 genes](data/introgression_fix1/qtl_genes/SPAD2_qtl_genes.tsv) | 128.7–141.5 (64) | 0.01 |
| Q6 | chr7 123.3 Mb | 119.8–130.7 | 4.5 | -2.98 ± 0.70 | **4.7** | 0.2 | Zd -1.08 (7), Zl -1.00 (5), Zv -1.70 (37), Zx -3.53 (53) | [198 genes](data/introgression_fix1/qtl_genes/SPAD2_qtl_genes.tsv) | 122.2–126.2 (56) | 0.04 |
| Q7 | chr9 96.8 Mb | 86.3–117.0 | 4.6 | +2.59 ± 0.61 | **4.7** | 1.8 | Zd +2.35 (18), Zh -6.69 (5), Zl +3.16 (11), Zv +0.58 (43), Zx +3.54 (80) | [505 genes](data/introgression_fix1/qtl_genes/SPAD2_qtl_genes.tsv) | 93.2–112.9 (129) | 0.01 |
| Q8 | chr10 148.2 Mb | 146.7–150.2 | 4.4 | +2.15 ± 0.55 | **4.1** | 0.3 | Zd +1.73 (38), Zh +2.43 (5), Zl +3.37 (13), Zv +0.82 (32), Zx +2.99 (75) | [182 genes](data/introgression_fix1/qtl_genes/SPAD2_qtl_genes.tsv) | 148.2–149.5 (52) | 0.10 |

Bonferroni over the 8 QTL: −log10 p > 2.20 (bold). Top windows: the fewest windows holding half of the QTL's fine-mapping probability (section 3.4.3). Every gene in each interval, with its GO annotation, the k-mer GWAS −log10 p of its window and its window's PIP: [`data/introgression_fix1/qtl_genes/SPAD2_qtl_genes.tsv`](data/introgression_fix1/qtl_genes/SPAD2_qtl_genes.tsv).

</details>

#### 3.6.6 Plant height (PH)

![PH: k-mer ancestry GWAS and joint QTL model](figures/fix1/joint_qtl_PH.png)

15 regions, 10 hold in the joint model (Bonferroni over 15: −log10 p > 2.48). Effects are in cm per teosinte copy and alike across taxa (no "taxa differ" test passes). Seven of the first model's eight held QTL are found again (chr2, chr3, chr4, chr5 118 Mb, chr7, chr9, chr10 *ZmCCT10*); chr1 201.8 Mb is replaced by chr1 184.1 Mb (Q1, −3.83 cm), and two QTL are new: chr5 2.6 Mb (Q6, +3.05 cm; just below the threshold in the first model) and chr10 6.5 Mb (Q14, +3.27 cm, added by the conditional scan). The strongest are Q13 chr9 139.8 Mb (+4.64 cm, joint −log10 p 10.5) and Q12 chr7 137.8 Mb (below).

| QTL | Known genes | Fine-mapping PIP | Effect (joint) | Joint −log10 p | Taxa differ | What it says |
|---|---|---|---|---|---|---|
| **Q12 chr7 137.8 Mb** | Zm00001eb316740 (tsh4) | 0.017 (CS) | -5.10 cm | 10.2 | 0.3 (no) | No significant difference between taxa (Zd -5.01, Zh -5.94, Zl +0.38, Zv -4.58, Zx -5.78). Same region as SL Q10. |
| **Q3 chr3 13.2 Mb** | Zm00001eb123060 (ra2) | 0.057 (**top**) | -3.73 cm | 6.8 | 0.7 (no) | No significant difference between taxa (Zd -4.92, Zh +5.69, Zl -0.23, Zv -3.63, Zx -4.42). Same window as EH Q3. |
| **Q15 chr10 99.5 Mb** | Zm00001eb418700 (ZmCCT10) | 0.014 (CS) | +4.70 cm | 6.3 | 0.1 (no) | No significant difference between taxa (Zd +2.46, Zv +7.09, Zx +4.06). The ZmCCT10 region: later flowering, taller plants. |
| **Q4 chr4 198.4 Mb** | Zm00001eb197430 (ns2) | 0.027 (**top**) | +4.05 cm | 2.8 | 0.1 (no) | No significant difference between taxa (Zd +6.75, Zh +4.71, Zl +11.94, Zv +5.04, Zx +3.72). Same region as EH Q4, BL Q5 and SL Q3. |

**Candidate genes:**
- **Q12, chr7 137.8 Mb:** Zm00001eb316740 (*tsh4*, *tassel sheath4*; 138.54 Mb, in a window with −log10 p 9.9; the top window is at 137.75 Mb, −log10 p 10.5; in the first model it was at 139.05 Mb). *tsh4* encodes an SBP-box (SPL) transcription factor and is one of the genes that miR156 represses (the dominant *Corngrass1* mutant overexpresses miR156). It works in the developing tassel and ear: *Fine-mapping (3.4.3): PIP 0.017 (highest window 0.08); in the 95% credible set, but outside the top windows (137.5–139.1 Mb, 26 genes).*
  - **It suppresses the bract and promotes branching.** In the maize inflorescence, each branch forms with a small leaf (the bract) below it, and that leaf normally stays tiny. In *tsh4* mutants the bracts grow out into leaves and the tassel has fewer branches. TSH4 binds the promoters of *tsh1* (a GATA zinc-finger gene) and *liguleless2* and keeps them on in the suppressed bract; *tsh1*, *tsh4* and *lg2* act together to keep bracts small and branch meristems indeterminate, i.e. growing into long branches (Xiao et al. 2022).
  - **It works with its two closest relatives.** *unbranched2* and *unbranched3* act redundantly with *tsh4*: the triple mutant makes very few tassel branches, more tillers, and up to twice as many leaves above the top ear, yet flowers at the same time as the wild type. The three genes slow the rate at which new side primordia form, so the meristem has time to renew its stem cells (Chuck et al. 2014).
  - **It responds to shade.** UB2, UB3 and TSH4 are hubs of a network of over 2,400 genes during early tassel development. Simulated shade makes tassels and ears develop early, but the *ub2 ub3* double mutant and the *tsh4* single mutant hardly respond; the three proteins directly regulate *barren inflorescence2* (*BIF2*, auxin signaling) and *ZmTCP30* (Kong et al. 2023).
  - **It is the target of RA2.** RA2 (our Q3 candidate) binds the *tsh4* promoter and represses it; ZmELF3.1 makes this repression stronger (Xie et al. 2023). *tsh4* acts against the *ramosa* genes in branch determinacy and is epistatic to them: the *ra2 tsh4* double mutant has as few tassel branches as *tsh4* alone (Xie et al. 2023; Xiao et al. 2022).

  In our lines the teosinte allele at Q12 makes plants 5.10 cm shorter per copy, alike across taxa (−4.6 to −5.9 cm in four taxa; *luxurians* +0.4 cm with 9 carriers). Our plant height is measured from the ground to the tip of the tassel, so the tassel is part of it. None of these papers measured plant height or tassel length in *tsh4* mutants, so the link to PH is not shown; the leaves above the top ear (LAE), which the triple mutant has more of, show no signal at *tsh4* in our lines (−log10 p 1.0). Two other genes sit at the peak and can't be separated from *tsh4* with our tract lengths (all three are within 1.2 Mb of the top window, −log10 p 9.9–10.4), and in this model the fine-mapping favours them slightly over *tsh4*:
  - Zm00001eb316780 (138.64 Mb), a leucine-rich repeat protein that may be the maize homolog of *Arabidopsis FLOR1*, an LRR protein that binds MADS-box transcription factors and promotes flowering (Torti et al. 2012, as cited in Zhang et al. 2025). In a maize shoot transcriptome atlas, the maize *FLOR1* homolog rises in expression while the ear elongates fast (Zhang et al. 2025); its effect on plant height is not known. *Fine-mapping (3.4.3): PIP 0.061 (highest window 0.08); in the top windows (26 genes).*
  - Zm00001eb316880 (*wakl38*, wall-associated receptor kinase-like 38; 138.94 Mb). It is one of 13 candidate genes for ear shank length from BSA-seq and RNA-seq in a cross of a long-shank (L135) and a short-shank (WL134) inbred: it lies in their chr7 QTL *qESL10* and is expressed more in the long-shank parent (An et al. 2026). The two parents also differ in plant height, but the gene's function has not been tested. *Fine-mapping (3.4.3): PIP 0.084 (the highest in the interval); in the top windows (26 genes).*

  **Does *tsh4* interact with *ra2* in our lines?** Since RA2 represses *tsh4*, we tested whether the teosinte segments at the two genes (PH Q3 and Q12) interact for plant height. With the pool state at each gene (B73 = 0, SEG = 0.5, DONOR = 1) and taxon as covariates, in an ordinary linear model on the PH BLUEs of 1,257 pools, each segment shortens plants on its own (*ra2* −3.4 cm and *tsh4* −4.9 cm per copy, both p < 10⁻⁴), and the interaction is not significant (p = 0.27). The 24 pools that carry both segments are 15.9 cm shorter than the 971 pools that carry neither, more than the two single effects added up (5.3 + 5.3 cm), but not significantly. So the two segments add up for plant height; with only 24 double carriers, the test can't rule out a modest interaction. Script: `final_script/ra2_tsh4_interaction.py`; output `data/introgression_fix1/checks/ra2_tsh4.txt`.
- **Q3, chr3 13.2 Mb:** Zm00001eb123060 (*ra2*, *ramosa2*; 12.83 Mb, in a window with −log10 p 6.3; the top window is 0.4 Mb away). *ra2* encodes a LATERAL ORGAN BOUNDARY (LOB) domain transcription factor. It is one of the three *ramosa* genes that keep the short branches of the maize inflorescence (the spikelet pairs) short: *Fine-mapping (3.4.3): PIP 0.057 (highest window 0.09); in the top windows (12.8–13.5 Mb, 19 genes; 3 genes in the first model).*
  - In *ra2* mutants, spikelet pairs are replaced by longer branches that carry single spikelets or groups of several spikelets, the tassel has more secondary branches, the branches stand more upright, and the ears branch. *ra2* is switched on briefly in the small group of cells that will become each new branch meristem of the tassel and the ear (long branches, spikelet pairs and spikelets), but not in the vegetative shoot meristem. It acts upstream of the other branching genes and is needed for normal expression of *ra1* (Bortiri et al. 2006).
  - Xie et al. (2023) showed how RA2 limits tassel branching: it binds the promoter of *tassel sheath4* (*TSH4*, Zm00001eb316740), an SBP-box gene that promotes branching, and represses it. CRISPR *ra2* knockouts had about four times as many tassel branches as the wild type, and plants overexpressing RA2 about half as many. The circadian clock protein ZmELF3.1 (Zm00001eb159890) binds both RA2 and ZmELF4.2 (Zm00001eb311850) and makes RA2 bind and repress *TSH4* more strongly; *Zmelf3.1* mutants flower later and have 28–57% more tassel branches. This links the flowering clock to tassel branching.

  Neither paper reports an effect of *ra2* on plant height; the mutants were described for their inflorescences. In our lines the teosinte allele at Q3 makes plants 3.73 cm shorter per copy, alike across taxa, and the same window lowers the ear (EH Q3, −2.15 cm per copy). We did not score tassel branching, so we cannot check whether this segment changes the tassel the way *ra2* mutants do. *ra2* is the best-known gene in the interval (12.4–18.7 Mb), not a proven cause of the height effect.
- **Q4, chr4 198.4 Mb:** Zm00001eb197430 (*ns2*, *narrow sheath2*; Zm00001d052598 in B73 v4; 196.23 Mb, in a window with −log10 p 6.8). *ns2* encodes a WUSCHEL-related homeobox protein (WOX3); its *Arabidopsis* ortholog is *PRESSED FLOWER*, which is needed for the lateral sepals to form. *ns2* and its duplicate *ns1* together give the classic *narrow sheath* mutant (Nardmann et al. 2004). *Fine-mapping (3.4.3): PIP 0.027 (highest window 0.033); in the top windows (195.5–198.7 Mb, 40 genes).*
  - Wang X et al. (2023) took 12 tassel branch number (TBN) QTL mapped in 10 RIL populations from 14 elite inbreds, and ranked the genes in each interval with QTG-Miner: expression differences in young tassels between pairs of RILs that differ only at the QTL, sequence variants that change the protein, and a multi-omics gene network. They then tested the top genes with EMS or CRISPR mutants and confirmed seven genes for seven QTL.
  - *ns2* is the gene of *qTBN4-2*: *ns2* mutants had significantly more tassel branches than their wild-type siblings. The authors note that *ns2* mutants are short (brachytic), with narrow sheaths and shorter internodes, and propose that NS2, like *lrs1* (the *lg2* paralog behind *qTBN8-2*), affects tassel branching by marking the boundaries between organs.

  In our lines the teosinte allele at Q4 makes plants taller (+4.05 cm per copy; joint −log10 p 2.8; +3.7 to +11.9 cm in the five taxa). The same segment raises the ear (EH Q4, +3.15 cm, the strongest EH QTL), shortens the leaf blade (BL Q5, −1.20 cm) and lengthens the sheath (SL Q3). Loss of *ns2* shortens internodes, so a teosinte allele with more NS2 activity would fit taller plants, but we have no expression data to test this. *ns2* is 2.2 Mb from the PH top window, and the interval (194.5–201.4 Mb) holds 164 genes.

The chr6 peaks (Q9–Q11), Q5 chr4 218.6 Mb and Q7 chr5 64.0 Mb (*td1*, Zm00001eb228140) lose their signal in the joint model. The chr5 2.6 Mb QTL (Q6), shared with EH (Q5) and SL (Q4), now holds, and its probability is concentrated: its best window (2.6–2.7 Mb) has PIP 0.53.

<details><summary>All PH QTL with per-taxon effects and all genes in the interval</summary>

| QTL | Position | Interval (Mb) | Scan −log10 p | Effect per copy (joint) | Joint −log10 p | Taxa differ −log10 p | Per-taxon effects (carrier plots) | Genes in interval | Top windows, Mb (genes) | Best window PIP |
|---|---|---|---|---|---|---|---|---|---|---|
| Q1 | chr1 184.1 Mb | 183.6–186.1 | 5.3 | -3.83 cm ± 0.86 | **5.1** | 0.4 | Zd -3.37 (18), Zh -0.07 (6), Zl +0.24 (9), Zv -4.46 (21), Zx -4.17 (72) | [63 genes](data/introgression_fix1/qtl_genes/PH_qtl_genes.tsv) | 184.0–185.7 (30) | 0.07 |
| Q2 | chr2 2.2 Mb | 0.0–3.4 | 4.9 | -2.64 cm ± 0.59 | **5.1** | 0.2 | Zd -4.00 (43), Zh -2.98 (9), Zl -1.91 (15), Zv -2.31 (63), Zx -2.21 (111) | [165 genes](data/introgression_fix1/qtl_genes/PH_qtl_genes.tsv) | 1.8–2.7 (45) | 0.14 |
| Q3 | chr3 13.2 Mb | 12.4–18.7 | 6.5 | -3.73 cm ± 0.71 | **6.8** | 0.7 | Zd -4.92 (28), Zh +5.69 (14), Zl -0.23 (10), Zv -3.63 (27), Zx -4.42 (87) | [169 genes](data/introgression_fix1/qtl_genes/PH_qtl_genes.tsv) | 12.8–13.5 (19) | 0.09 |
| Q4 | chr4 198.4 Mb | 194.5–201.4 | 7.0 | +4.05 cm ± 1.28 | **2.8** | 0.1 | Zd +6.75 (16), Zh +4.71 (7), Zl +11.94 (6), Zv +5.04 (33), Zx +3.72 (75) | [164 genes](data/introgression_fix1/qtl_genes/PH_qtl_genes.tsv) | 195.5–198.7 (40) | 0.03 |
| Q5 | chr4 218.6 Mb | 215.2–223.7 | 4.9 | +0.23 cm ± 1.38 | 0.1 | 0.1 | Zd -4.68 (15), Zl +11.94 (6), Zv -1.92 (32), Zx +0.68 (69) | [142 genes](data/introgression_fix1/qtl_genes/PH_qtl_genes.tsv) | – | 0.00 |
| Q6 | chr5 2.6 Mb | 2.4–3.4 | 4.2 | +3.05 cm ± 0.74 | **4.5** | 0.5 | Zd +4.69 (32), Zh +3.90 (5), Zl +0.75 (11), Zv +1.75 (44), Zx +3.07 (73) | [81 genes](data/introgression_fix1/qtl_genes/PH_qtl_genes.tsv) | 2.6–2.7 (11) | 0.53 |
| Q7 | chr5 64.0 Mb | 60.0–69.1 | 4.7 | -2.24 cm ± 1.13 | 1.3 | 0.1 | Zd -8.03 (8), Zl -4.07 (12), Zv -1.93 (37), Zx -3.66 (51) | [217 genes](data/introgression_fix1/qtl_genes/PH_qtl_genes.tsv) | 60.8–68.5 (113) | 0.02 |
| Q8 | chr5 118.0 Mb | 102.1–121.4 | 6.2 | -3.46 cm ± 1.15 | **2.6** | 0.0 | Zd +3.98 (8), Zl -3.08 (8), Zv -3.82 (28), Zx -2.53 (50) | [115 genes](data/introgression_fix1/qtl_genes/PH_qtl_genes.tsv) | 113.0–120.7 (44) | 0.02 |
| Q9 | chr6 50.6 Mb | 49.4–51.3 | 6.4 | +5.94 cm ± 2.15 | 2.3 | 0.2 | Zd +2.06 (14), Zh +6.06 (6), Zl +3.00 (10), Zv -1.84 (29), Zx +8.40 (41) | [17 genes](data/introgression_fix1/qtl_genes/PH_qtl_genes.tsv) | 49.4–51.2 (5) | 0.12 |
| Q10 | chr6 67.2 Mb | 64.2–68.7 | 5.0 | -3.16 cm ± 2.12 | 0.9 | 0.3 | Zd +1.86 (13), Zh +4.59 (6), Zl +4.14 (10), Zv -1.01 (30), Zx -4.34 (50) | [65 genes](data/introgression_fix1/qtl_genes/PH_qtl_genes.tsv) | 64.2–68.2 (14) | 0.10 |
| Q11 | chr6 105.0 Mb | 102.9–106.2 | 6.5 | +3.26 cm ± 1.14 | 2.4 | 2.3 | Zd +3.34 (10), Zh -2.93 (7), Zl -1.10 (8), Zv +9.56 (28), Zx +2.04 (45) | [84 genes](data/introgression_fix1/qtl_genes/PH_qtl_genes.tsv) | 104.6–105.6 (27) | 0.08 |
| Q12 | chr7 137.8 Mb | 136.3–140.5 | 10.5 | -5.10 cm ± 0.78 | **10.2** | 0.3 | Zd -5.01 (12), Zh -5.94 (8), Zl +0.38 (9), Zv -4.58 (45), Zx -5.78 (69) | [103 genes](data/introgression_fix1/qtl_genes/PH_qtl_genes.tsv) | 137.5–139.1 (26) | 0.08 |
| Q13 | chr9 139.8 Mb | 139.0–141.0 | 10.1 | +4.64 cm ± 0.70 | **10.5** | 0.5 | Zd +3.88 (26), Zl +2.08 (25), Zv +6.23 (41), Zx +4.19 (87) | [61 genes](data/introgression_fix1/qtl_genes/PH_qtl_genes.tsv) | 139.6–139.9 (5) | 0.21 |
| Q14 (cond.) | chr10 6.5 Mb | 5.0–15.2 | 8.3 | +3.27 cm ± 0.82 | **4.2** | 2.4 | Zd +7.74 (25), Zh -13.61 (5), Zl +5.96 (11), Zv +2.98 (34), Zx +2.64 (60) | [234 genes](data/introgression_fix1/qtl_genes/PH_qtl_genes.tsv) | – | 0.07 |
| Q15 | chr10 99.5 Mb | 93.0–99.9 | 11.0 | +4.70 cm ± 0.94 | **6.3** | 0.1 | Zd +2.46 (24), Zv +7.09 (23), Zx +4.06 (66) | [111 genes](data/introgression_fix1/qtl_genes/PH_qtl_genes.tsv) | 97.9–99.6 (33) | 0.04 |

Bonferroni over the 15 QTL: −log10 p > 2.48 (bold). Top windows: the fewest windows holding half of the QTL's fine-mapping probability (section 3.4.3). Every gene in each interval, with its GO annotation, the k-mer GWAS −log10 p of its window and its window's PIP: [`data/introgression_fix1/qtl_genes/PH_qtl_genes.tsv`](data/introgression_fix1/qtl_genes/PH_qtl_genes.tsv).

</details>

#### 3.6.7 Ear height (EH)

![EH: k-mer ancestry GWAS and joint QTL model](figures/fix1/joint_qtl_EH.png)

16 regions, 9 hold (Bonferroni over 16: −log10 p > 2.51). Ear height follows flowering: several QTL hold known flowering genes in their intervals. Six of the first model's seven held QTL are found again (chr2 25 Mb, chr3 13 Mb, chr4 197 Mb, chr5 2.6 Mb, chr7 184 Mb, chr8 126–129 Mb); chr7 124 Mb (Q8) no longer holds, and chr2 83.5 Mb, chr6 105 Mb and chr10 99.5 Mb (*ZmEMF1a*, just below the cut before) now hold.

| QTL | Known genes | Fine-mapping PIP | Effect (joint) | Joint −log10 p | Taxa differ | What it says |
|---|---|---|---|---|---|---|
| **Q4 chr4 196.6 Mb** | Zm00001eb197430 (ns2) | 0.026 (**top**) | +3.15 cm | 8.1 | 0.1 (no) | No significant difference between taxa (Zd +2.88, Zh +0.18, Zl +6.76, Zv +2.49, Zx +3.62). The strongest EH QTL; shared with PH Q4, BL Q5 and SL Q3. |
| **Q11 chr8 129.2 Mb** | Zm00001eb353250 (ZCN8) | 0.016 (**top**) | +2.82 cm | 6.0 | 0.3 (no) | No significant difference between taxa (Zd +2.22, Zh +5.38, Zl +1.48, Zv +2.98, Zx +2.41). *ZCN8* is the maize florigen (126.68 Mb, inside the interval 125.7–133.7 Mb). |
| **Q10 chr7 184.4 Mb** | Zm00001eb331630 (PRR37-like); Zm00001eb331880 (dlf1) | 0.018 (outside CS); 0.076 (**top**) | +2.15 cm | 5.7 | 0.4 (no) | No significant difference between taxa (Zd +2.95, Zl −0.62, Zv +1.63, Zx +2.13). Same region as BL Q7. |
| **Q3 chr3 13.2 Mb** | Zm00001eb123060 (ra2) | 0.047 (**top**) | −2.15 cm | 5.1 | 1.8 (no) | No significant difference between taxa (Zd −5.33, Zh +0.25, Zl +0.11, Zv −3.38, Zx −1.22). Same window as PH Q3. |
| **Q1 chr2 25.4 Mb** | — | 0.024 (best window) | +2.06 cm | 5.0 | 1.3 (no) | Higher ears (Zh +2.15, Zv +2.18, Zx +2.73; Zd −1.09). *ckx12* and *GA2ox2* (29.1–29.3 Mb), in this QTL's interval in the first model, are now just outside it (16.2–28.7 Mb). |
| **Q5 chr5 2.6 Mb** | — | 0.244 (best window) | +2.06 cm | 4.5 | 0.6 (no) | Higher ears in every taxon. Same region as PH Q6 and SL Q4. |
| **Q2 chr2 83.5 Mb** | — | 0.010 (best window) | +1.55 cm | 2.7 | 0.6 (no) | Wide interval (79.6–98.4 Mb). |
| **Q7 chr6 105.0 Mb** | — | 0.078 (best window) | +2.26 cm | 2.6 | 3.2 (yes) | The only QTL where taxa differ: *parviglumis* +6.64 cm and *diploperennis* +3.61 cm, *luxurians* −2.81 cm. |
| **Q16 chr10 99.5 Mb** | Zm00001eb417870 (ZmEMF1a); Zm00001eb418700 (ZmCCT10) | 0.006 (CS); 0.008 (CS) | +2.41 cm | 2.6 | 0.0 (no) | Holds just above the cut (Zd −2.13, Zv +3.19, Zx +3.29). |

Not held: Q12 chr9 116.8 Mb (+1.57 cm, joint −log10 p 2.1; *ZmCCT9*, Zm00001eb391230; Huang et al. 2018) and Q8 chr7 124.0 Mb (−1.64 cm, 1.8).

**Candidate gene in Q16 (chr10 99.5 Mb):** Zm00001eb417870 (*ZmEMF1a*; Zm00001d024813 in B73 v4; 90.88 Mb, window −log10 p 6.6; interval 88.8–99.6 Mb). It is an EMBRYONIC FLOWER1 (EMF1)-like Polycomb protein; EMF1 is named after the *Arabidopsis* mutant that flowers right after germination. Li et al. (2026) cloned it as the gene behind the maize small-kernel mutant *sks1* (three EMS alleles): *Fine-mapping (3.4.3): PIP 0.006 (highest window 0.03); in the 95% credible set, but outside the top windows (93.3–99.6 Mb, 52 genes).*
- ZmEMF1a binds subunits of both Polycomb complexes (PRC1 and PRC2) and places the repressive histone marks H3K27me3 (at 6 days after pollination) and H2Aub1 (at 10 days) on its targets in the endosperm. Without it, endosperm cells keep dividing and fail to expand and differentiate, so the kernels are small.
- Its main targets are MADS-box genes: 25 MADS genes lose H3K27me3 and are switched on in the mutant.
- **The link to ear height:** besides small kernels, *sks1* plants are significantly shorter, with **lower ears**, shorter internodes and a shorter tassel spindle.
- Its paralog *ZmEMF1b* (Zm00001d010894, on chr8) is expressed mostly in tassels.

In our lines the teosinte segment raises the ear by 2.41 cm per copy at Q16 (joint −log10 p 2.6, just above the cut of 2.51; in the first model 2.5, just below), and every window from 89 to 100 Mb is significant for EH (−log10 p 6.0–7.2), so genes in the interval can't be told apart with our tract lengths. The interval also holds *ZmCCT10* (Zm00001eb418700, 96.18 Mb), the flowering gene behind DTS Q2, DTA Q6 and PH Q15, so later flowering may add to the ear-height effect; and two less-known genes: another EMF1-like protein (Zm00001eb418530, 95.20 Mb) and *mads85* (Zm00001eb417540, 89.03 Mb), related to rice *OsMADS87*, which controls seed size and heat sensitivity in young endosperm (Chen C et al. 2016).

<details><summary>All EH QTL with per-taxon effects and all genes in the interval</summary>

| QTL | Position | Interval (Mb) | Scan −log10 p | Effect per copy (joint) | Joint −log10 p | Taxa differ −log10 p | Per-taxon effects (carrier plots) | Genes in interval | Top windows, Mb (genes) | Best window PIP |
|---|---|---|---|---|---|---|---|---|---|---|
| Q1 | chr2 25.4 Mb | 16.2–28.7 | 6.5 | +2.06 cm ± 0.47 | **5.0** | 1.3 | Zd -1.09 (31), Zh +2.15 (5), Zl +0.64 (21), Zv +2.18 (48), Zx +2.73 (107) | [373 genes](data/introgression_fix1/qtl_genes/EH_qtl_genes.tsv) | 18.1–26.3 (86) | 0.02 |
| Q2 | chr2 83.5 Mb | 79.6–98.4 | 5.4 | +1.55 cm ± 0.51 | **2.7** | 0.6 | Zd +0.61 (31), Zh +3.86 (7), Zl +3.66 (14), Zv +3.37 (41), Zx +0.83 (101) | [145 genes](data/introgression_fix1/qtl_genes/EH_qtl_genes.tsv) | 81.3–95.6 (52) | 0.01 |
| Q3 | chr3 13.2 Mb | 12.4–19.7 | 5.0 | -2.15 cm ± 0.48 | **5.1** | 1.8 | Zd -5.33 (28), Zh +0.25 (14), Zl +0.11 (10), Zv -3.38 (27), Zx -1.22 (87) | [182 genes](data/introgression_fix1/qtl_genes/EH_qtl_genes.tsv) | 12.6–17.0 (32) | 0.07 |
| Q4 | chr4 196.6 Mb | 194.7–201.4 | 7.2 | +3.15 cm ± 0.55 | **8.1** | 0.1 | Zd +2.88 (16), Zh +0.18 (7), Zl +6.76 (6), Zv +2.49 (38), Zx +3.62 (75) | [161 genes](data/introgression_fix1/qtl_genes/EH_qtl_genes.tsv) | 195.5–198.7 (25) | 0.07 |
| Q5 | chr5 2.6 Mb | 1.0–2.9 | 4.5 | +2.06 cm ± 0.50 | **4.5** | 0.6 | Zd +3.45 (32), Zh +3.28 (5), Zl +3.77 (11), Zv +1.09 (44), Zx +1.62 (73) | [146 genes](data/introgression_fix1/qtl_genes/EH_qtl_genes.tsv) | 2.4–2.7 (34) | 0.24 |
| Q6 | chr6 50.6 Mb | 48.6–51.9 | 4.1 | +1.17 cm ± 0.70 | 1.0 | 0.3 | Zd -0.05 (14), Zh +0.69 (6), Zl +0.89 (10), Zv -1.14 (29), Zx +1.70 (41) | [27 genes](data/introgression_fix1/qtl_genes/EH_qtl_genes.tsv) | 49.3–51.3 (17) | 0.04 |
| Q7 | chr6 105.0 Mb | 102.9–106.5 | 5.5 | +2.26 cm ± 0.75 | **2.6** | **3.2** | Zd +3.61 (10), Zh -0.10 (7), Zl -2.81 (8), Zv +6.64 (28), Zx +0.91 (45) | [91 genes](data/introgression_fix1/qtl_genes/EH_qtl_genes.tsv) | 104.6–105.6 (27) | 0.08 |
| Q8 | chr7 124.0 Mb | 122.6–131.5 | 4.6 | -1.64 cm ± 0.68 | 1.8 | 0.5 | Zd +0.30 (11), Zh -5.22 (6), Zl -2.44 (7), Zv -1.08 (45), Zx -1.85 (58) | [169 genes](data/introgression_fix1/qtl_genes/EH_qtl_genes.tsv) | 122.9–126.0 (48) | 0.03 |
| Q9 | chr7 136.3 Mb | 135.2–140.5 | 4.2 | -1.11 cm ± 0.64 | 1.1 | 0.4 | Zd -2.70 (13), Zh -1.73 (8), Zl +3.20 (9), Zv -0.20 (46), Zx -1.77 (65) | [137 genes](data/introgression_fix1/qtl_genes/EH_qtl_genes.tsv) | 136.3–139.4 (73) | 0.03 |
| Q10 | chr7 184.4 Mb | 183.4–185.8 | 6.0 | +2.15 cm ± 0.45 | **5.7** | 0.4 | Zd +2.95 (32), Zl -0.62 (9), Zv +1.63 (48), Zx +2.13 (95) | [82 genes](data/introgression_fix1/qtl_genes/EH_qtl_genes.tsv) | 184.2–185.2 (26) | 0.08 |
| Q11 | chr8 129.2 Mb | 125.7–133.7 | 6.3 | +2.82 cm ± 0.58 | **6.0** | 0.3 | Zd +2.22 (23), Zh +5.38 (11), Zl +1.48 (8), Zv +2.98 (30), Zx +2.41 (62) | [179 genes](data/introgression_fix1/qtl_genes/EH_qtl_genes.tsv) | 126.2–130.1 (43) | 0.05 |
| Q12 | chr9 116.8 Mb | 115.2–122.2 | 5.3 | +1.57 cm ± 0.60 | 2.1 | 1.1 | Zd +3.13 (22), Zl -3.26 (21), Zv +1.61 (47), Zx +2.11 (89) | [166 genes](data/introgression_fix1/qtl_genes/EH_qtl_genes.tsv) | 116.5–120.2 (42) | 0.03 |
| Q13 | chr9 139.8 Mb | 127.8–146.1 | 4.6 | +1.19 cm ± 0.58 | 1.4 | 0.1 | Zd +1.66 (26), Zl +2.77 (24), Zv +0.45 (41), Zx +0.99 (87) | [493 genes](data/introgression_fix1/qtl_genes/EH_qtl_genes.tsv) | – | 0.01 |
| Q14 | chr10 6.8 Mb | 5.0–9.5 | 5.4 | +1.79 cm ± 0.65 | 2.3 | 1.4 | Zd +3.79 (23), Zh -4.75 (5), Zl +4.76 (11), Zv +0.99 (31), Zx +1.12 (63) | [122 genes](data/introgression_fix1/qtl_genes/EH_qtl_genes.tsv) | 5.0–7.9 (56) | 0.05 |
| Q15 | chr10 66.8 Mb | 63.0–67.8 | 5.0 | +0.00 cm ± 0.93 | 0.0 | 0.2 | Zd +2.76 (24), Zv +0.62 (22), Zx -0.47 (55) | [49 genes](data/introgression_fix1/qtl_genes/EH_qtl_genes.tsv) | – | 0.00 |
| Q16 | chr10 99.5 Mb | 88.8–99.6 | 7.2 | +2.41 cm ± 0.80 | **2.6** | 0.0 | Zd -2.13 (24), Zv +3.19 (23), Zx +3.29 (66) | [187 genes](data/introgression_fix1/qtl_genes/EH_qtl_genes.tsv) | 93.3–99.6 (52) | 0.03 |

Bonferroni over the 16 QTL: −log10 p > 2.51 (bold). Top windows: the fewest windows holding half of the QTL's fine-mapping probability (section 3.4.3). Every gene in each interval, with its GO annotation, the k-mer GWAS −log10 p of its window and its window's PIP: [`data/introgression_fix1/qtl_genes/EH_qtl_genes.tsv`](data/introgression_fix1/qtl_genes/EH_qtl_genes.tsv).

</details>

#### 3.6.8 Ear number (EN; number of nodes with ears)

![EN: k-mer ancestry GWAS and joint QTL model](figures/fix1/joint_qtl_EN.png)

17 regions plus one QTL added by the conditional scans (Q9 chr4 183.6 Mb, open triangle); 12 of the 18 hold (Bonferroni over 18: −log10 p > 2.56). EN counts the stem nodes that bear an ear. **Every QTL that holds is positive: teosinte alleles add ear-bearing nodes**, as expected from teosinte, which bears ears at many nodes.

The scan is almost the same as in the first model, but chr2, chr4 and chr7 are significant over long stretches, and within such a stretch the lead window can move between two nearly equal peaks. So several QTL sit at a different place than before: chr2 19.0 Mb (Q2; −log10 p 18.3, against 17.9 at the first model's 22.4 Mb), chr2 205.3 Mb (Q3; first model 202.4 Mb), chr4 119.3 Mb and 183.6 Mb (Q8, Q9; first model 120.2 and 167.0 Mb, which still has −log10 p 14.6 in the scan) and chr7 165.3 Mb (Q12; first model 143.6 Mb, still 10.1 in the scan). The QTL at chr1 22.9 Mb (Q1), chr6 96 and 173 Mb (Q10, Q11), chr8 36.7 Mb (Q13) and chr10 148.6 Mb (Q18; first model 144.0 Mb) are the same or close.

| QTL | Known genes | Fine-mapping PIP | Effect (joint) | Joint −log10 p | Taxa differ | What it says |
|---|---|---|---|---|---|---|
| **Q1 chr1 22.9 Mb** | Zm00001eb007670 (ZmD14) | 0.012 (outside CS) | +0.20 | 16.3 | 2.1 (no) | No significant difference between taxa (Zd +0.21, Zh +0.15, Zl −0.12, Zv +0.24, Zx +0.22). *ZmD14* is the strigolactone receptor (shoot branching). |
| **Q2 chr2 19.0 Mb** | — | 0.14 (best window) | +0.16 | 15.5 | 0.5 (no) | The strongest EN scan peak (18.3); alike in the taxa (Zd +0.07 to Zh, Zx +0.19). |
| **Q3 chr2 205.3 Mb** | — | 0.23 (best window) | +0.16 | 15.1 | 3.5 (yes) | Differs between taxa: *parviglumis* +0.22 and *mexicana* +0.19, *diploperennis* −0.09, *luxurians* +0.02. Same window as Prolif Q6. |
| **Q12 chr7 165.3 Mb** | — | 0.19 (best window) | +0.15 | 11.2 | 0.5 (no) | Alike in the taxa (Zd +0.11, Zv +0.16, Zx +0.16). Next to the ASI Q10 peak (165.6 Mb). |

Not holding: Q4 chr3 12.4 Mb (+0.07, joint 1.4) holds *ra2* (Zm00001eb123060). *Rap2.7* (Zm00001eb355240) and *lg4a* (Zm00001eb354880), in a chr8 QTL of the first model that did not hold, and the clock-coupling gene *ZGT* (Zm00001eb020350; first model's chr1 81.3 Mb QTL) are not in any QTL interval in this model.

<details><summary>All EN QTL with per-taxon effects and all genes in the interval</summary>

| QTL | Position | Interval (Mb) | Scan −log10 p | Effect per copy (joint) | Joint −log10 p | Taxa differ −log10 p | Per-taxon effects (carrier plots) | Genes in interval | Top windows, Mb (genes) | Best window PIP |
|---|---|---|---|---|---|---|---|---|---|---|
| Q1 | chr1 22.9 Mb | 22.6–23.7 | 15.6 | +0.20 ± 0.02 | **16.3** | 2.1 | Zd +0.21 (15), Zh +0.15 (10), Zl -0.12 (13), Zv +0.24 (30), Zx +0.22 (58) | [34 genes](data/introgression_fix1/qtl_genes/EN_qtl_genes.tsv) | 22.9–23.1 (8) | 0.30 |
| Q2 | chr2 18.9 Mb | 18.1–19.7 | 18.3 | +0.16 ± 0.02 | **15.5** | 0.5 | Zd +0.07 (24), Zh +0.19 (5), Zl +0.16 (19), Zv +0.13 (46), Zx +0.19 (108) | [54 genes](data/introgression_fix1/qtl_genes/EN_qtl_genes.tsv) | 18.5–19.0 (15) | 0.14 |
| Q3 | chr2 205.2 Mb | 204.4–205.5 | 16.2 | +0.16 ± 0.02 | **15.1** | **3.5** | Zd -0.09 (26), Zl +0.02 (11), Zv +0.22 (50), Zx +0.19 (107) | [39 genes](data/introgression_fix1/qtl_genes/EN_qtl_genes.tsv) | 205.0–205.3 (16) | 0.23 |
| Q4 | chr3 12.4 Mb | 10.9–13.9 | 5.8 | +0.07 ± 0.03 | 1.4 | 0.0 | Zd +0.03 (29), Zh +0.09 (13), Zl +0.15 (10), Zv +0.13 (29), Zx +0.01 (92) | [73 genes](data/introgression_fix1/qtl_genes/EN_qtl_genes.tsv) | 11.0–13.2 (25) | 0.05 |
| Q5 | chr3 20.9 Mb | 17.5–22.8 | 5.2 | +0.02 ± 0.04 | 0.2 | 0.1 | Zd +0.10 (23), Zh -0.04 (12), Zl -0.09 (10), Zv +0.00 (25), Zx +0.02 (86) | [99 genes](data/introgression_fix1/qtl_genes/EN_qtl_genes.tsv) | – | 0.00 |
| Q6 | chr3 37.0 Mb | 27.9–43.4 | 4.6 | +0.02 ± 0.03 | 0.2 | 0.6 | Zd -0.06 (24), Zh -0.04 (11), Zl -0.10 (8), Zv -0.01 (25), Zx +0.08 (73) | [261 genes](data/introgression_fix1/qtl_genes/EN_qtl_genes.tsv) | – | 0.00 |
| Q7 | chr3 175.1 Mb | 171.4–176.8 | 5.1 | +0.07 ± 0.02 | **2.8** | 2.1 | Zd -0.06 (24), Zh +0.04 (10), Zl +0.10 (14), Zv +0.03 (36), Zx +0.13 (77) | [131 genes](data/introgression_fix1/qtl_genes/EN_qtl_genes.tsv) | 173.8–175.4 (42) | 0.05 |
| Q8 | chr4 119.2 Mb | 119.2–121.3 | 14.8 | +0.18 ± 0.03 | **8.7** | 0.6 | Zh -0.00 (6), Zv +0.21 (32), Zx +0.17 (50) | [12 genes](data/introgression_fix1/qtl_genes/EN_qtl_genes.tsv) | 119.2–119.5 (1) | 0.19 |
| Q9 (cond.) | chr4 183.6 Mb | 179.8–212.2 | 13.0 | +0.14 ± 0.02 | **7.6** | 0.1 | Zd +0.05 (14), Zh +0.00 (9), Zl +0.04 (6), Zv +0.12 (40), Zx +0.16 (76) | [852 genes](data/introgression_fix1/qtl_genes/EN_qtl_genes.tsv) | 180.3–212.2 (97) | 0.08 |
| Q10 | chr6 96.2 Mb | 94.2–99.2 | 5.2 | +0.11 ± 0.03 | **4.7** | 0.1 | Zd +0.10 (10), Zh +0.16 (7), Zl +0.12 (7), Zv +0.10 (30), Zx +0.12 (43) | [119 genes](data/introgression_fix1/qtl_genes/EN_qtl_genes.tsv) | 95.6–96.9 (25) | 0.06 |
| Q11 | chr6 173.2 Mb | 171.5–175.3 | 4.7 | +0.10 ± 0.02 | **4.4** | 0.6 | Zd +0.04 (16), Zl +0.10 (9), Zv +0.06 (35), Zx +0.15 (55) | [192 genes](data/introgression_fix1/qtl_genes/EN_qtl_genes.tsv) | 171.8–173.3 (42) | 0.06 |
| Q12 | chr7 165.2 Mb | 165.1–166.6 | 11.0 | +0.15 ± 0.02 | **11.2** | 0.5 | Zd +0.11 (17), Zl +0.01 (13), Zv +0.16 (44), Zx +0.16 (77) | [56 genes](data/introgression_fix1/qtl_genes/EN_qtl_genes.tsv) | 165.2–165.5 (7) | 0.19 |
| Q13 | chr8 36.6 Mb | 29.3–41.0 | 7.7 | +0.12 ± 0.03 | **3.7** | 1.2 | Zd +0.01 (28), Zh -0.00 (11), Zl -0.04 (6), Zv +0.20 (46), Zx +0.19 (60) | [140 genes](data/introgression_fix1/qtl_genes/EN_qtl_genes.tsv) | 32.6–37.1 (26) | 0.04 |
| Q14 | chr8 98.2 Mb | 95.0–106.9 | 5.1 | +0.04 ± 0.03 | 0.7 | 0.0 | Zd +0.03 (19), Zh +0.05 (8), Zl +0.10 (5), Zv -0.04 (38), Zx +0.01 (74) | [207 genes](data/introgression_fix1/qtl_genes/EN_qtl_genes.tsv) | – | 0.00 |
| Q15 | chr9 19.4 Mb | 18.2–20.1 | 7.5 | +0.10 ± 0.03 | **3.8** | 0.4 | Zd +0.08 (26), Zh +0.34 (8), Zl +0.08 (22), Zv +0.15 (51), Zx +0.06 (94) | [86 genes](data/introgression_fix1/qtl_genes/EN_qtl_genes.tsv) | 18.8–19.5 (21) | 0.13 |
| Q16 | chr9 70.0 Mb | 67.8–70.3 | 5.5 | -0.05 ± 0.07 | 0.3 | 0.5 | Zd +0.08 (25), Zh -0.47 (7), Zl +0.04 (21), Zv +0.01 (45), Zx +0.07 (90) | [29 genes](data/introgression_fix1/qtl_genes/EN_qtl_genes.tsv) | 67.8–70.3 (5) | 0.15 |
| Q17 | chr9 81.0 Mb | 78.8–82.9 | 5.8 | +0.07 ± 0.06 | 0.6 | 0.8 | Zd +0.07 (27), Zh +0.07 (6), Zl +0.13 (21), Zv +0.01 (51), Zx +0.06 (84) | [41 genes](data/introgression_fix1/qtl_genes/EN_qtl_genes.tsv) | – | 0.23 |
| Q18 | chr10 148.6 Mb | 148.2–149.3 | 5.5 | +0.10 ± 0.02 | **6.1** | 0.2 | Zd +0.10 (39), Zh +0.10 (7), Zl +0.05 (14), Zv +0.11 (40), Zx +0.09 (85) | [60 genes](data/introgression_fix1/qtl_genes/EN_qtl_genes.tsv) | 148.5–148.7 (7) | 0.34 |

Bonferroni over the 18 QTL: −log10 p > 2.56 (bold). Top windows: the fewest windows holding half of the QTL's fine-mapping probability (section 3.4.3). Every gene in each interval, with its GO annotation, the k-mer GWAS −log10 p of its window and its window's PIP: [`data/introgression_fix1/qtl_genes/EN_qtl_genes.tsv`](data/introgression_fix1/qtl_genes/EN_qtl_genes.tsv).

</details>

#### 3.6.9 Prolificacy (Prolif; most ears on one node)

![Prolif: k-mer ancestry GWAS and joint QTL model](figures/fix1/joint_qtl_Prolif.png)

13 regions plus two QTL added by the conditional scans (Q5 chr2 34.4 Mb and Q11 chr4 114.3 Mb, open triangles); 9 of the 15 hold (Bonferroni over 15: −log10 p > 2.48). Prolificacy is the largest number of ears on a single node (teosinte clusters several ears per node). **Largely the same loci as EN, all positive** (teosinte adds ears per node): the same loci raise both the number of ear-bearing nodes and the number of ears per node. Eight of the first model's ten held QTL are found again; chr7 142.3 Mb no longer holds (1.4), and the first model's chr2 191.3 Mb QTL is now at chr2 205.2 Mb (Q6), the same window as EN Q3 (chr2 is significant from 101 to 220 Mb, and the lead moved between two nearly equal peaks; 191.3 Mb still has −log10 p 11.3 in the scan).

| QTL | Known genes | Fine-mapping PIP | Effect (joint) | Joint −log10 p | Taxa differ | What it says |
|---|---|---|---|---|---|---|
| **Q1 chr1 23.2 Mb** | — | 0.57 (best window; top windows 23.2–23.3 Mb, 2 genes) | +0.20 | 16.2 | 1.0 (no) | The strongest QTL (scan 19.4); the same peak as EN Q1. Fine-mapping puts more than half of the probability on one window. |
| **Q7 chr3 11.4 Mb** | Zm00001eb122680 (CYP734A6) | 0.158 (**top**) | +0.15 | 12.0 | 0.4 (no) | No significant difference between taxa (Zd +0.18, Zh +0.12, Zl +0.08, Zv +0.23, Zx +0.12). |
| **Q4 chr2 9.2 Mb** | — | 0.07 (best window) | +0.08 | 5.4 | 0.7 (no) | No significant difference between taxa (Zd +0.11, Zh +0.13, Zl +0.12, Zv +0.03, Zx +0.06). Its interval (6.5–9.8 Mb) no longer reaches *ckx12* and *GA2ox2*. |
| **Q5 chr2 34.4 Mb** (added) | Zm00001eb076990 (ckx12); Zm00001eb077080 (GA2ox2) | 0.001 (outside CS); <0.001 (outside CS) | +0.08 | 5.0 | 0.7 (no) | No significant difference between taxa (Zd +0.06, Zl +0.01, Zv +0.13, Zx +0.10). Added by the conditional scan: a second QTL on the chr2 short arm; its 29.7-Mb interval has no credible set. |
| **Q3 chr1 214.1 Mb** | Zm00001eb040280 (CKX10) | 0.010 (CS) | +0.12 | 6.2 | 0.6 (no) | No significant difference between taxa (Zd +0.05, Zh −0.09, Zl +0.02, Zv +0.13, Zx +0.15). *CKX10* is a cytokinin oxidase. Stronger than in the first model (joint 3.1). |

The other held QTL are Q6 chr2 205.2 Mb (+0.12; taxa differ 3.5: mostly *mexicana*), Q11 chr4 114.3 Mb (+0.14, added by the conditional scan), Q12 chr4 171.3 Mb (+0.15) and Q15 chr9 32.0 Mb (+0.09).

<details><summary>All Prolif QTL with per-taxon effects and all genes in the interval</summary>

| QTL | Position | Interval (Mb) | Scan −log10 p | Effect per copy (joint) | Joint −log10 p | Taxa differ −log10 p | Per-taxon effects (carrier plots) | Genes in interval | Top windows, Mb (genes) | Best window PIP |
|---|---|---|---|---|---|---|---|---|---|---|
| Q1 | chr1 23.2 Mb | 22.6–23.3 | 19.4 | +0.20 ± 0.02 | **16.2** | 1.0 | Zd +0.14 (15), Zh +0.07 (10), Zl -0.07 (13), Zv +0.26 (28), Zx +0.22 (57) | [24 genes](data/introgression_fix1/qtl_genes/Prolif_qtl_genes.tsv) | 23.2–23.3 (2) | 0.57 |
| Q2 | chr1 51.4 Mb | 48.5–54.0 | 4.8 | +0.03 ± 0.02 | 0.7 | 1.8 | Zd +0.05 (25), Zh +0.12 (11), Zl +0.08 (14), Zv -0.02 (24), Zx -0.03 (64) | [118 genes](data/introgression_fix1/qtl_genes/Prolif_qtl_genes.tsv) | – | 0.00 |
| Q3 | chr1 214.1 Mb | 211.4–216.5 | 5.6 | +0.12 ± 0.02 | **6.2** | 0.6 | Zd +0.05 (7), Zh -0.09 (7), Zl +0.02 (6), Zv +0.13 (31), Zx +0.15 (73) | [114 genes](data/introgression_fix1/qtl_genes/Prolif_qtl_genes.tsv) | 213.9–216.1 (38) | 0.08 |
| Q4 | chr2 9.2 Mb | 6.5–9.8 | 7.0 | +0.08 ± 0.02 | **5.4** | 0.7 | Zd +0.11 (37), Zh +0.13 (10), Zl +0.12 (17), Zv +0.03 (52), Zx +0.06 (98) | [151 genes](data/introgression_fix1/qtl_genes/Prolif_qtl_genes.tsv) | 6.5–9.4 (39) | 0.07 |
| Q5 (cond.) | chr2 34.4 Mb | 22.9–52.6 | 9.1 | +0.08 ± 0.02 | **5.0** | 0.7 | Zd +0.06 (33), Zl +0.01 (18), Zv +0.13 (46), Zx +0.10 (93) | [754 genes](data/introgression_fix1/qtl_genes/Prolif_qtl_genes.tsv) | – | 0.04 |
| Q6 | chr2 205.2 Mb | 202.3–205.8 | 12.6 | +0.12 ± 0.02 | **10.7** | **3.5** | Zd -0.03 (26), Zl +0.01 (12), Zv +0.08 (50), Zx +0.18 (107) | [135 genes](data/introgression_fix1/qtl_genes/Prolif_qtl_genes.tsv) | 205.0–205.4 (20) | 0.15 |
| Q7 | chr3 11.4 Mb | 11.0–11.7 | 17.2 | +0.15 ± 0.02 | **12.0** | 0.4 | Zd +0.18 (29), Zh +0.12 (13), Zl +0.08 (10), Zv +0.23 (31), Zx +0.12 (90) | [11 genes](data/introgression_fix1/qtl_genes/Prolif_qtl_genes.tsv) | 11.2–11.5 (4) | 0.26 |
| Q8 | chr3 73.5 Mb | 63.7–81.3 | 5.4 | -0.08 ± 0.06 | 0.7 | 0.6 | Zd +0.00 (25), Zh -0.04 (13), Zl -0.09 (8), Zv -0.07 (20), Zx +0.02 (65) | [107 genes](data/introgression_fix1/qtl_genes/Prolif_qtl_genes.tsv) | – | 0.03 |
| Q9 | chr3 119.8 Mb | 113.5–137.3 | 6.3 | +0.08 ± 0.06 | 0.9 | 0.8 | Zd +0.07 (26), Zh -0.04 (13), Zl -0.09 (8), Zv +0.01 (26), Zx +0.03 (74) | [341 genes](data/introgression_fix1/qtl_genes/Prolif_qtl_genes.tsv) | – | 0.00 |
| Q10 | chr3 172.3 Mb | 170.1–181.1 | 5.3 | +0.06 ± 0.02 | 2.4 | 1.5 | Zd +0.03 (24), Zh +0.04 (12), Zl +0.03 (12), Zv +0.00 (35), Zx +0.11 (83) | [278 genes](data/introgression_fix1/qtl_genes/Prolif_qtl_genes.tsv) | 171.7–180.2 (96) | 0.02 |
| Q11 (cond.) | chr4 114.3 Mb | 77.3–135.6 | 12.2 | +0.14 ± 0.03 | **5.2** | 0.1 | Zh +0.04 (6), Zv +0.16 (32), Zx +0.15 (50) | [503 genes](data/introgression_fix1/qtl_genes/Prolif_qtl_genes.tsv) | 100.6–121.7 (56) | 0.02 |
| Q12 | chr4 171.3 Mb | 168.5–172.1 | 16.5 | +0.15 ± 0.03 | **8.9** | 1.6 | Zd +0.23 (11), Zh -0.05 (5), Zl +0.08 (6), Zv +0.14 (28), Zx +0.17 (67) | [106 genes](data/introgression_fix1/qtl_genes/Prolif_qtl_genes.tsv) | 171.1–171.9 (13) | 0.07 |
| Q13 | chr7 136.2 Mb | 135.6–139.1 | 4.9 | +0.05 ± 0.03 | 1.3 | **4.2** | Zd -0.15 (14), Zh +0.36 (9), Zl -0.04 (8), Zv +0.08 (46), Zx +0.05 (64) | [98 genes](data/introgression_fix1/qtl_genes/Prolif_qtl_genes.tsv) | – | 0.03 |
| Q14 | chr7 142.2 Mb | 141.8–146.0 | 5.4 | +0.06 ± 0.03 | 1.4 | 1.3 | Zd +0.09 (16), Zh -0.05 (8), Zl +0.04 (10), Zv +0.02 (42), Zx +0.05 (73) | [111 genes](data/introgression_fix1/qtl_genes/Prolif_qtl_genes.tsv) | 142.2–145.4 (50) | 0.04 |
| Q15 | chr9 32.0 Mb | 25.4–32.7 | 5.0 | +0.09 ± 0.02 | **5.4** | 0.1 | Zd +0.11 (24), Zh -0.01 (8), Zl +0.14 (22), Zv +0.07 (46), Zx +0.09 (93) | [156 genes](data/introgression_fix1/qtl_genes/Prolif_qtl_genes.tsv) | 31.6–32.5 (10) | 0.07 |

Bonferroni over the 15 QTL: −log10 p > 2.48 (bold). Top windows: the fewest windows holding half of the QTL's fine-mapping probability (section 3.4.3). Every gene in each interval, with its GO annotation, the k-mer GWAS −log10 p of its window and its window's PIP: [`data/introgression_fix1/qtl_genes/Prolif_qtl_genes.tsv`](data/introgression_fix1/qtl_genes/Prolif_qtl_genes.tsv).

</details>

#### 3.6.10 Sheath length (SL; 2023 only)

![SL: k-mer ancestry GWAS and joint QTL model](figures/fix1/joint_qtl_SL.png)

13 regions, 11 hold (Bonferroni over 13: −log10 p > 2.41; 12 in the first model, where chr5 92 Mb and chr10 151 Mb held; chr5 161.8 Mb is new here). Sheath length shares several regions with plant height, as expected: the sheath wraps the internode, and longer internodes make both taller plants and longer sheaths.

<details><summary>All SL QTL with per-taxon effects and all genes in the interval</summary>

| QTL | Position | Interval (Mb) | Scan −log10 p | Effect per copy (joint) | Joint −log10 p | Taxa differ −log10 p | Per-taxon effects (carrier plots) | Genes in interval | Top windows, Mb (genes) | Best window PIP |
|---|---|---|---|---|---|---|---|---|---|---|
| Q1 | chr1 69.7 Mb | 67.6–71.1 | 5.2 | +0.29 ± 0.06 | **5.7** | 0.1 | Zd +0.25 (18), Zh +0.58 (10), Zl -0.25 (7), Zv +0.40 (15), Zx +0.25 (49) | [74 genes](data/introgression_fix1/qtl_genes/SL_qtl_genes.tsv) | 69.6–70.4 (10) | 0.12 |
| Q2 | chr2 6.0 Mb | 5.8–7.4 | 7.0 | -0.23 ± 0.04 | **7.2** | 1.3 | Zd -0.43 (41), Zh -0.17 (11), Zl -0.66 (13), Zv -0.16 (54), Zx -0.14 (94) | [92 genes](data/introgression_fix1/qtl_genes/SL_qtl_genes.tsv) | 5.9–6.5 (26) | 0.16 |
| Q3 | chr4 198.8 Mb | 183.0–205.5 | 5.6 | +0.26 ± 0.05 | **6.0** | 0.4 | Zd +0.14 (15), Zl +0.10 (5), Zv +0.14 (31), Zx +0.36 (68) | [641 genes](data/introgression_fix1/qtl_genes/SL_qtl_genes.tsv) | 183.4–199.7 (158) | 0.02 |
| Q4 | chr5 2.5 Mb | 1.8–2.7 | 5.8 | +0.25 ± 0.05 | **6.5** | 1.3 | Zd +0.15 (30), Zh +0.91 (5), Zl -0.18 (12), Zv +0.20 (42), Zx +0.29 (70) | [76 genes](data/introgression_fix1/qtl_genes/SL_qtl_genes.tsv) | 2.5–2.7 (22) | 0.39 |
| Q5 | chr5 49.4 Mb | 47.9–52.0 | 5.2 | -0.18 ± 0.10 | 1.2 | 1.4 | Zd -0.28 (10), Zl -0.20 (12), Zv -0.34 (29), Zx +0.03 (41) | [72 genes](data/introgression_fix1/qtl_genes/SL_qtl_genes.tsv) | 48.2–49.9 (18) | 0.04 |
| Q6 | chr5 64.2 Mb | 63.3–67.4 | 5.0 | -0.03 ± 0.11 | 0.1 | 0.2 | Zd +0.12 (8), Zl -0.15 (11), Zv -0.17 (34), Zx -0.12 (41) | [96 genes](data/introgression_fix1/qtl_genes/SL_qtl_genes.tsv) | – | 0.00 |
| Q7 | chr5 161.8 Mb | 157.2–164.8 | 7.0 | -0.27 ± 0.07 | **3.8** | 0.8 | Zd -0.59 (7), Zl +0.01 (8), Zv -0.16 (35), Zx -0.25 (49) | [129 genes](data/introgression_fix1/qtl_genes/SL_qtl_genes.tsv) | 157.6–161.8 (27) | 0.04 |
| Q8 | chr6 116.2 Mb | 113.6–118.1 | 4.9 | -0.26 ± 0.06 | **4.7** | 0.6 | Zd -0.28 (12), Zh -0.51 (7), Zl -0.04 (9), Zv -0.21 (22), Zx -0.17 (38) | [117 genes](data/introgression_fix1/qtl_genes/SL_qtl_genes.tsv) | 115.1–116.5 (41) | 0.05 |
| Q9 | chr6 167.8 Mb | 166.6–171.7 | 4.5 | -0.22 ± 0.06 | **3.6** | 0.1 | Zd -0.32 (17), Zh -0.45 (8), Zl -0.33 (8), Zv -0.23 (29), Zx -0.16 (43) | [256 genes](data/introgression_fix1/qtl_genes/SL_qtl_genes.tsv) | 167.4–171.0 (77) | 0.05 |
| Q10 | chr7 141.7 Mb | 139.6–141.9 | 9.2 | -0.34 ± 0.05 | **9.9** | 0.8 | Zd -0.26 (12), Zh +0.03 (7), Zl -0.60 (9), Zv -0.39 (38), Zx -0.35 (70) | [38 genes](data/introgression_fix1/qtl_genes/SL_qtl_genes.tsv) | 140.6–141.8 (20) | 0.15 |
| Q11 | chr9 157.2 Mb | 155.9–160.5 | 7.5 | -0.24 ± 0.04 | **7.3** | 2.2 | Zd -0.63 (20), Zh -0.17 (10), Zl -0.32 (13), Zv -0.27 (44), Zx -0.15 (94) | [263 genes](data/introgression_fix1/qtl_genes/SL_qtl_genes.tsv) | – | 0.10 |
| Q12 | chr10 25.6 Mb | 25.0–30.7 | 7.1 | +0.23 ± 0.07 | **3.1** | 2.0 | Zd +0.58 (22), Zv +0.10 (21), Zx +0.11 (56) | [71 genes](data/introgression_fix1/qtl_genes/SL_qtl_genes.tsv) | 25.3–30.0 (26) | 0.04 |
| Q13 | chr10 129.7 Mb | 125.9–131.1 | 7.2 | +0.20 ± 0.06 | **3.2** | 0.1 | Zd +0.21 (31), Zh +0.43 (8), Zl -0.12 (5), Zv +0.23 (33), Zx +0.28 (61) | [158 genes](data/introgression_fix1/qtl_genes/SL_qtl_genes.tsv) | 127.7–129.9 (55) | 0.05 |

Bonferroni over the 13 QTL: −log10 p > 2.41 (bold). Top windows: the fewest windows holding half of the QTL's fine-mapping probability (section 3.4.3). Every gene in each interval, with its GO annotation, the k-mer GWAS −log10 p of its window and its window's PIP: [`data/introgression_fix1/qtl_genes/SL_qtl_genes.tsv`](data/introgression_fix1/qtl_genes/SL_qtl_genes.tsv).

</details>

#### 3.6.11 Blade length (BL; 2023 only)

![BL: k-mer ancestry GWAS and joint QTL model](figures/fix1/joint_qtl_BL.png)

6 regions plus one QTL added by the conditional scans (Q1 chr2 9.2 Mb, open triangle); 6 of the 7 hold (Bonferroni over 7: −log10 p > 2.15). **chr2 has two QTL** (Q1 and Q2), which were one region in the scan. The six held QTL are those of the first model plus chr3 184.4 Mb (Q3, just below the cut before), with chr4 238.3 Mb not holding (1.3).

| QTL | Known genes | Fine-mapping PIP | Effect (joint) | Joint −log10 p | Taxa differ | What it says |
|---|---|---|---|---|---|---|
| **Q2 chr2 35.1 Mb** | — | 0.20 (best window; top windows 35.0–35.3 Mb, 9 genes) | −1.52 | 7.8 | 2.0 (no) | No significant difference between taxa (Zd −0.71, Zl +1.59, Zv −1.57, Zx −1.93). Not in *luxurians*. The peak moved from 39.8 Mb (first model; still −log10 p 10.4 in the scan), so *FT-INTERACTING PROTEIN 1* (Zm00001eb079900, 40.14 Mb) is no longer in the interval (34.7–36.0 Mb). |
| **Q1 chr2 9.2 Mb** (added) | — | 0.07 (best window; no credible set) | −1.35 | 6.8 | 0.0 (no) | Shorter blades in every taxon (Zd −1.76, Zh −1.31, Zl −2.90, Zv −1.40, Zx −1.04). |
| **Q7 chr7 185.2 Mb** | Zm00001eb331630 (PRR37-like); Zm00001eb331880 (dlf1) | 0.021 (CS); 0.066 (**top**) | −1.19 | 4.9 | 1.0 (no) | No significant difference between taxa (Zd −0.05, Zl +0.51, Zv −1.49, Zx −1.62). Same region as EH Q10. |
| **Q5 chr4 198.1 Mb** | Zm00001eb197430 (ns2) | 0.019 (**top**) | −1.20 | 2.8 | 0.1 (no) | The PH/EH/SL chr4 region (Zd +2.12, Zh −2.91, Zl +0.93, Zv −0.03, Zx −1.80). |
| **Q4 chr3 204.2 Mb** | Zm00001eb153000 (CCD8); Zm00001eb153190 (ZCN12) | 0.011 (CS); 0.015 (**top**) | +0.92 | 2.2 | 0.4 (no) | No significant difference between taxa (Zd +0.78, Zh +0.39, Zl +0.04, Zv +0.90, Zx +1.22). |
| **Q3 chr3 184.4 Mb** | — | 0.06 (best window) | +0.94 | 2.1 | 0.7 (no) | Longer blades; holds just above the cut. *liguleless2* (Zm00001eb147220), named here in the first model, is outside this model's interval (182.9–185.8 Mb). |

The chr4 238.3 Mb peak (Q6, −log10 p 4.9 in the scan) does not hold in the joint model (1.3): given chr4 198 Mb and the other QTL, it is not an independent QTL.

<details><summary>All BL QTL with per-taxon effects and all genes in the interval</summary>

| QTL | Position | Interval (Mb) | Scan −log10 p | Effect per copy (joint) | Joint −log10 p | Taxa differ −log10 p | Per-taxon effects (carrier plots) | Genes in interval | Top windows, Mb (genes) | Best window PIP |
|---|---|---|---|---|---|---|---|---|---|---|
| Q1 (cond.) | chr2 9.2 Mb | 2.1–11.2 | 9.2 | -1.35 ± 0.26 | **6.8** | 0.0 | Zd -1.76 (35), Zh -1.31 (10), Zl -2.90 (15), Zv -1.40 (51), Zx -1.04 (90) | [484 genes](data/introgression_fix1/qtl_genes/BL_qtl_genes.tsv) | – | 0.07 |
| Q2 | chr2 35.1 Mb | 34.7–36.0 | 10.7 | -1.52 ± 0.27 | **7.8** | 2.0 | Zd -0.71 (32), Zl +1.59 (16), Zv -1.57 (50), Zx -1.93 (91) | [26 genes](data/introgression_fix1/qtl_genes/BL_qtl_genes.tsv) | 35.0–35.3 (9) | 0.20 |
| Q3 | chr3 184.4 Mb | 182.9–185.8 | 4.7 | +0.94 ± 0.35 | **2.1** | 0.7 | Zd +0.80 (12), Zh -1.24 (6), Zl +1.11 (15), Zv +2.03 (28), Zx +0.88 (67) | [85 genes](data/introgression_fix1/qtl_genes/BL_qtl_genes.tsv) | 183.0–184.5 (40) | 0.06 |
| Q4 | chr3 204.2 Mb | 195.5–206.5 | 4.8 | +0.92 ± 0.33 | **2.2** | 0.4 | Zd +0.78 (19), Zh +0.39 (9), Zl +0.04 (14), Zv +0.90 (40), Zx +1.22 (58) | [307 genes](data/introgression_fix1/qtl_genes/BL_qtl_genes.tsv) | 202.4–205.9 (84) | 0.03 |
| Q5 | chr4 198.1 Mb | 193.8–207.2 | 6.7 | -1.20 ± 0.38 | **2.8** | 0.1 | Zd +2.12 (15), Zh -2.91 (6), Zl +0.93 (5), Zv -0.03 (30), Zx -1.80 (66) | [339 genes](data/introgression_fix1/qtl_genes/BL_qtl_genes.tsv) | 194.0–198.4 (69) | 0.02 |
| Q6 | chr4 238.3 Mb | 236.9–240.2 | 4.9 | -0.78 ± 0.39 | 1.3 | 0.7 | Zd -3.58 (17), Zl -4.61 (6), Zv -1.89 (32), Zx -0.19 (62) | [68 genes](data/introgression_fix1/qtl_genes/BL_qtl_genes.tsv) | 237.7–239.6 (27) | 0.04 |
| Q7 | chr7 185.2 Mb | 183.5–185.8 | 5.5 | -1.19 ± 0.27 | **4.9** | 1.0 | Zd -0.05 (30), Zl +0.51 (8), Zv -1.49 (39), Zx -1.62 (87) | [80 genes](data/introgression_fix1/qtl_genes/BL_qtl_genes.tsv) | 184.3–185.4 (27) | 0.09 |

Bonferroni over the 7 QTL: −log10 p > 2.15 (bold). Top windows: the fewest windows holding half of the QTL's fine-mapping probability (section 3.4.3). Every gene in each interval, with its GO annotation, the k-mer GWAS −log10 p of its window and its window's PIP: [`data/introgression_fix1/qtl_genes/BL_qtl_genes.tsv`](data/introgression_fix1/qtl_genes/BL_qtl_genes.tsv).

</details>

#### 3.6.12 Blade width (BW; 2023 only)

![BW: k-mer ancestry GWAS and joint QTL model](figures/fix1/joint_qtl_BW.png)

12 regions, 7 hold (Bonferroni over 12: −log10 p > 2.38). **All effects are negative: teosinte alleles make narrower leaf blades**, as expected from teosinte's narrow leaves. Five of the first model's six held QTL are found again (chr4 237.6 Mb, chr5 39–41 Mb, chr6 146.4 Mb, chr7 136 Mb, chr9 0.8–1.8 Mb); chr4 166.3 Mb does not hold (1.9, as before), and chr9 159.8 Mb and chr10 149.7 Mb are new.

| QTL | Known genes | Fine-mapping PIP | Effect (joint) | Joint −log10 p | Taxa differ | What it says |
|---|---|---|---|---|---|---|
| **Q4 chr5 41.2 Mb** | — | 0.07 (best window) | −0.33 | 10.7 | 2.3 (no) | The strongest BW QTL; mostly *mexicana* (Zx −0.51, Zl −0.25, Zv −0.16, Zd +0.06). |
| **Q5 chr6 146.4 Mb** | — | 0.06 (best window) | −0.28 | 7.8 | 0.0 (no) | Alike in the taxa (Zd −0.25, Zv −0.28, Zx −0.31). |
| **Q10 chr9 0.8 Mb** | — | 0.26 (best window, 0.7–0.9 Mb, no gene) | −0.24 | 6.9 | 1.3 (no) | Mostly *mexicana* (Zx −0.34, Zd −0.15, Zv −0.06). |
| **Q3 chr4 237.6 Mb** | — | 0.12 (best window) | −0.24 | 6.6 | 0.0 (no) | Alike in the taxa (Zd −0.22, Zl −0.31, Zv −0.34, Zx −0.21). |
| **Q6 chr7 135.7 Mb** | — | 0.07 (best window) | −0.18 | 4.7 | 0.3 (no) | Alike in the taxa (Zd −0.10, Zh −0.20, Zl −0.09, Zv −0.28, Zx −0.15). Its interval (130.7–136.3 Mb) no longer reaches *tsh4* (Zm00001eb316740, 138.54 Mb), which was in it in the first model. |

The other held QTL are Q11 chr9 159.8 Mb (−0.15) and Q12 chr10 149.7 Mb (−0.16); both have their probability on a few windows (best window PIP 0.28 and 0.30). The chr6 133.1 Mb peak (*fea4*, Zm00001eb280500; −log10 p 4.86) is above the threshold, but the scan puts it in the same significant region as chr6 146.4 Mb and does not make it a QTL of its own.

<details><summary>All BW QTL with per-taxon effects and all genes in the interval</summary>

| QTL | Position | Interval (Mb) | Scan −log10 p | Effect per copy (joint) | Joint −log10 p | Taxa differ −log10 p | Per-taxon effects (carrier plots) | Genes in interval | Top windows, Mb (genes) | Best window PIP |
|---|---|---|---|---|---|---|---|---|---|---|
| Q1 | chr4 166.3 Mb | 165.1–167.6 | 5.6 | -0.16 ± 0.06 | 1.9 | 0.2 | Zd -0.00 (7), Zv -0.07 (30), Zx -0.18 (64) | [64 genes](data/introgression_fix1/qtl_genes/BW_qtl_genes.tsv) | 165.5–166.4 (18) | 0.09 |
| Q2 | chr4 182.7 Mb | 180.3–183.4 | 4.9 | -0.04 ± 0.06 | 0.3 | 0.4 | Zd -0.08 (13), Zh +1.01 (8), Zl +0.91 (5), Zv -0.03 (31), Zx -0.07 (68) | [83 genes](data/introgression_fix1/qtl_genes/BW_qtl_genes.tsv) | – | 0.00 |
| Q3 | chr4 237.6 Mb | 237.0–238.7 | 8.7 | -0.24 ± 0.05 | **6.6** | 0.0 | Zd -0.22 (19), Zl -0.31 (6), Zv -0.34 (29), Zx -0.21 (62) | [32 genes](data/introgression_fix1/qtl_genes/BW_qtl_genes.tsv) | 237.1–237.7 (11) | 0.12 |
| Q4 | chr5 41.2 Mb | 33.3–43.8 | 9.4 | -0.33 ± 0.05 | **10.7** | 2.3 | Zd +0.06 (12), Zl -0.25 (12), Zv -0.16 (29), Zx -0.51 (46) | [216 genes](data/introgression_fix1/qtl_genes/BW_qtl_genes.tsv) | 36.8–42.3 (41) | 0.07 |
| Q5 | chr6 146.4 Mb | 144.6–150.7 | 9.4 | -0.28 ± 0.05 | **7.8** | 0.0 | Zd -0.25 (12), Zv -0.28 (28), Zx -0.31 (49) | [157 genes](data/introgression_fix1/qtl_genes/BW_qtl_genes.tsv) | 146.4–148.7 (29) | 0.06 |
| Q6 | chr7 135.7 Mb | 130.7–136.3 | 4.4 | -0.18 ± 0.04 | **4.7** | 0.3 | Zd -0.10 (11), Zh -0.20 (8), Zl -0.09 (7), Zv -0.28 (44), Zx -0.15 (60) | [138 genes](data/introgression_fix1/qtl_genes/BW_qtl_genes.tsv) | 132.7–136.3 (33) | 0.07 |
| Q7 | chr8 24.8 Mb | 22.5–30.8 | 4.5 | -0.09 ± 0.09 | 0.5 | 1.1 | Zd -0.04 (31), Zh +0.03 (11), Zl +0.41 (5), Zv +0.07 (47), Zx -0.34 (54) | [172 genes](data/introgression_fix1/qtl_genes/BW_qtl_genes.tsv) | – | 0.00 |
| Q8 | chr8 61.6 Mb | 60.3–64.6 | 4.5 | -0.04 ± 0.11 | 0.1 | 1.2 | Zd -0.04 (25), Zh -0.18 (9), Zl +0.18 (6), Zv -0.16 (43), Zx +0.08 (57) | [56 genes](data/introgression_fix1/qtl_genes/BW_qtl_genes.tsv) | – | 0.00 |
| Q9 | chr8 74.5 Mb | 64.8–84.1 | 4.3 | -0.08 ± 0.09 | 0.5 | 1.3 | Zd -0.03 (25), Zh -0.48 (11), Zl +0.18 (6), Zv -0.13 (39), Zx -0.03 (61) | [335 genes](data/introgression_fix1/qtl_genes/BW_qtl_genes.tsv) | – | 0.02 |
| Q10 | chr9 0.8 Mb | 0.7–1.8 | 5.8 | -0.24 ± 0.04 | **6.9** | 1.3 | Zd -0.15 (19), Zv -0.06 (30), Zx -0.34 (73) | [5 genes](data/introgression_fix1/qtl_genes/BW_qtl_genes.tsv) | 0.7–0.9 (0) | 0.26 |
| Q11 | chr9 159.8 Mb | 158.9–160.2 | 4.3 | -0.15 ± 0.04 | **4.8** | 0.2 | Zd -0.15 (24), Zh -0.26 (6), Zl -0.31 (14), Zv -0.07 (35), Zx -0.14 (92) | [69 genes](data/introgression_fix1/qtl_genes/BW_qtl_genes.tsv) | 159.6–159.9 (14) | 0.28 |
| Q12 | chr10 149.7 Mb | 149.4–150.2 | 4.3 | -0.16 ± 0.04 | **4.9** | 0.3 | Zd -0.04 (37), Zh -0.74 (5), Zl -0.33 (12), Zv -0.21 (39), Zx -0.15 (79) | [55 genes](data/introgression_fix1/qtl_genes/BW_qtl_genes.tsv) | 149.6–149.9 (15) | 0.30 |

Bonferroni over the 12 QTL: −log10 p > 2.38 (bold). Top windows: the fewest windows holding half of the QTL's fine-mapping probability (section 3.4.3). Every gene in each interval, with its GO annotation, the k-mer GWAS −log10 p of its window and its window's PIP: [`data/introgression_fix1/qtl_genes/BW_qtl_genes.tsv`](data/introgression_fix1/qtl_genes/BW_qtl_genes.tsv).

</details>

#### 3.6.13 Regions shared between traits

| Region | Traits (sign of the teosinte effect) | Known genes in the intervals |
|---|---|---|
| chr10 ~92–106 Mb | DTS +, DTA +, PH +, EH +, NBR + | Zm00001eb418700 (*ZmCCT10*) |
| chr10 ~146–150 Mb | SPAD +, EN +, BW −, PredN + (its best region, below its threshold) | Zm00001eb432440 (FRS6-like), Zm00001eb432850 (chlorophyll a/b-binding) |
| chr5 ~217–222 Mb | DTS +, DTA +, ASI + (mostly *luxurians*) | Zm00001eb257680 (HDR1-like), Zm00001eb255670 (*mads73*) |
| chr9 ~1–14 Mb | BW −, NBR −; DTA − and DTS − (above the scan threshold for DTA, not holding) | Zm00001eb372960 (*MADS74*) |
| chr1 ~22–24 Mb | EN +, Prolif + (NBR −, not holding) | Zm00001eb007670 (*ZmD14*) |
| chr3 ~11–19 Mb | Prolif +, PH −, EH − (EN +, not holding) | Zm00001eb123060 (*ra2*), Zm00001eb122680 (CYP734A6) |
| chr4 ~100–135 Mb | EN +, Prolif + (conditional QTL) | |
| chr4 ~171 Mb | Prolif + (EN scan peak at 167 Mb, BW − not holding) | |
| chr2 ~6–52 Mb | Prolif + (two QTL), BL − (two QTL) | Zm00001eb077080 (*GA2ox2*), Zm00001eb076990 (*ckx12*) |
| chr2 ~18–30 Mb | ASI −, EH +, EN + | Zm00001eb077080 (*GA2ox2*), Zm00001eb076990 (*ckx12*) (in the ASI interval) |
| chr2 ~202–206 Mb | EN +, Prolif + (same window) | |
| chr7 ~131–146 Mb | PH −, SL −, BW − | Zm00001eb316740 (*tsh4*) |
| chr7 ~165 Mb | ASI −, EN + | Zm00001eb323880 (*bhlh90*) |
| chr4 ~194–207 Mb | PH +, EH +, SL +, BL −, ASI − | Zm00001eb197430 (*ns2*) |
| chr5 ~2–3 Mb | PH +, EH +, SL + | |
| chr8 ~126–129 Mb | EH +, NBR + | Zm00001eb353250 (*ZCN8*) |
| chr7 ~183–186 Mb | EH +, BL − (DTS + in 2025) | Zm00001eb331880 (*dlf1*), Zm00001eb331630 (PRR37-like) |
| chr9 ~139–141 Mb | PH + (NBR + in the first model) | Zm00001eb397240 (*PRR73*) |

All genes in every region, with their GO terms and the −log10 p of the window they lie in, are in `data/introgression_fix1/GO/genes_GO_all_<TRAIT>.tsv` (one file per trait).

#### 3.6.14 *mexicana* split by race (SUBTAXON)

The same two analyses were rerun with *mexicana* split into its four races (Chalco, Mesa Central, Durango, Nobogame; section 2.1) instead of one Zx group: the fixed effect, the permutation strata and the per-taxon effects all use the race. Figures: `figures/fix1/subtaxon/joint_qtl_<TRAIT>.png`; tables: `data/introgression_fix1/assoc_gridlmm_spats_subtaxon/` and `data/introgression_fix1/joint_qtl/noperm_subtaxon/`.

| Trait | Threshold, taxon / race | QTL, taxon / race | Hold in joint model, taxon / race |
|---|---|---|---|
| DTS | 5.50 / 5.19 | 2 / 3 | 2 / 3 |
| DTA | 5.00 / 5.19 | 6 / 6 | 3 / 3 |
| ASI | 4.42 / 4.30 | 10 / 10 | 6 / 6 |
| PH | 4.17 / 4.10 | 15 / 13 | 10 / 11 |
| EH | 4.11 / 4.19 | 16 / 15 | 9 / 9 |
| EN | 4.49 / 4.43 | 18 / 17 | 12 / 12 |
| Prolif | 4.72 / 4.77 | 15 / 13 | 9 / 9 |
| NBR | 4.04 / 4.31 | 8 / 5 | 3 / 2 |
| SL | 4.20 / 4.17 | 13 / 13 | 11 / 11 |
| BL | 4.15 / 4.33 | 7 / 8 | 6 / 5 |
| BW | 4.21 / 4.27 | 12 / 11 | 7 / 6 |

**What changes:**
- **The QTL are essentially the same.** Most differences are a peak moving by one or two windows (e.g. PH chr7 137.8 → 141.7 Mb, PH chr4 198.4 → 197.8 Mb) or a broad region being split differently at a slightly different threshold. The races run adds DTS chr1 84.5 Mb (joint −log10 p 5.1), which is DTA Q1 and suggestive for DTS in the taxon run. So splitting *mexicana* does not create or remove signals; treating it as one taxon was not hiding race-specific QTL.
- **DTS / DTA chr5 221.8 Mb:** within *mexicana*, the effect comes from **Chalco** (+1.07 ± 0.41 d for DTS, 8 carrier plots); Mesa Central (+0.21), Nobogame (+0.19) and Durango (+0.09) are near zero. Chalco is the highland race (median 2,390 m). This is the same QTL where *luxurians* and *diploperennis* have the largest effects.
- **DTA chr9 10.9 Mb** (teosinte flowers earlier; not holding in this model): Durango −0.70 d and Nobogame −0.68 d, Mesa Central −0.10, Chalco +0.15.
- ***ZmCCT10* (chr10 97.8 Mb):** all four races delay flowering (DTS: Durango +1.02, Chalco +0.95, Nobogame +0.94, Mesa Central +0.61 d).
- **Caution:** a race has 6–40 carrier plots at a QTL, and Durango (3 accessions) and Nobogame (2) come from very few donor plants, so a race effect can be one accession's allele.

<details><summary>Figures with the races (DTS, DTA)</summary>

![DTS by race](figures/fix1/subtaxon/joint_qtl_DTS.png)

![DTA by race](figures/fix1/subtaxon/joint_qtl_DTA.png)

</details>

#### 3.6.15 Other traits

Scanned with the same pipeline; figures in `figures/fix1/joint_qtl_<TRAIT>.png`, tables in the same folders as above.

| Trait | What it is | Plots | 5% threshold | QTL | Hold in joint model | Strongest QTL (−log10 p) |
|---|---|---|---|---|---|---|
| PredN | leaf nitrogen predicted from leaf spectra (scans with Mahalanobis distance > 3 dropped; mean of the rest per plot); 2025, field B5 block 2, one plot per line; H2 0.25 (`final_script/28b_spats_nitrogen.R`) | 698 | 4.23 | 0 | 0 | chr10 146.9 Mb (2.8) |
| LAE | number of leaves above the ear; 2025 | 1248 | 4.23 | 9 | 4 | chr2 34.6 Mb (9.1) |
| Kinki | how much the stem bends, score 0–2; 2023 | 1118 | 4.23 | 6 | 4 | chr2 30.1 Mb (11.5) |
| StPu | stem pubescence (present/absent; line mean) | 1260 | 6.91 | 1 | 1 | chr4 240.8 Mb (8.7) |
| StPi | stem pigmentation (present/absent; line mean) | 1260 | 4.72 | 9 | 4 | chr2 18.1 Mb (11.1) |

StPu and StPi are scored present/absent, so their line means are mostly 0 with a few high values, and they agree poorly between years (r = 0.16 and 0.09). With such a skewed trait, a few plots can give very small p-values by chance, which is why the permutation threshold for StPu (6.91) is far above the others; the permutations account for this. A model for binary traits (e.g. a logistic mixed model) would suit them better.

### 3.7 What is new here

- **Local ancestry from k-mers in pooled, ultra-low-coverage samples.** Each sample is a pool of 6 BC2S3 plants sequenced at ~0.8×. The k-mer HMM (Section 2) turns presence and absence of B73 single-copy k-mers into B73 / SEG / DONOR calls per window, with a pool-aware emission model. k-mer GWAS methods use k-mers directly as markers in deeply sequenced individuals; here k-mers are turned into ancestry first.
- **The same ancestry calls give the genotype and the kinship.** The LOCO kinship is built from the ancestry dosage, so relatedness is measured in shared donor tracts, with no SNP calling or imputation.
- **Scan and joint model on the same QTL.** A k-mer ancestry GWAS finds the peaks; a joint QTL model on those peaks gives each QTL's effect adjusted for the others, its effect in each teosinte taxon, and a test that the taxa differ. The dominance test comes from the SEG state.

The statistics themselves (mixed models, LOCO, permutations, R/qtl2) are standard and cited above; what is new is putting k-mer ancestry from pooled skim-seq into them. To be checked against the literature before publication.

### 3.8 Limitations

- **Resolution** is set by tract length: regions are Mb wide, and genes inside them are candidates only.
- **λ_GC** is high under LOCO; the permutation threshold controls the genome-wide error, but linked windows on the same chromosome are not independent evidence, and two nearby peaks may be one QTL.
- **The pool is not the plot.** The genotype is 6 plants of a plot; the phenotype is the plot mean. For SEG plots, plants differ; the dosage of 1 is the expectation.
- **Donor reference.** DONOR levels start from MaizeGDB assemblies and are then learned from our pools, but not from our donor plants themselves (sections 2.6.6 and 2.13).
- **The QTL list depends on which window leads a long region.** Several traits are significant over tens of Mb (e.g. EN and Prolif on chr2 and chr4), and the joint model starts from one window per region (plus conditional QTL). Where two peaks in one region are nearly equal, a small change in the calls moves the lead (e.g. EN chr4 167.0 → 119.3 Mb, Prolif chr2 191.3 → 205.2 Mb between the two models, with the GWAS scans themselves correlated at r ≥ 0.97). A QTL that moved this way is a property of the region, not a new locus.
- **Replication:** DTS was checked in each year separately (section 3.6.1: both main QTL replicate). The same per-year runs exist for every trait (`*_y2023`, `*_y2025`) but are not yet described here; none within taxon.
- **Small taxa:** *luxurians* (Zl) and *huehuetenangensis* (Zh) have 5–25 carrier plots at most QTL, so their per-taxon effects have wide CIs. The large *luxurians*-specific effects (DTS/DTA chr5 221.8 Mb) should be checked against how many accessions and F1 families the carriers come from.
- **How peaks are split into regions.** A region is a run of windows above the threshold, joined when less than 5 Mb apart, so one region can hold more than one QTL. The conditional scans of the joint model (section 3.4) test this. They added 6 QTL in the taxon run (BL chr2 9.2; EN chr4 183.6; PH chr10 6.5; Prolif chr2 34.4, chr4 114.3; StPi chr2 163.0 Mb) and 4 in the race run, and found that other peaks inside regions (e.g. BW chr6 133 Mb, ASI chr3 96 Mb) are not independent QTL.
- **GO terms:** some of the 15,577 genes in the regions have no GO annotation in the lookup; they are listed in the gene tables without GO terms.

---

## Section 4 — References

All references cited in this README, in alphabetical order.

- Aguirre-Liguori JA, Gaut BS, Jaramillo-Correa JP, Tenaillon MI, Montes-Hernández S, García-Oliva F, Hearne SJ, Eguiarte LE (2019) Divergence with gene flow is driven by local adaptation to temperature and soil phosphorus concentration in teosinte subspecies (*Zea mays parviglumis* and *Zea mays mexicana*). *Molecular Ecology* 28:2814–2830. doi:10.1111/mec.15098
- Ahmed HI, et al. (2023) Einkorn genomics sheds light on history of the oldest domesticated wheat. *Nature* 620:830–838. (IBSpy used to find einkorn introgressions in bread wheat.)
- Ashikari M, Sakakibara H, Lin S, Yamamoto T, Takashi T, Nishimura A, Angeles ER, Qian Q, Kitano H, Matsuoka M (2005) Cytokinin oxidase regulates rice grain production. *Science* 309:741–745. (Cited from Ni et al. 2024.)
- An H, Li K, Liu X, Wu Y, He J, Wang Y, Gao Z (2026) Identification of QTLs and new candidate genes affecting ear shank length via BSA-seq and transcriptomic analysis in maize. *Front Plant Sci* 17:1768852. doi:10.3389/fpls.2026.1768852
- Baum LE, Petrie T, Soules G, Weiss N (1970) A maximization technique occurring in the statistical analysis of probabilistic functions of Markov chains. *Ann Math Stat* 41:164–171.
- Bolger AM, Lohse M, Usadel B (2014) Trimmomatic: a flexible trimmer for Illumina sequence data. *Bioinformatics* 30:2114–2120.
- Bortiri E, Chuck G, Vollbrecht E, Rocheford T, Martienssen R, Hake S (2006) *ramosa2* encodes a LATERAL ORGAN BOUNDARY domain protein that determines the fate of stem cells in branch meristems of maize. *Plant Cell* 18:574–585. doi:10.1105/tpc.105.039032
- Brutnell TP, Sawers RJH, Mant A, Langdale JA (1999) BUNDLE SHEATH DEFECTIVE2, a novel protein required for post-translational regulation of the *rbcL* gene of maize. *Plant Cell* 11:849–864.
- Broman KW, Gatti DM, Simecek P, et al. (2019) R/qtl2: software for mapping quantitative trait loci with high-dimensional data and multiparent populations. *Genetics* 211:495–502.
- Broman KW, Sen Ś (2009) *A Guide to QTL Mapping with R/qtl*. Springer, New York.
- Brown CT, Irber L (2016) sourmash: a library for MinHash sketching of DNA. *J Open Source Softw* 1:27.
- Campos-Martin R, Schmickler S, Goel M, Schneeberger K, Tresch A (2023) Reliable genotyping of recombinant genomes using a robust hidden Markov model. *Plant Physiol* 192:821–836. doi:10.1093/plphys/kiad191
- Chen C, Begcy K, Liu K, Folsom JJ, Wang Z, Zhang C, Walia H (2016) Heat stress yields a unique MADS box transcription factor in determining seed size and thermal sensitivity. *Plant Physiol* 171:606–622. doi:10.1104/pp.15.01992
- Chen S, Gao S, Wang D, Liu J, Ren Y, Wang Z, Wei X, Wang Q, Huang X (2024) FKF1b controls reproductive transition associated with adaptation to geographical distribution in maize. *J Integr Plant Biol* 66:943–955. doi:10.1111/jipb.13639
- Chen Q, et al. (2019) TeoNAM: a nested association mapping population for domestication and agronomic trait analysis in maize. *Genetics* 213:1065–1078.
- Chuck GS, Brown PJ, Meeley R, Hake S (2014) Maize SBP-box transcription factors *unbranched2* and *unbranched3* affect yield traits by regulating the rate of lateral primordia initiation. *PNAS* 111:18775–18780. doi:10.1073/pnas.1407401112
- Churchill GA, Doerge RW (1994) Empirical threshold values for quantitative trait mapping. *Genetics* 138:963–971.
- Danecek P, Bonfield JK, Liddle J, et al. (2021) Twelve years of SAMtools and BCFtools. *GigaScience* 10:giab008.
- Devlin B, Roeder K (1999) Genomic control for association studies. *Biometrics* 55:997–1004.
- Haldane JBS (1919) The combination of linkage values and the calculation of distances between the loci of linked factors. *J Genet* 8:299–309.
- Hrytsenko Y, Daniels NM, Schwartz RS (2025) Determining population structure from k-mer frequencies. *PeerJ* 13:e18939. doi:10.7717/peerj.18939
- Huang C, Sun H, Xu D, et al. (2018) ZmCCT9 enhances maize adaptation to higher latitudes. *PNAS* 115:E334–E341.
- Hung H-Y, Shannon LM, Tian F, et al. (2012) ZmCCT and the genetic basis of day-length adaptation underlying the postdomestication spread of maize. *PNAS* 109:E1913–E1921.
- Irber L, Brooks PT, Reiter T, et al. (2022) Lightweight compositional analysis of metagenomes with FracMinHash and minimum metagenome covers. *bioRxiv* doi:10.1101/2022.01.11.475838
- Ito S, Song YH, Josephson-Day AR, Miller RJ, Breton G, Olmstead RG, Imaizumi T (2012) FLOWERING BHLH transcriptional activators control expression of the photoperiodic flowering regulator CONSTANS in *Arabidopsis*. *PNAS* 109:3582–3587. doi:10.1073/pnas.1118876109
- Kong D, Li C, Xue W, Wei H, Ding H, Hu G, Zhang X, Zhang G, Zou T, Xian Y, Wang B, Zhao Y, Liu Y, Xie Y, Xu M, Wu H, Liu Q, Wang H (2023) UB2/UB3/TSH4-anchored transcriptional networks regulate early maize inflorescence development in response to simulated shade. *Plant Cell* 35:717–737. doi:10.1093/plcell/koac352
- Kokot M, Długosz M, Deorowicz S (2017) KMC 3: counting and manipulating k-mer statistics. *Bioinformatics* 33:2759–2761.
- Korneliussen TS, Albrechtsen A, Nielsen R (2014) ANGSD: Analysis of Next Generation Sequencing Data. *BMC Bioinformatics* 15:356. doi:10.1186/s12859-014-0356-4
- Lander ES, Waterman MS (1988) Genomic mapping by fingerprinting random clones: a mathematical analysis. *Genomics* 2:231–239. (Coverage formula behind g.)
- Lee JH, Yoo SJ, Park SH, Hwang I, Lee JS, Ahn JH (2007) Role of SVP in the control of flowering time by ambient temperature in *Arabidopsis*. *Genes Dev* 21:397–402.
- Li H (2013) Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. *arXiv* 1303.3997.
- Li M, Sun M, et al., Zhang A, Lu C (2026) Uncovering the role of the PPR protein PHOTOSYSTEM ONE BIOGENESIS FACTOR6 in splicing chloroplast group II introns. *J Integr Plant Biol*. doi:10.1111/jipb.70289
- Li X, Zhang M, Huang Y, Li M, Zhu H, Li S, Shi D, Xu L, Wang Q, Zhang J, et al., Wang J, Liu H (2026) Polycomb protein ZmEMF1a restricts endosperm proliferation and directs differentiation via stage-specific H2Aub1 and H3K27me3 landscapes in maize. *New Phytol* 252:1130–1148. doi:10.1111/nph.71543
- Lippert C, Listgarten J, Liu Y, Kadie CM, Davidson RI, Heckerman D (2011) FaST linear mixed models for genome-wide association studies. *Nat Methods* 8:833–835.
- Listgarten J, Lippert C, Kadie CM, Davidson RI, Eskin E, Heckerman D (2012) Improved linear mixed models for genome-wide association studies. *Nat Methods* 9:525–526.
- Marín-González E, Matías-Hernández L, Aguilar-Jaramillo AE, Lee JH, Ahn JH, Suárez-López P, Pelaz S (2015) SHORT VEGETATIVE PHASE up-regulates TEMPRANILLO2 floral repressor at low ambient temperatures. *Plant Physiol* 169:1214–1224. doi:10.1104/pp.15.00570
- Méndez-Vigo B, Martínez-Zapater JM, Alonso-Blanco C (2013) The flowering repressor SVP underlies a novel *Arabidopsis thaliana* QTL interacting with the genetic background. *PLoS Genet* 9:e1003289. doi:10.1371/journal.pgen.1003289
- Nardmann J, Ji J, Werr W, Scanlon MJ (2004) The maize duplicate genes *narrow sheath1* and *narrow sheath2* encode a conserved homeobox gene function in a lateral domain of shoot apical meristems. *Development* 131:2827–2839.
- Neupane SP, et al. (2025) Genetic basis of Fusarium ear rot resistance and productivity traits in a heterozygous multi-parent recombinant inbred intercross (RIX) maize population. *BMC Plant Biol* 25:639. (MaizeGDB annotation of *bhlh90*.)
- Ni J, Tang D, Chen Z, Yang S, Wang X, Liu Z, Deng W, Wu H, You C, Yang J, Meng P, Bao R, Rong T, Liu J (2024) Identification and segregation of two closely linked major QTLs for kernel row number in advanced maize–teosinte populations. *Theor Appl Genet* 137:271. doi:10.1007/s00122-024-04780-0
- Oakey H, Verbyla A, Pitchford W, Cullis B, Kuchel H (2006) Joint modeling of additive and non-additive genetic line effects in single field trials. *Theor Appl Genet* 113:809–819.
- Odell SG, Hudson AI, Praud S, Dubreuil P, Tixier M-H, Ross-Ibarra J, Runcie DE (2022) Modeling allelic diversity of multiparent mapping populations affects detection of quantitative trait loci. *G3* 12:jkac011.
- Picard toolkit. Broad Institute. https://broadinstitute.github.io/picard/
- Poretsky E, Cagirici HB, Andorf CM, Sen TZ (2024) Harnessing the predicted maize pan-interactome for putative gene function prediction and prioritization of candidate genes for important traits. *G3* 14:jkae059. doi:10.1093/g3journal/jkae059
- Portwood JL 2nd, Woodhouse MR, Cannon EK, et al. (2019) MaizeGDB 2018: the maize multi-genome genetics and genomics database. *Nucleic Acids Res* 47:D1146–D1154. (Source of the teosinte assemblies.)
- Quiroz-Chávez J (2023) Wheat haplotype diversity by a k-mer based approach. PhD thesis, University of East Anglia. (IBSpy.)
- Rabiner LR (1989) A tutorial on hidden Markov models and selected applications in speech recognition. *Proc IEEE* 77:257–286.
- Rahman A, Hallgrímsdóttir I, Eisen M, Pachter L (2018) Association mapping from sequencing reads using k-mers. *eLife* 7:e32920.
- Rodríguez-Álvarez MX, Boer MP, van Eeuwijk FA, Eilers PHC (2018) Correcting for spatial heterogeneity in plant breeding experiments with P-splines. *Spatial Statistics* 23:52–71.
- Roose JL, Frankel LK, Bricker TM (2014) The PsbP domain protein 1 functions in the assembly of lumenal domains in photosystem I. *J Biol Chem* 289:23776–23785. doi:10.1074/jbc.M114.589085
- Rowan BA, Patel V, Weigel D, Schneeberger K (2015) Rapid and inexpensive whole-genome genotyping-by-sequencing for crossover localization and fine-scale genetic mapping. *G3* 5:385–398.
- Runcie DE, Crawford L (2019) Fast and flexible linear mixed models for genome-wide genetics. *PLoS Genet* 15:e1007978.
- sabre: barcode demultiplexing and trimming tool for FastQ files. https://github.com/najoshi/sabre
- Selvanayagam S, Quiroz-Chavez J, Ramirez-Gonzalez RH, Uauy C, Smit S, Schranz ME (2026) KCFtools: rapid alignment-free method for introgression screening and GWAS using k-mer profiles. *Bioinformatics*. doi:10.1093/bioinformatics/btag354
- Smith LH, Langdale JA, Chollet R (1998) A functional Calvin cycle is not indispensable for the light activation of C4 phosphoenolpyruvate carboxylase kinase and its target enzyme in the maize mutant *bundle sheath defective2-mutable1*. *Plant Physiol* 118:191–197.
- Smith MW, O'Brien SJ (2005) Mapping by admixture linkage disequilibrium: advances, limitations and guidelines. *Nat Rev Genet* 6:623–632.
- Stitzer MC, Seetharam AS, et al., Hufford MB (2025) Extensive genome evolution distinguishes maize within a stable tribe of grasses. *bioRxiv* doi:10.1101/2025.01.22.633974. (PanAnd assemblies.)
- Studer AJ, Wang H, Doebley JF (2017) Selection during maize domestication targeted a gene network controlling plant and inflorescence architecture. *Genetics* 207:755–765. doi:10.1534/genetics.117.300071
- Sun D, Chen S, Cui Z, Lin J, Liu M, Jin Y, Zhang A, Gao Y, Cao H, Ruan Y (2022) Genome-wide association study reveals the genetic basis of brace root angle and diameter in maize. *Front Genet* 13:963852. doi:10.3389/fgene.2022.963852
- Sun X, Zhang Z, Wu J, Cui X, Feng D, Wang K, Xu M, Zhou L, Han X, Gu X, Lu T (2016) The *Oryza sativa* regulator HDR1 associates with the kinase OsK4 to control photoperiodic flowering. *PLoS Genet* 12:e1005927. doi:10.1371/journal.pgen.1005927
- Tang W, Wang W, Chen D, Ji Q, Jing Y, Wang H, Lin R (2012) Transposase-derived proteins FHY3/FAR1 interact with PHYTOCHROME-INTERACTING FACTOR1 to regulate chlorophyll biosynthesis by modulating *HEMB1* during deetiolation in *Arabidopsis*. *Plant Cell* 24:1984–2000.
- Teixeira PF, Kmiec B, Branca RMM, Murcha MW, Byzia A, Ivanova A, Whelan J, Drag M, Lehtiö J, Glaser E (2016) A multi-step peptidolytic cascade for amino acid recovery in chloroplasts. *Nat Chem Biol*. doi:10.1038/nchembio.2227
- Viterbi AJ (1967) Error bounds for convolutional codes and an asymptotically optimum decoding algorithm. *IEEE Trans Inf Theory* 13:260–269.
- Voichek Y, Weigel D (2020) Identifying genetic variants underlying phenotypic variation in plants without complete genomes. *Nat Genet* 52:534–540.
- Waditee-Sirisattha R, Shibato J, Rakwal R, Sirisattha S, Hattori A, Nakano T, Takabe T, Tsujimoto M (2011a) The *Arabidopsis* aminopeptidase LAP2 regulates plant growth, leaf longevity and stress response. *New Phytol* 191:958–969. doi:10.1111/j.1469-8137.2011.03758.x
- Waditee-Sirisattha R, Hattori A, Shibato J, Rakwal R, Sirisattha S, Takabe T, Tsujimoto M (2011b) Role of the *Arabidopsis* leucine aminopeptidase 2. *Plant Signal Behav* 6:1581–1583. doi:10.4161/psb.6.10.17105
- Wang D, Gao Z, Dong M, et al. (2026) Integrated GWAS and Fst analysis identify flowering-time regulatory genes in maize. *Plant Genome* 19:e70185. doi:10.1002/tpg2.70185
- Wang G, Sarkar A, Carbonetto P, Stephens M (2020) A simple new approach to variable selection in regression, with application to genetic fine mapping. *J R Stat Soc Series B* 82:1273–1300. doi:10.1111/rssb.12388
- Wang H, Studer AJ, Zhao Q, Meeley R, Doebley JF (2015) Evidence that the origin of naked kernels during maize domestication was caused by a single amino acid substitution in *tga1*. *Genetics* 200:965–974. doi:10.1534/genetics.115.175752
- Wang X, Li J, Han L, Liang C, Li J, Shang X, Miao X, Luo Z, Zhu W, Li Z, Li T, Qi Y, Li H, Lu X, Li L (2023) QTG-Miner aids rapid dissection of the genetic base of tassel branch number in maize. *Nat Commun* 14:5232. doi:10.1038/s41467-023-41022-1
- Xiao Y, Guo J, Dong Z, Richardson A, Patterson E, Mangrum S, Bybee S, Bertolini E, Bartlett M, et al., Whipple C (2022) Boundary domain genes were recruited to suppress bract growth and promote branching in maize. *Sci Adv* 8:eabm6835. doi:10.1126/sciadv.abm6835
- Xie Y, Zhao Y, Chen L, Wang Y, Xue W, Kong D, Li C, Zhou L, Li H, Zhao Y, Wang B, Xu M, Zhao B, Bilska-Kos A, Wang H (2023) ZmELF3.1 integrates the RA2-TSH4 module to repress maize tassel branching. *New Phytol*. doi:10.1111/nph.19329
- Zhang J, Verbraeken L, Sprenger H, et al., Maere S, Nelissen H (2025) High-resolution transcriptional atlas of growing maize shoot organs throughout plant development under well-watered and drought conditions. *bioRxiv*. doi:10.1101/2025.03.12.642568
- Zhong T, Mullens A, Morales L, Swarts KL, Stafstrom WC, He Y, Sermons SM, Yang Q, Lopez-Zuniga LO, Rucker E, Thomason WE, Nelson RJ, Jamann T, Balint-Kurti PJ, Holland JB (2025) A maize near-isogenic line population designed for gene discovery and characterization of allelic effects. *Plant J* 122:e70228. doi:10.1111/tpj.70228
- Yang J, Weedon MN, Purcell S, et al. (2011) Genomic inflation factors under polygenic inheritance. *Eur J Hum Genet* 19:807–812.
- Yang J, Zaitlen NA, Goddard ME, Visscher PM, Price AL (2014) Advantages and pitfalls in the application of mixed-model association methods. *Nat Genet* 46:100–106.
- Yang Q, Li Z, Li W, et al. (2013) CACTA-like transposable element in *ZmCCT* attenuated photoperiod sensitivity and accelerated the postdomestication spread of maize. *PNAS* 110:16969–16974.
- Yin Y, Yan Z, Guan J, Huo Y, Wang T, Li T, Cui Z, Ma W, Wang X, Chen W (2023) Two interacting basic helix-loop-helix transcription factors control flowering time in rice. *Plant Physiol* 192:205–221. doi:10.1093/plphys/kiad077
