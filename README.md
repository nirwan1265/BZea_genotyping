# BZea: k-mer introgression calling and QTL mapping in maize × teosinte BC2S3 pools

This repository goes from raw paired-end reads of **pooled, low-coverage (~0.8×)** BC2S3 samples to **teosinte introgression calls** and **QTL** for flowering, height, ear and leaf traits. B73 was crossed to 5 teosinte taxa, backcrossed twice to B73 and selfed three times; each sample is a pool of 6 plants from one BC2S3 plot.

- **Introgressions are called from k-mers**, not from mapped SNPs: for each pool, which B73 single-copy 31-mers its reads contain, per 100-kb window, decoded with an HMM that models the pool: each window is **B73** (all 6 plants AA), **DONOR** (all BB) or **SEG** (segregating: the plants are a mix of AA, AB and BB).
- **QTL are mapped on that k-mer ancestry**: a genome scan (k-mer ancestry GWAS, GridLMM with leave-one-chromosome-out kinship) finds the peaks, and a joint QTL model (R/qtl2) fits them together, with an effect for each teosinte taxon.

> **Pipeline and every command:** [`final_script/`](final_script/README.md). Results tables: [`data/introgression_spread/`](data/introgression_spread/). Figures: [`figures/`](figures/).

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
- **Sequenced sample:** one leaf each was taken from 6 different BC2S3 plants of the plot and pooled into one sample. These 6 plants are siblings, not copies of one plant: each got its own two gametes from the BC2S2 mother. So where the mother is heterozygous, each plant is AA, AB or BB with chance ¼, ½, ¼.
- **Pedigree names:** in `Zv.0200_P3_P2_P1.1.1.1`, `Zv.0200` is the teosinte accession, `P3` the F1 plant, `P2` the BC1 plant, `P1` the BC2 plant, and `.1.1.1` the selfed line. Pools that share `Zv.0200_P3` descend from one F1 plant (an **F1 family**: they all carry the same teosinte haplotype); pools that share `Zv.0200_P3_P2_P1` come from the same BC2 plant.

## Expected genotype proportions (per locus)

**Notation, used everywhere in this README and in `final_script/`:**
- **A** = the B73 allele, **B** = the teosinte (donor) allele
- **AA** = homozygous B73, **AB** = heterozygous, **BB** = homozygous teosinte

**Assumptions:** B73 is AA, the teosinte donor contributes one B allele to the F1, and there is no selection or segregation distortion. (In practice one plant was selected at each selfing step, so a departure from these numbers can point to selection; see [Results](#28-results).) These are the same for every donor taxon: Zd (*diploperennis*), Zh (*huehuetenangensis*), Zl (*luxurians*), Zv (*parviglumis*) and Zx (*mexicana*).

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

These are exactly the Step 3 numbers. This is how the HMM results are reported at both levels: the pool states as called (compared with Step 4), and converted to plant genotypes with this table (compared with Step 3). Heterozygosity (AB) can only be estimated this way, as half of SEG; no single plant is genotyped.

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
- **Expected DONOR level** per window from the teosinte genome assemblies (MaizeGDB / PanAnd).
- **HMM** with three states, B73 / SEG / DONOR, genetic-map transitions and pool-aware emissions (steps 20–21). Outputs: per-pool state per window, tracts, `pool_summary.tsv`, `family_summary.tsv`, chromosome paintings, interactive figures.
- **Population structure from k-mers** (sourmash; step 32).
- **SNP route** (ANGSD sites and read counts, GL-HMM; steps 01–15), kept as a cross-check.

### 3. k-mer ancestry GWAS and joint QTL mapping (Section 3)
- **Phenotypes:** SpATS spatially corrected BLUEs per year (2023, 2025), combined across years (step 28).
- **k-mer ancestry GWAS:** GridLMM, taxon as fixed effect, **LOCO kinship** from the ancestry dosage, 5% genome-wide threshold from 1,000 permutations (steps 25–27); rerun with *mexicana* split by race.
- **Joint QTL model:** the peaks fitted together in R/qtl2: effect of each QTL adjusted for the others, its effect in each taxon, and a test that the taxa differ (step 29).
- **Genes and GO terms** in each QTL interval, side-by-side figures (steps 30–31).

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

**What we added:** low-coverage *pools* instead of single plants, a three-state pool model (B73 / SEG / DONOR) with per-window DONOR levels from assemblies, and the HMM on the NAM genetic map.

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
5. **Allow for our donor plants differing from the assembly plants (step 20).** Where the two assemblies of a taxon disagree on d, the model is told that the DONOR level there is uncertain, so it doesn't flip between SEG and DONOR on small differences (see τ_w under [Emissions](#261-emissions-how-well-does-r-fit-each-state)).
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

The HMM combines both to find the most likely series of states along each chromosome. It has the same three-state structure as the SNP HMM in step 14 (`14_hmm_pooled.py`); see Rabiner (1989) for HMMs in general. Only the observation is different: a k-mer score per window instead of SNP reads.

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
- **Where to:** if the state changes, it moves to one of the other two states in proportion to their prior (e.g. from B73, more often to DONOR (9.4%) than to SEG (6.25%)).

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
| d_w (DONOR level) | donor taxon × window | MaizeGDB assemblies of that taxon | `19_kmer_donor_reference.sh` |
| τ_w (uncertainty of d_w) | donor taxon × window | disagreement of the taxon's two assemblies | `20_kmer_hmm.py` (`--d_spread`) |
| a (scale on d) | pool | best likelihood over 0.5–1.8 | `20_kmer_hmm.py` |
| ρ (switch rate) | whole data set | best total likelihood | `21_kmer_collect.py` |
| Prior | fixed | BC2S3 expectation | `config.sh` |
| Plants per pool | fixed | 6 | `config.sh` (`KMER_PLANTS`) |
| Outlier share | fixed | 2% | `20_kmer_hmm.py` (`--eps`) |

"Donor taxon" is the pool's teosinte (Zd, Zh, Zl, Zv or Zx, from the sample sheet). Every pool of that taxon uses the same d_w and τ_w, taken from the MaizeGDB assemblies of the taxon; nothing here comes from our BC1 or BC2 plants.

**Why the prior is fixed.** It only sets the state each chromosome starts in, and which state a switch goes to. It doesn't set how much of the genome ends up in each state: with about 2,000 windows per chromosome, the data decide. DONOR came out at 5.0% although the prior expects 9.4%, so the prior isn't forcing the result.

#### 2.6.5 How long are the tracts?

Tract length depends on how many meioses broke the donor chromosome up: here the F1, BC1 and BC2 plants, then two rounds of selfing before the BC2S2 mother of each pool. Each crossover cuts a donor segment, so after these meioses most segments are many Mb long, but a few are short. Short ones come mostly from regions with high recombination (near chromosome ends, where 500 kb can be more than 1 cM) and from crossovers in different generations falling close together.

To check that short tracts are expected and not just a feature of our data, the crosses were simulated on the NAM genetic map (F1 → BC1 → BC2 → BC2S1 → BC2S2, crossovers placed at random along the map, 400 lines per chromosome, no selection), and the lengths of the introgressed tracts (SEG and DONOR joined) were compared with our calls (1,292 pools):

| Introgressed tracts | Median length | < 0.5 Mb | < 1 Mb | < 2 Mb | < 5 Mb |
|---|---|---|---|---|---|
| Simulated BC2S2 mothers (3,506 tracts) | 11.3 Mb | 4.4% | 8.4% | 15.3% | 30.8% |
| Our pools, k-mer HMM (12,662 tracts) | 9.8 Mb | 2.6% | 6.5% | 13.7% | 33.2% |

- **Short tracts are expected from the crossing scheme alone**, not only in our data: about 1 in 12 introgressed tracts should be under 1 Mb.
- We find slightly fewer under 1 Mb than the simulation. That fits the HMM's detection limit: a piece needs about 4 windows (400 kb) to be called ([Transitions](#262-transitions-how-often-can-the-state-change)), so the shortest tracts are missed.
- **Any single short tract can still be an error**, e.g. a window run where our donor plant differs from the assemblies, or a structural variant. They are a small part of the genome: tracts under 1 Mb hold 0.16% of the introgressed sequence (0.11% in the simulation). See [Limitations](#213-limitations-and-future-work) for how they can be checked.

The simulation does not include crossover interference (which makes two crossovers close together in one meiosis rarer) or selection, so it is a rough guide, not an exact expectation. Script: `final_script/tract_length_sim.py` (seed 1).

### 2.7 How k-mers solve mapping bias, low depth and pooling

Back to the three problems of the SNP route ([above](#22-the-problem-mapping-reads-and-calling-snps-at-08x-in-pools)):

| Problem (above) | How the k-mer approach avoids it |
|---|---|
| **1. Donor reads are lost** | No mapping. A pool's raw reads are counted as k-mers, and only the presence of B73 k-mers is checked, so donor reads can't be filtered away. |
| **2. Too few reads per site** | No per-site calls. Each 100-kb window adds up ~32,000 k-mers, so its score r is measured precisely from thousands of positions at once. |
| **3. Pools aren't diploid** | The model is built for pools: DONOR = d, SEG = d + (1 − d)·g(f), averaged over the donor fractions a 6-plant family can have (Binomial(12, ½)). |

**Telling DONOR from SEG.** It no longer depends on allele ratios at single sites, which reference bias distorts. It depends on the window's **level**: in the Zd example above, DONOR ≈ 0.27 and SEG ≈ 0.66, with B73 = 1. The DONOR level is measured from teosinte genomes, not guessed.

### 2.8 Results

1,292 donor pools (the 93 low-depth pools are left out; see [Noisy pools](#211-noisy-pools-and-resequencing)).

| Taxon | Pools | Donor alleles | SEG (pool) | DONOR (pool) | AB (plant) | BB (plant) | AB/(AB+BB) |
|---|---|---|---|---|---|---|---|
| Zd | 227 | 5.6% | 5.1% | 3.0% | 2.6% | 4.3% | 37% |
| Zh | 61 | 8.2% | 5.7% | 5.3% | 2.9% | 6.7% | 30% |
| Zl | 98 | 7.2% | 6.8% | 3.8% | 3.4% | 5.5% | 38% |
| Zv | 333 | 8.4% | 6.5% | 5.1% | 3.2% | 6.8% | 32% |
| Zx | 573 | 9.0% | 6.5% | 5.8% | 3.2% | 7.4% | 30% |
| **All** | **1,292** | **8.1%** | **6.2%** | **5.0%** | **3.1%** | **6.5%** | **32%** |
| *Expected BC2S3* | | *12.5%* | *6.25%* | *9.4%* | *3.1%* | *10.9%* | *22%* |

> **What AB/(AB+BB) means.** Take only the plants that carry at least one teosinte allele at a locus: these are the AB plants and the BB plants. AB/(AB+BB) is the fraction of them that are **AB** (heterozygous). The rest, BB/(AB+BB), are homozygous teosinte. Example: out of 100 plants that carry teosinte at a locus, how many are AB.
> - **Expected for BC2S3:** AB = 3.13% and BB = 10.94% of the genome, so AB/(AB+BB) = 3.13 / (3.13 + 10.94) = **22%**: about 22 AB and 78 BB plants out of 100.
> - **Observed:** AB = 3.12% and BB = 6.51%, so AB/(AB+BB) = 3.12 / (3.12 + 6.51) = **32%**: about 32 AB and 68 BB out of 100.
> - **Why it's useful:** it doesn't depend on how much teosinte there is in total. If lines lost teosinte at random (drift, or segments lost before BC2), AB and BB would drop by the same factor and the ratio would stay at 22%. It goes up only if BB plants are lost more than AB plants.

![Observed vs expected, all taxa](final_script/figures/results_spread/fig1_overall.png)

The figure has three groups, because these numbers are different kinds of measures:
- **Pool states (SEG, DONOR):** what the HMM calls, per pool and window.
- **Plant genotypes (AB, BB):** what the pool states mean for single plants. They are **not measured** directly: a pool can't say which plant is AB. They are computed from the pool states using Mendel: a SEG region holds ¼ AA, ½ AB and ¼ BB plants, so **AB = SEG/2** and **BB = DONOR + SEG/4**. So AB lives inside SEG.
- **Donor alleles:** the total share of teosinte alleles, = DONOR + SEG/2.

![Observed vs expected by taxon](final_script/figures/results_spread/fig2_by_taxon.png)

![Per family](final_script/figures/results_spread/fig3_by_family.png)

![Per chromosome](final_script/figures/results_spread/fig4_by_chromosome.png)

![Donor share per pool](final_script/figures/results_spread/fig5_pool_share.png)

Figures from step 24 (`24_kmer_figures.py`); the per-family values are in `final_script/figures/results_spread/family_values.tsv`.

**What this shows:**
- **Het is as expected:** SEG (6.23%) and AB (3.12%) match the BC2S3 expectation; the F1-family bootstrap finds no difference (p = 0.93).
- **Homozygous donor is about 60% of the expectation** (BB 6.5% vs 10.9%; p < 0.0002 in every taxon). That loss is the whole shortfall in donor alleles (8.1% vs 12.5%).
- **So a larger share of the remaining donor is het:** 32% of donor-carrying plants are AB, against 22% expected. This is not extra het: AB is exactly as expected (3.1% of the genome). The share goes up only because BB went down, so AB is a bigger part of a smaller total. This holds on every chromosome (29–35%).
- **A likely reason: selection during selfing.** At every selfing step one plant per row was chosen and advanced (ear to row), and the plants with the best-looking ears were the ones kept. That selection works against BB more than against AB:
  - **Selecting against teosinte.** A plant homozygous for a teosinte segment (BB) shows the teosinte form of the traits that segment controls (e.g. smaller or harder ears), with no B73 copy to mask it. Such plants don't look as good as B73-like plants, so they are less often chosen. An AB plant still has a B73 copy, which often hides the teosinte effect, so it is chosen more often than a BB plant.
  - **Selecting for vigour, like hybrid vigour.** An AB plant is a small "hybrid" between B73 and teosinte at that segment, and AB plants can be more vigorous than either homozygote. Picking the best-looking ears is, in part, picking for that vigour, which keeps AB.
  - **Hidden harmful alleles.** Teosinte is outcrossing, so it carries recessive harmful alleles that do no harm in AB plants but show up in BB plants (the same reason maize inbreeding reduces vigour). This also removes BB plants.

  The result is fewer BB plants while AB stays as expected, which is exactly what raises AB/(AB+BB) and lowers the total donor share. TeoNAM (Chen et al. 2019), a maize × teosinte population made the same way (backcross, then selfing), also reported more heterozygosity than expected, as did our lab's earlier lines.
- **F1 families:** only 2 of 98 F1 families reach the expected 12.5% donor alleles (Zx.0390_P2, Zx.0130_P1); most sit between 5% and 11%. AB/(AB+BB) is above 22% in most F1 families of every taxon.
- **Small taxa:** Zd (6 F1 families), Zl (5) and Zh (2) have few F1 families, so their intervals are less reliable. Zd's SEG (5.1%) is somewhat below the expectation (p = 0.0004).

**One sample, as an example** (PN12_SID1106, Zv). The state of each 100-kb window along the 10 chromosomes:

![Example sample](final_script/figures/results_spread/example_PN12_SID1106.png)

On such a plot there are only three colours, the three pool states. AB isn't drawn, because a pool can't show which plants are AB: inside a yellow (SEG) region about half the 6 plants are AB, a quarter AA and a quarter BB. All samples can be browsed in the [interactive viewer](#212-interactive-figures-step-23).

The full statistics, F1 family results and the method history are in [`final_script/README.md`](final_script/README.md).

### 2.9 Validation

- **B73 controls are clean:** all 10 have ≤ 0.56% non-B73 windows (donor share ≤ 0.28%). This is a fair test, because the k-mer route never used the controls to filter sites.
- **Known tracts:** they are found with the expected levels. For example, a Zd DONOR tract on chr2 has r ≈ 0.20 against an assembly value of 0.26.
- **Agreement with the SNP route:** donor shares agree where the SNP route works. Where they differ, the cause is on the SNP side: lost donor reads, κ collapsing, or noisy pools.
- **Independent k-mer check (section 2.10):** teosinte-specific k-mers in the reads give the pool's donor taxon (matches the pedigree in 93.9% of pools, 99.7% when the assignment is clear) and its donor share, which agrees with the HMM's across pools (Spearman ρ 0.91–0.98 within taxa).

### 2.10 Population structure from k-mers (step 32)

The HMM calls say *whether* a window carries teosinte, but not *which* teosinte, because they are built only from B73 k-mers that are present or missing. A PCA of the ancestry calls therefore does not separate the five taxa (taxon explains ≤ 2.5% of the variance of each of the first 10 PCs; the PCs follow tracts shared within families).

To see which teosinte each pool carries, we used the k-mers that the donor DNA brings in and B73 does not have, which differ between teosinte taxa. Following the idea of reading population structure from k-mer content without alignment (Hrytsenko et al. 2025), step 32 uses **sourmash** (Brown & Irber 2016; Irber et al. 2022) to compare every pool with the teosinte assemblies:
1. **Sketch** (32a): k = 31, FracMinHash with scaled = 1,000 (every k-mer whose hash is in the lowest 1/1,000 of the hash space; the same k-mers are kept in every sample, so sketches compare directly). One sketch for B73, one per assembly (Zd Gigi, Zd Momo, Zh RIMHU001, Zl RIL003, Zv TIL01, Zv TIL11, Zx TIL18, Zx TIL25), one per pool from its trimmed reads.
2. **Taxon-specific k-mers** (32b): for each assembly, its k-mers minus B73's; for each taxon, the union over its assemblies minus the k-mers of the other four taxa.
3. **Profile per pool** (32b): the share of each taxon's specific k-mers found in the pool's reads, f_taxon. Its largest value is the pool's **k-mer taxon**, which is checked against the pedigree taxon. The share of all teosinte k-mers found, f_teo, measures how much donor DNA the pool has, independently of the HMM.
4. **PCA** (32c): each pool's f values divided by their sum (depth and donor amount cancel; only *which* teosinte is left), centred log-ratio, `prcomp`. Points coloured by pedigree taxon; arrows = taxa.

The results below are for the 1,385 teosinte donor pools; the 29 landrace pools are left out.

![PCA of the teosinte k-mer profiles](data/introgression_spread/sourmash/pca/kmer_pca.png)

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
  > This fits the known gene flow between the two subspecies (Aguirre-Liguori et al. 2019): higher-elevation *parviglumis* is closer to the elevation range of *mexicana*. They are also close to *mexicana* on the map: 23, 25 and 45 km from the nearest Chalco collection site, among the five closest of the 36 *parviglumis* accessions (median 97 km). A labelling error would more likely affect single pools than most pools of an accession. These three are worth checking against the seed records. In the k-mer HMM, these pools also fit a lower DONOR level than the Zv assemblies give (scale factor a: median 0.8, against 0.9 for the other Zv pools), as expected if their donor DNA is partly *mexicana*-like (Zx keeps fewer B73 k-mers than Zv: d = 0.37 vs 0.43).
  >
  > One other explanation should be ruled out first. The taxon-specific k-mers come from only two assemblies per taxon. If the two *parviglumis* assemblies (TIL01, TIL11) were collected far from these accessions, these accessions would share fewer of their specific k-mers and look relatively more *mexicana*-like even without gene flow. Where TIL01, TIL11, TIL18 and TIL25 were collected should be checked.
- **In *diploperennis*,** the mismatches come mostly from Zd.0030 and Zd.0040 (both from Milpillas, 1,390 m, listed as race "Zea-diploperennis-Na"): their assignments are weaker (median ratio 1.8 and 1.5, against 3.0 and 3.5 for Zd.0010 and Zd.0020), so they are less like the two Zd assemblies.

![Per accession: pools whose k-mer taxon differs from the pedigree taxon, by elevation](data/introgression_spread/sourmash/pca/kmer_mismatch_elevation.png)

Next, we asked whether this *mexicana*-like signal follows elevation across all *parviglumis* accessions, not only in the pools that were assigned to the wrong taxon. For this, each pool got a ***mexicana* index**: of the taxon-specific k-mers it shares with either *mexicana* or *parviglumis*, the fraction that are *mexicana*, f_Zx / (f_Zx + f_Zv). An index near 0 means the pool looks like *parviglumis*, near 1 like *mexicana*. We took the median over the pools of each accession and plotted it against the accession's elevation (right panel below; the left and middle panels are the PCA coloured by elevation).
- **Within *parviglumis*, the index rises with elevation** (Spearman ρ = 0.49, 36 accessions, p = 0.003), and the highest accessions reach the values of *mexicana*. So the two subspecies form a continuous gradient, not two separate groups.
- **The index also rises from west to east** (with longitude, ρ = 0.73). Elevation still adds to it once longitude is accounted for (partial ρ = 0.45, p = 0.006), but with 36 accessions, geography and elevation cannot be fully separated.
- **Within *mexicana*, the index does not change with elevation** (ρ = 0.12, 37 accessions).
- **How much teosinte DNA a pool carries does not depend on its donor's elevation** (HMM donor share, ρ = 0.00), so the gradient is about *which* teosinte, not how much.

![k-mer PCA coloured by elevation, and the mexicana index along elevation](data/introgression_spread/sourmash/pca/kmer_pca_elevation.png)

We then looked within *mexicana*, at its four races: Chalco, Mesa Central, Durango and Nobogame (section 2.1). The taxon-level profile cannot separate them, since all four are *mexicana*. Instead, we used each pool's share of the k-mers of each of the 8 assemblies separately, and compared the two *mexicana* assemblies, TIL18 and TIL25 (figure below; collection sites from `data/metadata_BZea.csv`, where 75 of the 126 teosinte longitudes lack the minus sign; the scripts correct this, since all sites are in Mexico and Guatemala):
- **A, collection sites:** Chalco comes from the Valley of Mexico in the southeast of the range, Mesa Central from the central plateau to its west, and Durango and Nobogame from the north.
- **B, PCA of the 577 *mexicana* pools:** Chalco separates from the other races along PC2, and Mesa Central lies between Chalco and the northern races.
- **C, TIL18 vs TIL25:** Chalco pools share relatively more k-mers with TIL18 (median log ratio +0.20), while Mesa Central (−0.07), Durango (−0.18) and Nobogame (−0.16) are closer to TIL25.
- **D, *mexicana* index along elevation:** the *mexicana* races and the *parviglumis* accessions on one axis.

So Chalco differs from the other three races in its k-mer content, as expected from geography. It is also the race that carries the late-flowering allele at chr5 221.8 Mb (section 3.6.13). Which collections TIL18 and TIL25 come from should be checked before reading more into the direction of panel C.

![The four mexicana races](data/introgression_spread/sourmash/pca/mexicana_races.png)

Finally, the k-mer profile gives an independent check of the HMM's donor share. The HMM measures the donor share from *B73* k-mers that are missing; sourmash measures it from *teosinte* k-mers that are present (f_teo). They use different k-mers from the same reads. After correcting f_teo for sequencing depth (residual of log f_teo on log sketch size, within taxon), the two agree closely across pools: Spearman ρ = 0.91 (Zd), 0.98 (Zh), 0.97 (Zl), 0.94 (Zv) and 0.97 (Zx) (1,292 pools, NOISY pools left out).

![Teosinte k-mers vs the HMM donor share](data/introgression_spread/sourmash/pca/kmer_vs_hmm.png)

Per-pool shares of each taxon's k-mers: `data/introgression_spread/sourmash/pca/kmer_taxon_shares.png`. Tables: `kmer_taxon_profile.tsv` (per pool), `pca/kmer_taxon_agreement.tsv`, `pca/kmer_accession_agreement.tsv` (per accession, with elevation), `pca/kmer_vs_hmm.tsv`, `pca/kmer_pca_scores.tsv`.

```bash
N=$(awk -F'\t' 'NR>1 && ($2=="control"||$2=="donor")' ${D11}/pools.tsv | wc -l)
J=$(sbatch --parsable --array=1-$(( 9 + (N + 9) / 10 ))%60 32a_sourmash_sketch.sh)   # needs: pip install sourmash
sbatch --dependency=afterok:$J 32b_kmer_taxon_profile.sh
```
Output in `${D32}` (copied to `data/introgression_spread/sourmash/`): `kmer_taxon_profile.tsv` (per pool: f per taxon, f_teo, containment per assembly, k-mer taxon, agreement) and `pca/`. The PCA step alone, on a laptop: `POOL_SUMMARY=<pool_summary.tsv> Rscript final_script/32c_kmer_pca.R kmer_taxon_profile.tsv pca data/metadata_BZea.csv` (needs R with data.table, ggplot2, patchwork).

### 2.11 Noisy pools and resequencing

93 of the 1,385 donor pools are flagged **NOISY** (σ > 0.15). They are left out of the results above. Almost all of them are on three plates with low sequencing depth:

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
| **PN11** | **93** | **34** | **0.130** | **150 M** | **Zl 18, Zv 16** |
| PN12 | 95 | 4 | 0.095 | 211 M | Zv 4 |
| **PN13** | **88** | **23** | **0.140** | **115 M** | **Zd 18, Zv 5** |
| PN14 | 94 | 1 | 0.090 | 207 M | Zd 1 |
| **PN15** | **93** | **24** | **0.138** | **109 M** | **Zv 22, Zd 2** |
| PN16 | 95 | 0 | 0.084 | 214 M | — |
| PN17 | 96 | 0 | 0.082 | 217 M | — |

**Recommendation: resequence PN11, PN13 and PN15.**
- **Depth:** PN13 and PN15 have about half the k-mers of a typical plate (109–115 M vs ~200 M), and PN11 about 75%. On PN11, 37% of pools are NOISY.
- **Taxa:** the noisy pools are concentrated in a few taxa, so leaving them out takes a large bite out of those taxa:
  - Zl: 18 of its 20 noisy pools are on PN11.
  - Zd: 18 of 22 are on PN13.
  - Zv: 43 of 47 are on PN11, PN13 and PN15.

**Other samples to check.** Besides depth, a few samples look wrong for other reasons. They are kept in the results but flagged; together they change the main numbers by less than 0.1 percentage point.
- **Donor pools with almost no donor.** 5 pools are as clean as the B73 controls (donor share ≤ 0.11%): PN4_SID334 and PN4_SID335 (both from BC2 plant `Zx.0020_P3_P1_P2`), PN14_SID1313 and PN14_SID1314 (both from `Zd.0010_P2_P5_P4`), and PN7_SID649. Losing every donor segment by chance is very unlikely, and the two pairs each come from one BC2 plant. So these BC2 plants may have been B73 (an outcross or a seed mix-up), or the samples may be mislabelled B73 or check plants. Worth checking against the field and seed records.
- **Sample sheet.** PN13_SID1225 is listed twice, as `Zd.0020_P3_P1_P2.5.1.1` and as PurpleCheck; it is left out. The sheet also has 21 PurpleCheck and 13 B73 samples. The checks are not donor lines and are left out. 10 of the B73 samples are the controls; the other 3 (PN1_SID85, PN3_SID236, PN5_SID468) are not in the k-mer run and should be checked.

The full statistics and the per-plate table are in [`final_script/README.md`](final_script/README.md).

### 2.12 Interactive figures (step 23)

[`final_script/figures/introgression_viewer_spread.html`](final_script/figures/introgression_viewer_spread.html) is one self-contained page with every pool's chromosome plot ([download it](https://raw.githubusercontent.com/nirwan1265/BZea_genotyping/main/final_script/figures/introgression_viewer_spread.html) and open it in a browser).
- **Drop-down menus:** taxon → F1 family → "All pools" or one pool. There is also a search box for a pool ID.
- **Family view:** each chromosome shows one track per pool of the F1 family, so the segments the pools share line up.
- **Pool view:** the 10 chromosomes of one pool.
- **Hover** a segment for its position, length and state. A table under the figure lists the F1 family's pools (σ, flag, donor share, SEG, DONOR, tracts).
- **Download SVG / PNG** saves the figure on screen, legend included, for papers and slides.
- B73 is blue, SEG yellow and DONOR red, the same colors as in `blocks_<GROUP>.pdf`.

To view it, download the file and open it in a browser. GitHub shows HTML files as code, not as a page.

### 2.13 Limitations and future work

**Short tracts aren't proven yet.** The k-mer HMM calls tracts down to about 400 kb. Short tracts are expected from the crossing scheme ([Tract lengths](#265-how-long-are-the-tracts)), so most are probably real, but single short calls haven't been checked with independent data.

**SEG vs DONOR isn't checked directly either.** A SEG call says the family still segregates (plants AA, AB or BB). No independent data has confirmed this split yet. The split is what gives the main result (het as expected, homozygous donor about half), so it should be checked independently.

**The DONOR level depends on the assemblies.** The model needs to know how much of each window's B73 k-mers a teosinte genome keeps (d_w). The best source would have been **the actual donor parents** of our lines, sequenced. We don't have them, so we use teosinte assemblies from MaizeGDB (PanAnd). Those are other plants of the same taxa, so in diverse taxa they can differ from our donors:
- Zv and Zx still switch between SEG and DONOR about twice as often as the other taxa, even after allowing for this uncertainty.
- Zh and Zl have one assembly each, so the uncertainty can't be measured for them.
- Zm (landraces) can't be done yet, because it needs a maize-based donor level (e.g. from NAM founder assemblies).

A first check of whether the short tracts are real uses **relatives**. A real introgression is inherited: pools that descend from the same plant received their donor segments from that plant, so a real tract in one pool should often be found at the same place in its relatives. A false tract (noise) appears at random places, so it should be found in relatives no more often than in unrelated pools.

The relatives are read from the pedigree name. For example, in `Zv.0200_P3_P2_P1.1.1.1`:
- `Zv.0200` is the accession (the teosinte donor);
- `P3` is the F1 plant, `P2` the BC1 plant and `P1` the BC2 plant;
- `.1.1.1` are the three rounds of selfing that led to this plot.

So `Zv.0200_P3_P2_P1.1.1.1`, `Zv.0200_P3_P2_P1.2.1.1` and `Zv.0200_P3_P2_P1.3.1.1` all descend from the same BC2 plant, `Zv.0200_P3_P2_P1`, through different selfed lines.

The check works like this (1,292 pools; SEG and DONOR counted together as "introgression"):
1. Take one tract in one pool.
2. Take a group of other pools with a given relationship to it, e.g. all pools from the same BC2 plant.
3. Count how many of those pools also have an introgression at the position of the tract.
4. Repeat for every tract in every pool, separately for short (< 10 Mb) and long (≥ 10 Mb) tracts.

The table gives the result as a percentage: for example, 31.4% means that for a short tract in one pool, on average 31.4% of the other pools from the same BC2 plant have an introgression at the same place.

| The other pools | Short tracts (< 10 Mb) | Long tracts (≥ 10 Mb) |
|---|---|---|
| Same taxon, other accession (unrelated: the baseline) | 17.5% | 23.7% |
| Same accession, other F1 plant | 17.2% | 22.6% |
| Same F1 plant, other BC1 plant | 17.1% | 22.6% |
| Same BC1 plant, other BC2 plant | 17.4% | 24.4% |
| **Same BC2 plant** | **31.4%** | **43.8%** |

- **The baseline is not zero.** Introgressions are common and some genome regions carry them in many pools, so even unrelated pools share some by chance. That's why every row is compared with the baseline, not with zero.
- **Pools from the same BC2 plant share tracts about 1.8 times as often as unrelated pools** (31.4% vs 17.5% for short tracts, 43.8% vs 23.7% for long tracts).
- **Short tracts rise by the same factor as long tracts.** If the short tracts were mostly noise, they would not be inherited, and their sharing would stay at the baseline even among close relatives. Since they rise like the long tracts, most short tracts are probably real.
- **Pools from the same accession but a different BC2 plant do not share more than the baseline.** If the short tracts were an artifact of the accession's donor plant being different from the teosinte assemblies, they would appear at the same places in every pool of that accession. They don't, so that artifact is not what makes them.
- **Open question:** pools from the same BC1 plant but different BC2 plants should also share more than the baseline, about half the increase seen for the same BC2 plant, because each BC2 plant inherits about half of its BC1 parent's donor segments. They don't (24.4% vs 23.7% for long tracts). Why is not known yet.

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
4. **Donor plants.** If seed or tissue of the donor parents (or of the same accessions) can be found, 10–20x short-read sequencing gives their own d_w directly. No assembly is needed.
5. **Noisy plates.** Resequence PN11, PN13 and PN15.

---

## Section 3 — k-mer ancestry GWAS and joint QTL mapping

We ran **two analyses** on the same genotypes (k-mer teosinte ancestry per 100-kb window, Section 2) and the same phenotypes (two-year SpATS BLUEs):

| | **Analysis 1: k-mer ancestry GWAS** (step 26, GridLMM) | **Analysis 2: joint QTL model** (step 29.2, R/qtl2) |
|---|---|---|
| Question | Where in the genome are the QTL? | How many separate QTL, how large is each, and does the effect differ between teosinte taxa? |
| What is tested | Every window, one at a time | The QTL found by Analysis 1, plus any second QTL that a conditional scan finds inside one of its regions; all fitted together |
| Effect | One effect per teosinte copy, shared by all taxa | The effect of each QTL adjusted for the others, and its effect within each taxon |
| Significance | 5% genome-wide threshold from permutations | Bonferroni over the QTL of that trait (no genome search, so no permutations needed) |
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
| SPAD2 (leaf greenness, SPAD meter, 36 days after silking); 2023 only | 0.40 | | |
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

The table compares the two choices of K on the same data (DTS, 1,284 plots, taxon fixed), using **λ_GC** (the genomic inflation factor; Devlin & Roeder 1999). λ_GC is the median test statistic over all windows divided by the median expected if no window had any effect. λ_GC = 1 means the test statistics look as expected under no effect; above 1, they are larger than expected (inflated); below 1, smaller (deflated, the model is too conservative). In most GWAS a λ_GC well above 1 is read as a warning of uncorrected population structure.

| K | λ_GC | Best window chr10 (ZmCCT10 region, 94–98 Mb) | Best window chr5 (~221.7 Mb) |
|---|---|---|---|
| **LOCO** | 3.97 | −log10 p **9.1**, +1.0 day per teosinte copy | **9.7** |
| Whole genome | 0.65 | 2.9 | 6.1 |

With the whole-genome K, the chr10 flowering signal almost disappears and λ_GC falls **below 1**: the model over-corrects. LOCO is the fix.

**λ_GC under LOCO (~4) is not a sign of false positives here.** λ_GC assumes most markers are null and independent. In this population each QTL is tagged by hundreds of linked windows (tracts of several Mb, and with LOCO the tested chromosome's own background is not in K), and the traits are polygenic, so the median test statistic is expected to be high (Yang et al. 2011). Significance is therefore set by **permutation**, not by λ_GC.

#### 3.3.3 Model

GridLMM (Runcie & Crawford 2019), maximum likelihood, one chromosome at a time with its LOCO K:

`y = taxon (fixed) + b · dosage_w + u + e`,  `u ~ N(0, σ²_g K_LOCO)`,  `e ~ N(0, σ²_e I)`

- b = effect per teosinte copy (trait units, e.g. days); likelihood-ratio test for b (1 df).
- Dominance: the same model plus c · SEG_w; likelihood-ratio test for c.
- GridLMM fits the variance proportion on a grid instead of re-optimising it for every window, which makes ~21,000 tests with a 1,284 × 1,284 K fast (~3 min per trait on 4 cores). The code started from `scripts/Analyze_QTL_Founder2.R`.

We used **taxon** (5 levels: Zd, Zh, Zl, Zv, Zx) as the fixed effect, not the **F1 family**. An F1 family is all the plots that descend from one F1 plant (one teosinte plant × B73); the 1,284 plots come from 95 F1 families of 82 accessions. Plots of the same F1 family share many of their tracts (relatives check, Section 2, Limitations). With F1 family (95 levels) as a fixed effect, a tract carried by a whole family is confounded with the family mean and absorbed, and 94 extra parameters raise the threshold. Family relatedness is already in K. For DTS (829 plots, one year): the peaks on chr7, chr8/9 and chr10 were similar either way, but the 5% threshold rose from 3.8 to 5.0 and the chr3 peak (4.4) was lost. F1 family as a fixed effect is kept as a check (`FIXED_COL=F1FAMILY`).

#### 3.3.4 Significance: permutations

For each trait, phenotypes are shuffled **within taxon** (keeping the taxon means), the whole genome scan is rerun, and the largest −log10 p in the genome is kept. The 95th percentile of these 1,000 genome-wide maxima is the threshold. It controls the **family-wise error rate (FWER) at 5% across the whole genome**: if the trait had no QTL, a scan would put any window above it in only 5% of cases (Churchill & Doerge 1994). Because each permutation keeps the real correlation between neighbouring windows, the threshold accounts for it; Bonferroni over 21,309 windows (−log10 p = 5.63) would treat them as independent and be too strict. The thresholds come from **1,000 permutations** per trait. With *mexicana* split by race (section 3.6.13), phenotypes are shuffled within race.

#### 3.3.5 Regions and genes

Step 27 joins windows above the threshold into regions (gap < 5 Mb), keeps one region per locus, and gives each a support interval: the windows within 1.5 −log10 p of the top window (like a 1.5-LOD drop). Every B73 v5 gene overlapping the interval is listed from the GFF3. Output: `genes/regions_genes_<TRAIT>_ALL.tsv` (one row per gene, with the region's position, top −log10 p and effect). Step 30 (`final_script/30_genes_go.R`) adds each gene's window −log10 p and effect and its GO terms: `data/introgression_spread/GO/genes_GO_all_<TRAIT>.tsv`.

### 3.4 Analysis 2: joint QTL model with per-taxon effects

The scan in Analysis 1 tests each window on its own. A peak can therefore be a real QTL, or a "shadow" of a stronger QTL elsewhere on the chromosome that it shares tracts with. Analysis 2 takes the QTL from the scan and tests them **together**, to answer three questions for each one: is it still significant when the other QTL are in the model, how large is its effect, and is the effect the same in all teosinte taxa? It uses the package R/qtl2 (Broman et al. 2019), with the same plots, phenotypes, taxon and LOCO kinship as Analysis 1 (`final_script/29_joint_qtl.R`).

**Input QTL.** For each region of the scan, the top window and its support interval (`peaks_<TRAIT>_ALL.tsv`, step 27).

**Genotypes.** At each window, a plot's teosinte allele frequency p = dosage / 2 (B73 0, SEG 0.5, DONOR 1). R/qtl2 takes this as the probability of the teosinte allele (B73 allele: 1 − p). The kinship is LOCO, from the same values (`calc_kinship(probs, "loco")`).

**Step 1: more than one QTL in a region? (conditional scans)** A region of the scan is a run of windows above the threshold. With long tracts, one region can span tens of Mb and could hold two QTL. To check, the whole genome is scanned again with taxon and all QTL found so far as covariates:

`y = taxon + Σ_i γ_i · p_i + β · p_w + u + e` (for every window w; i = the QTL already in the model)

If the best window is still above the genome-wide threshold of Analysis 1, it is a QTL that the others do not explain: it is added to the list, and the scan is repeated, until no window passes. Two rules stop it from adding the wrong things: a window must be significant in the ordinary scan too (so this can split a region, but never create a new one), and it must not be in strong LD with a QTL already in the model (r² ≤ 0.5), since two nearly identical predictors give unstable, opposite effects. This is forward selection by conditional scans (Broman & Sen 2009). Close peaks with no independent signal stay one QTL. Added QTL are listed in `qtl2_<TRAIT>_added.tsv` and drawn as open triangles in the figures.

**Step 2: joint test and effect of each QTL.** Each QTL j is then fitted with taxon and **all the other QTL** as covariates (`fit1`):

`y = taxon + Σ_{i≠j} γ_i · p_i + μ + a_j · (2 p_j − 1) + u + e`,  `u ~ N(0, σ²_g K_LOCO)`,  `e ~ N(0, σ²_e I)`

- a_j is the **effect per teosinte copy**: going from B73 (p = 0) to DONOR (p = 1) changes the fitted value by 2a_j, i.e. a_j per copy (set with `contrasts = cbind(mu = c(1, 1), a = c(-1, 1))`), with its standard error.
- The **joint test** compares this model with the same model without QTL j: LOD = log10 of the likelihood ratio, converted to −log10 p with 1 df (χ² = 2 ln(10) · LOD).
- If QTL j's signal disappears once the others are in the model, it was a shadow of another QTL, not a separate one.

**Step 3: effect within each taxon.** The same model as step 2 (without the taxon term), fitted only to the plots of one taxon, gives that taxon's effect per copy. It is estimated only where the taxon has ≥ 5 carrier plots (plots with SEG or DONOR at the QTL).

**Step 4: do the taxa differ?** At the QTL window, the model with a separate effect for each taxon (`scan1(..., intcovar = taxon)`)

`y = taxon + μ + Σ_t a_t · 1[taxon = t] · (2 p − 1) + u + e`

is compared with the model with one shared effect a. The LOD difference is converted to −log10 p with (number of taxa − 1) df. A significant result means the teosinte allele at that QTL has a different effect depending on the taxon it came from.

**Significance.** The joint model does not search the genome: for a trait with k QTL it makes k joint tests and k "taxa differ" tests. Each is significant if −log10 p > −log10(0.05 / k) (Bonferroni over the k QTL of that trait; e.g. 2.08 for 6 QTL). This also controls the family-wise error rate at 5%, here over the k QTL.

Every statistical step is an R/qtl2 function; the script only converts the k-mer calls into R/qtl2's format.

#### 3.4.1 Why the joint model starts from the peaks, and not from all windows at once

A joint model needs a short list of QTL to fit together. It cannot take all 21,309 windows: that is far more predictors than the 1,284 plots, and neighbouring windows are nearly identical (one tract covers dozens of them), so their effects cannot be separated. The genome scan is the step that picks the list: one window per peak (the top window), with its support interval. The joint model then asks, for those few QTL together, which ones still matter once the others are in the model, how large each is, and in which taxa. This two-step design (scan, then fit the selected QTL jointly) is how QTL mapping is usually done; TeoNAM's joint linkage mapping chooses its QTL by stepwise selection instead (Chen et al. 2019). Using the scan's peaks keeps the two analyses on exactly the same QTL.

#### 3.4.2 Why some QTL show fewer taxa

A per-taxon effect needs carriers: plots of that taxon with the teosinte allele at the QTL. Taxa with fewer than 5 carrier plots at a QTL get no estimate. For example, DTS Q2 (chr10) has carriers from Zd (23), Zv (22) and Zx (62) but only 4 from *huehuetenangensis* and 2 from *luxurians*, so those two are not shown. Which taxa carry a tract depends on where the backcrosses happened to keep donor segments.

This is the mixed-model counterpart of TeoNAM's joint linkage mapping (Chen et al. 2019). Effects are estimated per taxon, not per accession: a typical window has ~90 carrier plots across 82 accessions, about one per accession. Details: `final_script/README.md`, section 12.

### 3.5 How to run

```bash
# Cluster: everything at once (inputs, scan, 1,000 permutations, peaks + genes, joint QTL),
# for taxon and for mexicana split by race; from final_script/, with the GWAS conda env active
bash submit_assoc_all.sh
# Laptop alternative for Analysis 1: one terminal per trait (scan + first 100 permutations), then up to 1,000
bash final_script/run_assoc_local.sh DTS
TOTAL=1000 bash final_script/run_perm_local.sh DTS
# Check: the same analysis on each year's BLUEs separately (default DTS; folders *_y2023, *_y2025)
bash submit_year_check.sh
# genes + GO, and the side-by-side figures (laptop)
D=data/introgression_spread
Rscript final_script/30_genes_go.R $D/assoc_gridlmm_spats/genes $D/GO/GO_terms_all_candidates.txt $D/GO
Rscript final_script/31_joint_qtl_figures.R $D/assoc_gridlmm_spats $D/joint_qtl/noperm $D/joint_qtl/noperm/figures
```
All results are in `data/introgression_spread/`:
- `assoc_gridlmm_spats/`: scans, summaries, 1,000 permutations (`perm_<TRAIT>_seed<N>.tsv`), Manhattan plots, `genes/` (peaks and every gene in each region)
- `joint_qtl/noperm/`: joint QTL tables, per-taxon effects, genes per QTL, and the figures
- `GO/`: every gene in the regions with GO terms and its window's −log10 p (`genes_GO_all_<TRAIT>.tsv`)
- `assoc_gridlmm_spats_subtaxon/`, `joint_qtl/noperm_subtaxon/`, `GO_subtaxon/`: the same with *mexicana* split by race

### 3.6 Results

**Overview (1,196–1,286 plots per trait; thresholds from 1,000 permutations):**

| Trait | Plots | 5% threshold (−log10 p) | QTL | Hold in joint model | Strongest QTL |
|---|---|---|---|---|---|
| DTS | 1284 | 5.25 | 2 | 2 | chr5 221.8 Mb (−log10 p 9.7, +1.00 d per copy) |
| DTA | 1284 | 5.14 | 3 | 3 | chr10 94.5 Mb (−log10 p 9.3, +0.92 d per copy) |
| ASI | 1284 | 4.35 | 12 | 6 | chr7 165.6 Mb (−log10 p 13.7, -0.36 d per copy) |
| PH | 1283 | 4.19 | 16 | 8 | chr10 97.0 Mb (−log10 p 10.8, +5.64 cm per copy) |
| EH | 1282 | 4.27 | 18 | 7 | chr4 198.4 Mb (−log10 p 7.6, +3.00 cm per copy) |
| EN | 1282 | 4.57 | 20 | 11 | chr2 22.4 Mb (−log10 p 18.6, +0.17 per copy) |
| Prolif | 1286 | 4.55 | 15 | 10 | chr1 23.2 Mb (−log10 p 20.3, +0.20 per copy) |
| NBR | 1274 | 4.21 | 8 | 4 | chr10 98.2 Mb (−log10 p 10.3, +0.17 per copy) |
| SL | 1196 | 4.06 | 13 | 12 | chr7 141.7 Mb (−log10 p 9.2, -0.32 per copy) |
| BL | 1197 | 4.14 | 6 | 5 | chr2 39.8 Mb (−log10 p 12.0, -1.58 per copy) |
| BW | 1197 | 4.22 | 7 | 6 | chr6 146.4 Mb (−log10 p 9.2, -0.29 per copy) |

"Hold in joint model" = QTL whose joint −log10 p passes Bonferroni over that trait's QTL. The other six traits are in [3.6.14](#3614-other-traits).

**Reading the figures.** Left: the k-mer ancestry GWAS; −log10 p per 100-kb window, dotted line = 5% permutation threshold, QTL numbered Q1, Q2 … by position; an open triangle marks a QTL added by the conditional scans (a second QTL inside a significant region, section 3.4). Middle: the same QTL in the joint model; effect per teosinte copy (trait units) with 95% CI, black = joint model across all taxa, colours = within each taxon (taxa with ≥ 5 carrier plots). Right: two columns of **−log10 p** (not effects): "Joint test" = is the QTL still significant with all other QTL in the model; "Taxa differ" = does its effect differ between taxa. For example, DTS Q1 "9.7 | 14.1" means p = 10⁻⁹·⁷ for the QTL and p = 10⁻¹⁴·¹ for the difference between taxa. Positive effect = the teosinte allele increases the trait (for DTS/DTA: later flowering).

#### 3.6.1 Days to silking (DTS)

![DTS: k-mer ancestry GWAS and joint QTL model](figures/kmerGWAS_joint_qtl_DTS.png)

Two QTL pass the 1,000-permutation threshold (−log10 p > 5.25), and both hold in the joint model (Bonferroni over 2 QTL: −log10 p > 1.60).

| QTL | Known genes | Effect (joint) | Joint −log10 p | Taxa differ | What it says |
|---|---|---|---|---|---|
| **Q1 chr5 221.8 Mb** | Zm00001eb257680 (HDR1-like) | +1.00 d | 9.7 | 14.1 (yes) | Differs between taxa (Zd +1.94, Zl +4.77, Zv +0.46, Zx +0.48). Strong, but mostly *luxurians* and *diploperennis*; Zv and Zx only about +0.5 d. |
| **Q2 chr10 94.5 Mb** | Zm00001eb418700 (ZmCCT10) | +1.01 d | 8.9 | 1.0 (no) | No significant difference between taxa (Zd +1.19, Zv +0.84, Zx +0.94). The ZmCCT10 (Zm00001eb418700) region: teosinte alleles delay flowering everywhere. |

**Suggestive peaks** (above −log10 p 4.5 but below the 5.25 threshold): chr1 85.8 Mb (5.1; mostly *luxurians*), chr9 10.9 Mb (4.9), chr5 206.4 Mb (4.9; holds Zm00001eb250710, *ZAG3*/*bde1*) and chr9 41.4 Mb (4.6). chr9 10.9 Mb is significant for DTA (below).

**Candidate genes:**
- **Q2, chr10:** *ZmCCT10* (**Zm00001eb418700**, annotated "GHD7_2"; 96.18 Mb), the main teosinte photoperiod gene; teosinte alleles delay flowering under long days (Hung et al. 2012; Yang et al. 2013).
- **Q1, chr5:** **Zm00001eb257680**, annotated as **Protein HEADING DATE REPRESSOR 1** (SnRK1-interacting protein; GO "regulation of flower development"), in the top window (221.7–221.8 Mb, 8 genes). Its functional ortholog in rice, *HDR1* (*Heading Date Repressor 1*), is a nuclear protein that suppresses flowering under long days: it upregulates *Hd1* and downregulates *Ehd1*, which lowers the florigens *Hd3a* and *RFT1*; it interacts with the kinase OsK4, another long-day repressor (Sun et al. 2016). The *hdr1* mutant flowers early in long days only, i.e. HDR1 acts through the photoperiod pathway, which fits a teosinte allele that delays silking. Not characterised for flowering in maize. Orthology still to be confirmed with a gene tree. The same window is the DTA Q1 peak.

<details><summary>All DTS QTL with per-taxon effects and all genes in the interval</summary>

| QTL | Position | Interval (Mb) | Scan −log10 p | Effect per copy (joint) | Joint −log10 p | Taxa differ −log10 p | Per-taxon effects (carrier plots) | Genes in interval |
|---|---|---|---|---|---|---|---|---|
| Q1 | chr5 221.8 Mb | 221.7–222.0 | 9.7 | +1.00 d ± 0.16 | **9.7** | **14.1** | Zd +1.94 (10), Zl +4.77 (13), Zv +0.46 (49), Zx +0.48 (50) | [16 genes](qtl_genes/DTS_qtl_genes.tsv) |
| Q2 | chr10 94.5 Mb | 92.8–102.3 | 9.1 | +1.01 d ± 0.17 | **8.9** | 1.0 | Zd +1.19 (23), Zv +0.84 (22), Zx +0.94 (62) | [145 genes](qtl_genes/DTS_qtl_genes.tsv) |

Bonferroni over the 2 QTL: −log10 p > 1.60 (bold). Every gene in each interval, with its GO annotation and the k-mer GWAS −log10 p of its window: [`qtl_genes/DTS_qtl_genes.tsv`](qtl_genes/DTS_qtl_genes.tsv).

</details>

#### 3.6.2 Days to anthesis (DTA)

![DTA: k-mer ancestry GWAS and joint QTL model](figures/kmerGWAS_joint_qtl_DTA.png)

Three QTL, all hold in the joint model (Bonferroni over 3: −log10 p > 1.78).

| QTL | Known genes | Effect (joint) | Joint −log10 p | Taxa differ | What it says |
|---|---|---|---|---|---|
| **Q3 chr10 94.5 Mb** | Zm00001eb418700 (ZmCCT10) | +0.92 d | 9.1 | 1.4 (no) | No significant difference between taxa (Zd +1.27, Zv +0.72, Zx +0.74). Same QTL as DTS Q2 (ZmCCT10, Zm00001eb418700), the strongest DTA QTL; all taxa delay. |
| **Q2 chr9 10.9 Mb** | Zm00001eb372960 (MADS74) | -0.54 d | 5.6 | 0.6 (no) | No significant difference between taxa (Zd -0.69, Zh -0.39, Zl -1.51, Zv -0.34, Zx -0.48). Teosinte flowers earlier. Suggestive for DTS (−log10 p 4.9). |
| **Q1 chr5 221.8 Mb** | Zm00001eb257680 (HDR1-like) | +0.65 d | 5.5 | 11.2 (yes) | Differs between taxa (Zd +1.47, Zl +3.67, Zv +0.21, Zx +0.19). Same QTL as DTS Q1 (Zm00001eb257680), mostly *luxurians* and *diploperennis*. It delays pollen shedding less than silking (DTS +1.00 d), so it also lengthens ASI. |

**Candidate genes:**
- **Q1, chr5 221.8 Mb, shared with DTS Q1:** Zm00001eb257680, the ortholog of rice *HEADING DATE REPRESSOR 1* (see DTS).
- **Q3, chr10 94.5 Mb, shared with DTS Q2:** Zm00001eb418700, *ZmCCT10* (see DTS).
- **Q2, chr9 10.9 Mb (significant for DTA, suggestive for DTS):** Zm00001eb372960 (*MADS74*; 8.36 Mb, inside the support interval 6.1–13.7 Mb). Here the teosinte allele makes plants shed pollen earlier (−0.54 d per copy), alike in every taxon. *MADS74* is a MADS-box transcription factor that PANTHER places in the subfamily of tomato *JOINTLESS*, an SVP-like MADS-box gene. If it acts like Arabidopsis *SVP*, it would be a flowering repressor, and an earlier-flowering teosinte allele would be a weaker repressor:
  - SVP delays flowering by repressing the florigen gene *FT*, binding directly to the CArG motifs in its sequence. It also links flowering to temperature: plants flower later at 16 °C than at 23 °C, but *svp* mutants flower at almost the same time at both (Lee et al. 2007). At 16 °C, SVP switches on a second repressor, *TEMPRANILLO2*, which also represses *FT* and *TWIN SISTER OF FT* (Marín-González et al. 2015).
  - Natural alleles of *SVP* change flowering time: the Arabidopsis flowering QTL *FAQ1* is *SVP*, where one amino-acid change makes a loss-of-function, early-flowering allele, and how large its effect is depends on the genetic background (Méndez-Vigo et al. 2013).
  - Whether *MADS74* is a true *SVP* ortholog, or a more distant member of the family, should be checked in its Gramene gene tree, as was done for *mads73* (ASI Q11).

<details><summary>All DTA QTL with per-taxon effects and all genes in the interval</summary>

| QTL | Position | Interval (Mb) | Scan −log10 p | Effect per copy (joint) | Joint −log10 p | Taxa differ −log10 p | Per-taxon effects (carrier plots) | Genes in interval |
|---|---|---|---|---|---|---|---|---|
| Q1 | chr5 221.8 Mb | 221.7–222.3 | 5.5 | +0.65 d ± 0.14 | **5.5** | **11.2** | Zd +1.47 (10), Zl +3.67 (13), Zv +0.21 (50), Zx +0.19 (50) | [26 genes](qtl_genes/DTA_qtl_genes.tsv) |
| Q2 | chr9 10.9 Mb | 6.1–13.7 | 5.5 | -0.54 d ± 0.11 | **5.6** | 0.6 | Zd -0.69 (25), Zh -0.39 (8), Zl -1.51 (22), Zv -0.34 (57), Zx -0.48 (97) | [203 genes](qtl_genes/DTA_qtl_genes.tsv) |
| Q3 | chr10 94.5 Mb | 92.8–105.9 | 9.3 | +0.92 d ± 0.15 | **9.1** | 1.4 | Zd +1.27 (23), Zv +0.72 (22), Zx +0.74 (62) | [197 genes](qtl_genes/DTA_qtl_genes.tsv) |

Bonferroni over the 3 QTL: −log10 p > 1.78 (bold). Every gene in each interval, with its GO annotation and the k-mer GWAS −log10 p of its window: [`qtl_genes/DTA_qtl_genes.tsv`](qtl_genes/DTA_qtl_genes.tsv).

</details>

#### 3.6.3 Anthesis–silking interval (ASI)

![ASI: k-mer ancestry GWAS and joint QTL model](figures/kmerGWAS_joint_qtl_ASI.png)

ASI = anthesis–silking interval, DTS − DTA. 12 regions, 6 hold in the joint model (Bonferroni over 12: −log10 p > 2.38). Many other chr2 and chr3 peaks drop in the joint model: they share one signal through linked tracts.

| QTL | Known genes | Effect (joint) | Joint −log10 p | Taxa differ | What it says |
|---|---|---|---|---|---|
| **Q12 chr7 165.6 Mb** | Zm00001eb323880 (bhlh90) | -0.36 d | 14.0 | 0.3 (no) | No significant difference between taxa (Zd -0.46, Zl -0.53, Zv -0.35, Zx -0.31). The strongest ASI QTL: teosinte shortens the interval. |
| **Q11 chr5 218.3 Mb** | Zm00001eb255670 (mads73) | +0.35 d | 8.8 | 4.1 (yes) | Differs between taxa (Zl +1.08, Zv +0.24, Zx +0.28). Lengthens ASI, mostly *luxurians*; next to the DTS/DTA chr5 221.8 Mb QTL. |
| **Q1 chr1 4.5 Mb** | Zm00001eb001510 (remorin) | +0.24 d | 6.6 | 0.0 (no) | No significant difference between taxa (Zd +0.19, Zh +0.21, Zl +0.14, Zv +0.26, Zx +0.24). Lengthens ASI. |
| **Q10 chr4 206.2 Mb** | Zm00001eb201580 (PAS domain protein, ADAGIO family) | -0.21 d | 4.7 | 0.3 (no) | No significant difference between taxa (Zd -0.41, Zh -0.38, Zl +0.47, Zv -0.26, Zx -0.15). Shortens ASI. |
| **Q2 chr2 26.1 Mb** | Zm00001eb076990 (CKX8); Zm00001eb077080 (GA2ox2) | -0.16 d | 3.6 | 0.1 (no) | No significant difference between taxa (Zd -0.19, Zh -0.46, Zl -0.20, Zv -0.19, Zx -0.13). Shortens ASI. |
| **Q9 chr4 52.0 Mb** | Zm00001eb175150 (tga1) | -0.24 d | 3.0 | 0.1 (no) | No significant difference between taxa (Zd -0.08, Zh -0.54, Zl -0.67, Zv -0.20, Zx -0.18). Shortens ASI. |

**Candidate genes:**
- **Q12, chr7 165.6 Mb (the strongest ASI QTL):** Zm00001eb323880 (*bhlh90*, bHLH transcription factor 90; 165.74 Mb, in a window with −log10 p 13.2 next to the top window). MaizeGDB annotates it as the ortholog of Arabidopsis *FLOWERING BHLH 4* (Neupane et al. 2025); Gramene lists Arabidopsis *FBH3* (AT1G51140) as its closest homolog and rice *OsFBH1* / *POH1* (Os08g0506700) as its closest rice gene. In Arabidopsis the four FBH proteins bind E-boxes in the *CONSTANS* promoter and activate *CO*: their overexpression causes early flowering regardless of photoperiod, and *CO* levels drop in the *fbh* quadruple mutant (Ito et al. 2012). In rice, POH1 forms a dimer with HBP1, both bind the *Hd1* promoter and activate *Hd1*; their overexpression delays flowering by suppressing *Ehd1*, *Hd3a* and *RFT1* (Yin et al. 2023).
- **Q11, chr5 218.3 Mb:** Zm00001eb255670 (*mads73*, MADS-box transcription factor; 217.89 Mb, in a window with −log10 p 9.4). In its Gramene gene tree, the Arabidopsis gene is *SHORT VEGETATIVE PHASE* (*SVP*, AT2G22540; also *AGL22*), a MADS-box repressor of flowering:
  - SVP represses the florigen gene *FT* by binding the CArG motifs in its sequence. It mediates the response of flowering to ambient temperature: plants flower later at 16 °C than at 23 °C, but *svp* mutants flower at almost the same time at both temperatures (Lee et al. 2007).
  - At low temperature (16 °C), SVP also switches on another floral repressor, *TEMPRANILLO2* (*TEM2*), which represses *FT* and *TWIN SISTER OF FT* (Marín-González et al. 2015).
  - Natural variation in *SVP* changes flowering time: the Arabidopsis flowering QTL *FAQ1* is *SVP*, where a single amino-acid change makes an early-flowering allele, and its effect depends on the genetic background (Méndez-Vigo et al. 2013).
- **Q1, chr1 4.5 Mb:** Zm00001eb001510 (remorin family protein; Zm00001d027403 in B73 v4; 4.45 Mb, in a window with −log10 p 6.3). In a maize panel, Wang et al. (2026) found a chr1 segment (v4 4,418,001–4,420,000 bp) both by GWAS for days to silking and days to tasseling and by a genetic-differentiation scan between subpopulations (Fst 0.53), and listed its seven genes, including this one, among 41 flowering-time candidates (their Table 2). In the same segment they point to a NAC transcription factor (*NAC78*, Zm00001d027395) and a histone H4 gene (Zm00001d027394).
- **Other known genes in the intervals:** Zm00001eb175150 (*tga1*; Q9), Zm00001eb201580 (PAS domain protein of the ADAGIO family, which includes the circadian clock F-box proteins; Q10), Zm00001eb077080 (*GA2ox2*) and Zm00001eb076990 (*CKX8*) (Q2), and Zm00001eb132420 (*LNK2*, a circadian clock gene) in Q6, which does not hold.

The chr3 95.8 Mb peak (−log10 p 5.0) lies inside the chr3 108.3 Mb region (88.5–129.4 Mb), and the conditional scan does not add it.

<details><summary>All ASI QTL with per-taxon effects and all genes in the interval</summary>

| QTL | Position | Interval (Mb) | Scan −log10 p | Effect per copy (joint) | Joint −log10 p | Taxa differ −log10 p | Per-taxon effects (carrier plots) | Genes in interval |
|---|---|---|---|---|---|---|---|---|
| Q1 | chr1 4.5 Mb | 3.7–5.4 | 6.4 | +0.24 d ± 0.05 | **6.6** | 0.0 | Zd +0.19 (15), Zh +0.21 (10), Zl +0.14 (12), Zv +0.26 (42), Zx +0.24 (68) | [73 genes](qtl_genes/ASI_qtl_genes.tsv) |
| Q2 | chr2 26.1 Mb | 24.8–29.4 | 7.4 | -0.16 d ± 0.04 | **3.6** | 0.1 | Zd -0.19 (35), Zh -0.46 (5), Zl -0.20 (21), Zv -0.19 (55), Zx -0.13 (105) | [147 genes](qtl_genes/ASI_qtl_genes.tsv) |
| Q3 | chr2 68.2 Mb | 60.7–77.1 | 5.3 | -0.12 d ± 0.07 | 1.0 | 0.0 | Zd -0.28 (33), Zh -0.14 (7), Zl -0.12 (15), Zv -0.05 (50), Zx -0.09 (87) | [172 genes](qtl_genes/ASI_qtl_genes.tsv) |
| Q4 | chr2 146.7 Mb | 100.9–150.2 | 4.7 | +0.04 d ± 0.07 | 0.2 | 0.0 | Zd +0.17 (24), Zh -0.14 (6), Zl -0.18 (14), Zv -0.02 (55), Zx +0.04 (111) | [547 genes](qtl_genes/ASI_qtl_genes.tsv) |
| Q5 | chr2 191.2 Mb | 185.7–193.2 | 5.2 | -0.14 d ± 0.05 | 2.2 | 0.2 | Zd -0.11 (17), Zh -0.22 (6), Zl -0.38 (16), Zv -0.16 (51), Zx -0.13 (112) | [179 genes](qtl_genes/ASI_qtl_genes.tsv) |
| Q6 | chr3 69.5 Mb | 61.9–91.1 | 5.0 | -0.06 d ± 0.13 | 0.2 | 0.5 | Zd -0.42 (25), Zh -0.19 (13), Zl -0.50 (7), Zv -0.30 (29), Zx -0.01 (68) | [159 genes](qtl_genes/ASI_qtl_genes.tsv) |
| Q7 | chr3 108.3 Mb | 106.5–121.4 | 5.2 | -0.17 d ± 0.13 | 0.7 | 0.5 | Zd -0.42 (24), Zh -0.17 (13), Zl -0.45 (8), Zv -0.32 (30), Zx -0.13 (73) | [175 genes](qtl_genes/ASI_qtl_genes.tsv) |
| Q8 | chr4 19.8 Mb | 15.2–24.3 | 4.6 | -0.10 d ± 0.07 | 0.9 | 0.5 | Zd +0.15 (11), Zh +0.09 (7), Zl -0.03 (7), Zv -0.32 (32), Zx -0.09 (53) | [163 genes](qtl_genes/ASI_qtl_genes.tsv) |
| Q9 | chr4 52.0 Mb | 46.5–52.4 | 7.5 | -0.24 d ± 0.07 | **3.0** | 0.1 | Zd -0.08 (6), Zh -0.54 (8), Zl -0.67 (6), Zv -0.20 (36), Zx -0.18 (48) | [70 genes](qtl_genes/ASI_qtl_genes.tsv) |
| Q10 | chr4 206.2 Mb | 202.3–225.3 | 5.4 | -0.21 d ± 0.05 | **4.7** | 0.3 | Zd -0.41 (17), Zh -0.38 (5), Zl +0.47 (6), Zv -0.26 (36), Zx -0.15 (77) | [409 genes](qtl_genes/ASI_qtl_genes.tsv) |
| Q11 | chr5 218.3 Mb | 217.4–218.6 | 9.4 | +0.35 d ± 0.06 | **8.8** | **4.1** | Zl +1.08 (11), Zv +0.24 (44), Zx +0.28 (40) | [47 genes](qtl_genes/ASI_qtl_genes.tsv) |
| Q12 | chr7 165.6 Mb | 164.8–166.1 | 13.7 | -0.36 d ± 0.05 | **14.0** | 0.3 | Zd -0.46 (18), Zl -0.53 (14), Zv -0.35 (50), Zx -0.31 (78) | [42 genes](qtl_genes/ASI_qtl_genes.tsv) |

Bonferroni over the 12 QTL: −log10 p > 2.38 (bold). Every gene in each interval, with its GO annotation and the k-mer GWAS −log10 p of its window: [`qtl_genes/ASI_qtl_genes.tsv`](qtl_genes/ASI_qtl_genes.tsv).

</details>

#### 3.6.4 Plant height (PH)

![PH: k-mer ancestry GWAS and joint QTL model](figures/kmerGWAS_joint_qtl_PH.png)

16 regions, 8 hold in the joint model (Bonferroni over 16: −log10 p > 2.51). Effects are in cm per teosinte copy and alike across taxa (no "taxa differ" test passes).

| QTL | Known genes | Effect (joint) | Joint −log10 p | Taxa differ | What it says |
|---|---|---|---|---|---|
| **Q15 chr9 139.8 Mb** | — | +4.73 cm | 10.9 | 0.6 (no) | No significant difference between taxa (Zd +3.87, Zl +1.83, Zv +6.73, Zx +4.18). Same region as NBR Q6. |
| **Q16 chr10 97.0 Mb** | Zm00001eb418700 (ZmCCT10) | +5.64 cm | 10.1 | 0.6 (no) | No significant difference between taxa (Zd +4.49, Zv +8.29, Zx +5.11). The ZmCCT10 region: later flowering, taller plants. |
| **Q14 chr7 139.1 Mb** | — | -4.97 cm | 9.8 | 0.1 (no) | No significant difference between taxa (Zd -5.08, Zh -6.45, Zl -0.91, Zv -4.47, Zx -5.73). Same region as SL Q9. |
| **Q3 chr3 12.8 Mb** | Zm00001eb123060 (ra2) | -3.77 cm | 7.2 | 0.8 (no) | No significant difference between taxa (Zd -5.16, Zh +5.38, Zl +0.26, Zv -3.20, Zx -4.30). Same window as EH Q4. |
| **Q1 chr1 201.8 Mb** | — | -4.17 cm | 6.1 | 0.1 (no) | No significant difference between taxa (Zd -4.82, Zl -1.55, Zv -5.08, Zx -3.29). |
| **Q2 chr2 2.5 Mb** | — | -2.59 cm | 4.9 | 0.3 (no) | No significant difference between taxa (Zd -3.78, Zh -6.14, Zl +0.01, Zv -2.51, Zx -2.04). |
| **Q4 chr4 198.4 Mb** | — | +4.11 cm | 2.9 | 0.1 (no) | No significant difference between taxa (Zd +7.56, Zh +10.21, Zl +14.45, Zv +5.63, Zx +3.49). Same region as EH Q5 and SL Q3. |
| **Q7 chr5 117.8 Mb** | — | -3.60 cm | 2.7 | 0.0 (no) | No significant difference between taxa (Zd +0.14, Zl -3.22, Zv -5.14, Zx -2.69). |

The chr6 peaks (Q8–Q13), Q5 chr4 218.7 Mb and Q6 chr5 64.2 Mb (*td1*, Zm00001eb228140) lose their signal in the joint model. A chr5 2.6 Mb peak (shared with EH and SL) is just below the threshold (−log10 p 4.18 vs 4.19).

<details><summary>All PH QTL with per-taxon effects and all genes in the interval</summary>

| QTL | Position | Interval (Mb) | Scan −log10 p | Effect per copy (joint) | Joint −log10 p | Taxa differ −log10 p | Per-taxon effects (carrier plots) | Genes in interval |
|---|---|---|---|---|---|---|---|---|
| Q1 | chr1 201.8 Mb | 200.4–203.8 | 6.0 | -4.17 cm ± 0.85 | **6.1** | 0.1 | Zd -4.82 (12), Zl -1.55 (8), Zv -5.08 (26), Zx -3.29 (78) | [89 genes](qtl_genes/PH_qtl_genes.tsv) |
| Q2 | chr2 2.5 Mb | 0.0–3.6 | 4.7 | -2.59 cm ± 0.59 | **4.9** | 0.3 | Zd -3.78 (42), Zh -6.14 (8), Zl +0.01 (14), Zv -2.51 (67), Zx -2.04 (109) | [171 genes](qtl_genes/PH_qtl_genes.tsv) |
| Q3 | chr3 12.8 Mb | 12.4–13.7 | 7.0 | -3.77 cm ± 0.70 | **7.2** | 0.8 | Zd -5.16 (29), Zh +5.38 (14), Zl +0.26 (10), Zv -3.20 (32), Zx -4.30 (90) | [41 genes](qtl_genes/PH_qtl_genes.tsv) |
| Q4 | chr4 198.4 Mb | 194.5–201.4 | 7.6 | +4.11 cm ± 1.29 | **2.9** | 0.1 | Zd +7.56 (17), Zh +10.21 (7), Zl +14.45 (6), Zv +5.63 (35), Zx +3.49 (75) | [164 genes](qtl_genes/PH_qtl_genes.tsv) |
| Q5 | chr4 218.7 Mb | 215.2–223.0 | 5.1 | +0.25 cm ± 1.38 | 0.1 | 0.1 | Zd -4.03 (15), Zl +14.45 (6), Zv -2.08 (34), Zx +0.93 (69) | [128 genes](qtl_genes/PH_qtl_genes.tsv) |
| Q6 | chr5 64.2 Mb | 60.6–69.1 | 5.3 | -2.15 cm ± 1.17 | 1.2 | 0.1 | Zd -2.76 (9), Zl -4.58 (12), Zv -0.90 (39), Zx -3.41 (52) | [201 genes](qtl_genes/PH_qtl_genes.tsv) |
| Q7 | chr5 117.8 Mb | 113.8–122.7 | 6.8 | -3.60 cm ± 1.17 | **2.7** | 0.0 | Zd +0.14 (8), Zl -3.22 (8), Zv -5.14 (31), Zx -2.69 (51) | [68 genes](qtl_genes/PH_qtl_genes.tsv) |
| Q8 | chr6 19.9 Mb | 15.8–21.1 | 5.1 | +2.02 cm ± 1.47 | 0.8 | 0.0 | Zd -0.97 (15), Zh +14.39 (7), Zl +2.35 (24), Zv -1.19 (26), Zx +1.56 (51) | [132 genes](qtl_genes/PH_qtl_genes.tsv) |
| Q9 | chr6 34.0 Mb | 33.0–48.4 | 4.6 | +0.03 cm ± 1.62 | 0.0 | 0.6 | Zd +1.11 (15), Zh -10.16 (8), Zl -7.55 (12), Zv +3.54 (29), Zx +1.52 (54) | [218 genes](qtl_genes/PH_qtl_genes.tsv) |
| Q10 | chr6 50.5 Mb | 49.3–51.2 | 5.5 | +3.51 cm ± 2.13 | 1.0 | 0.4 | Zd +5.36 (14), Zh -11.57 (6), Zl +4.68 (10), Zv -5.38 (30), Zx +6.50 (42) | [16 genes](qtl_genes/PH_qtl_genes.tsv) |
| Q11 | chr6 67.2 Mb | 64.1–69.1 | 4.6 | -5.71 cm ± 2.57 | 1.6 | 0.5 | Zd +9.21 (13), Zh -10.79 (6), Zl +10.37 (10), Zv -2.67 (31), Zx -11.41 (51) | [70 genes](qtl_genes/PH_qtl_genes.tsv) |
| Q12 | chr6 77.7 Mb | 76.9–78.0 | 5.7 | +4.06 cm ± 1.97 | 1.4 | 0.9 | Zd -9.49 (13), Zh +8.02 (6), Zl +10.37 (10), Zv +2.52 (22), Zx +7.65 (49) | [11 genes](qtl_genes/PH_qtl_genes.tsv) |
| Q13 | chr6 101.8 Mb | 100.2–106.5 | 5.8 | +2.71 cm ± 1.23 | 1.6 | 2.0 | Zd +4.84 (11), Zh -6.28 (7), Zl -4.14 (7), Zv +7.75 (30), Zx +2.16 (46) | [145 genes](qtl_genes/PH_qtl_genes.tsv) |
| Q14 | chr7 139.1 Mb | 137.1–142.0 | 10.4 | -4.97 cm ± 0.78 | **9.8** | 0.1 | Zd -5.08 (13), Zh -6.45 (8), Zl -0.91 (9), Zv -4.47 (51), Zx -5.73 (69) | [110 genes](qtl_genes/PH_qtl_genes.tsv) |
| Q15 | chr9 139.8 Mb | 132.7–143.5 | 10.3 | +4.73 cm ± 0.70 | **10.9** | 0.6 | Zd +3.87 (26), Zl +1.83 (25), Zv +6.73 (43), Zx +4.18 (88) | [285 genes](qtl_genes/PH_qtl_genes.tsv) |
| Q16 | chr10 97.0 Mb | 91.8–100.0 | 10.8 | +5.64 cm ± 0.87 | **10.1** | 0.6 | Zd +4.49 (23), Zv +8.29 (26), Zx +5.11 (64) | [136 genes](qtl_genes/PH_qtl_genes.tsv) |

Bonferroni over the 16 QTL: −log10 p > 2.51 (bold). Every gene in each interval, with its GO annotation and the k-mer GWAS −log10 p of its window: [`qtl_genes/PH_qtl_genes.tsv`](qtl_genes/PH_qtl_genes.tsv).

</details>

#### 3.6.5 Ear height (EH)

![EH: k-mer ancestry GWAS and joint QTL model](figures/kmerGWAS_joint_qtl_EH.png)

18 regions, 7 hold (Bonferroni over 18: −log10 p > 2.56). Ear height follows flowering: several QTL hold known flowering genes in their intervals.

| QTL | Known genes | Effect (joint) | Joint −log10 p | Taxa differ | What it says |
|---|---|---|---|---|---|
| **Q5 chr4 198.4 Mb** | — | +3.00 cm | 8.1 | 0.1 (no) | No significant difference between taxa (Zd +2.11, Zh -1.06, Zl +6.50, Zv +2.52, Zx +3.82). The strongest EH QTL; shared with PH and SL. |
| **Q10 chr7 184.2 Mb** | Zm00001eb331630 (PRR37-like); Zm00001eb331880 (dlf1) | +2.03 cm | 5.3 | 0.5 (no) | No significant difference between taxa (Zd +2.90, Zl -0.78, Zv +1.39, Zx +2.34). Same region as BL Q6. |
| **Q4 chr3 12.8 Mb** | Zm00001eb122520 (BZR1-like); Zm00001eb122680 (CYP734A6); Zm00001eb123060 (ra2) | -2.06 cm | 5.0 | 2.0 (no) | No significant difference between taxa (Zd -5.48, Zh +0.41, Zl +0.67, Zv -2.95, Zx -1.25). Same window as PH Q3; *ra2*, a BZR1 homolog and CYP734A6 (brassinosteroids) are in the interval. |
| **Q6 chr5 2.5 Mb** | — | +2.17 cm | 5.0 | 0.6 (no) | No significant difference between taxa (Zd +4.04, Zh +3.01, Zl +3.85, Zv +1.05, Zx +1.86). Shared with SL; just below the threshold for PH. |
| **Q9 chr7 124.0 Mb** | — | -2.39 cm | 4.8 | 0.5 (no) | No significant difference between taxa (Zd -0.22, Zh -6.32, Zl -0.22, Zv -1.35, Zx -3.07). |
| **Q1 chr2 25.8 Mb** | Zm00001eb076990 (CKX8); Zm00001eb077080 (GA2ox2) | +1.96 cm | 4.6 | 1.1 (no) | No significant difference between taxa (Zd -0.83, Zh +2.58, Zl +0.83, Zv +1.93, Zx +2.63). *GA2ox2* inactivates gibberellin. |
| **Q11 chr8 126.2 Mb** | Zm00001eb353250 (ZCN8) | +2.36 cm | 3.2 | 0.1 (no) | No significant difference between taxa (Zd -1.30, Zh +4.32, Zl +1.80, Zv +1.03, Zx +2.48). *ZCN8* is the maize florigen. |

Just below the cut: Q13 chr9 116.8 Mb (+1.63 cm, joint −log10 p 2.4; *ZmCCT9*, Zm00001eb391230; Huang et al. 2018) and Q18 chr10 97.0 Mb (+2.40 cm, 2.5; *ZmCCT10*). Q8 chr6 105 Mb is the only QTL where taxa differ (Zv +6.5 cm, Zd +4.2, others ≈ 0 or negative), but its joint test does not pass.

<details><summary>All EH QTL with per-taxon effects and all genes in the interval</summary>

| QTL | Position | Interval (Mb) | Scan −log10 p | Effect per copy (joint) | Joint −log10 p | Taxa differ −log10 p | Per-taxon effects (carrier plots) | Genes in interval |
|---|---|---|---|---|---|---|---|---|
| Q1 | chr2 25.8 Mb | 17.9–30.1 | 6.4 | +1.96 cm ± 0.47 | **4.6** | 1.1 | Zd -0.83 (33), Zh +2.58 (5), Zl +0.83 (21), Zv +1.93 (55), Zx +2.63 (106) | [367 genes](qtl_genes/EH_qtl_genes.tsv) |
| Q2 | chr2 83.2 Mb | 79.5–94.5 | 5.0 | +0.60 cm ± 1.01 | 0.3 | 0.4 | Zd +0.94 (31), Zh +4.97 (7), Zl +3.65 (14), Zv +1.74 (48), Zx -0.21 (101) | [123 genes](qtl_genes/EH_qtl_genes.tsv) |
| Q3 | chr2 135.8 Mb | 130.9–148.8 | 5.1 | +1.00 cm ± 0.99 | 0.5 | 0.3 | Zd -0.48 (27), Zh +5.14 (6), Zl +3.60 (13), Zv +1.22 (50), Zx +1.18 (108) | [198 genes](qtl_genes/EH_qtl_genes.tsv) |
| Q4 | chr3 12.8 Mb | 10.5–19.8 | 5.0 | -2.06 cm ± 0.47 | **5.0** | 2.0 | Zd -5.48 (29), Zh +0.41 (14), Zl +0.67 (10), Zv -2.95 (32), Zx -1.25 (90) | [227 genes](qtl_genes/EH_qtl_genes.tsv) |
| Q5 | chr4 198.4 Mb | 194.4–201.4 | 7.6 | +3.00 cm ± 0.52 | **8.1** | 0.1 | Zd +2.11 (17), Zh -1.06 (7), Zl +6.50 (6), Zv +2.52 (35), Zx +3.82 (75) | [168 genes](qtl_genes/EH_qtl_genes.tsv) |
| Q6 | chr5 2.5 Mb | 1.8–2.9 | 4.8 | +2.17 cm ± 0.50 | **5.0** | 0.6 | Zd +4.04 (33), Zh +3.01 (5), Zl +3.85 (11), Zv +1.05 (47), Zx +1.86 (75) | [90 genes](qtl_genes/EH_qtl_genes.tsv) |
| Q7 | chr6 77.7 Mb | 77.6–99.2 | 4.3 | +1.40 cm ± 0.73 | 1.3 | 0.7 | Zd -1.06 (13), Zh +1.17 (6), Zl +1.75 (10), Zv -0.92 (22), Zx +2.64 (49) | [358 genes](qtl_genes/EH_qtl_genes.tsv) |
| Q8 | chr6 105.0 Mb | 102.9–106.8 | 4.9 | +1.96 cm ± 0.75 | 2.0 | **3.5** | Zd +4.20 (10), Zh -0.54 (7), Zl -3.05 (8), Zv +6.52 (28), Zx -0.14 (44) | [101 genes](qtl_genes/EH_qtl_genes.tsv) |
| Q9 | chr7 124.0 Mb | 122.7–127.4 | 4.7 | -2.39 cm ± 0.56 | **4.8** | 0.5 | Zd -0.22 (11), Zh -6.32 (6), Zl -0.22 (6), Zv -1.35 (50), Zx -3.07 (60) | [97 genes](qtl_genes/EH_qtl_genes.tsv) |
| Q10 | chr7 184.2 Mb | 183.4–185.8 | 5.7 | +2.03 cm ± 0.45 | **5.3** | 0.5 | Zd +2.90 (32), Zl -0.78 (9), Zv +1.39 (53), Zx +2.34 (94) | [82 genes](qtl_genes/EH_qtl_genes.tsv) |
| Q11 | chr8 126.2 Mb | 125.7–133.5 | 6.7 | +2.36 cm ± 0.69 | **3.2** | 0.1 | Zd -1.30 (24), Zh +4.32 (11), Zl +1.80 (9), Zv +1.03 (37), Zx +2.48 (61) | [177 genes](qtl_genes/EH_qtl_genes.tsv) |
| Q12 | chr8 153.6 Mb | 148.8–166.1 | 4.4 | +0.94 cm ± 0.62 | 0.9 | 0.1 | Zd +4.67 (23), Zh +2.84 (7), Zl -0.91 (7), Zv +1.21 (45), Zx +0.83 (61) | [505 genes](qtl_genes/EH_qtl_genes.tsv) |
| Q13 | chr9 116.8 Mb | 114.8–122.2 | 5.3 | +1.63 cm ± 0.58 | 2.4 | 1.1 | Zd +3.11 (22), Zl -2.89 (21), Zv +1.79 (47), Zx +2.12 (89) | [171 genes](qtl_genes/EH_qtl_genes.tsv) |
| Q14 | chr9 142.1 Mb | 131.3–145.8 | 4.7 | +1.21 cm ± 0.56 | 1.5 | 0.0 | Zd +2.51 (25), Zl +3.37 (27), Zv +0.34 (42), Zx +0.94 (82) | [400 genes](qtl_genes/EH_qtl_genes.tsv) |
| Q15 | chr10 6.8 Mb | 5.0–10.2 | 5.3 | +1.73 cm ± 0.64 | 2.2 | 1.8 | Zd +4.24 (23), Zh -4.77 (6), Zl +3.98 (11), Zv +2.07 (34), Zx +1.06 (61) | [142 genes](qtl_genes/EH_qtl_genes.tsv) |
| Q16 | chr10 31.4 Mb | 25.0–33.3 | 4.7 | +1.38 cm ± 1.62 | 0.4 | 0.6 | Zd +0.42 (22), Zl -1.04 (5), Zv +1.01 (23), Zx -1.09 (63) | [86 genes](qtl_genes/EH_qtl_genes.tsv) |
| Q17 | chr10 63.5 Mb | 55.1–67.9 | 4.5 | -1.70 cm ± 1.74 | 0.5 | 0.4 | Zd +1.11 (24), Zl -1.04 (5), Zv -1.76 (25), Zx -0.32 (56) | [126 genes](qtl_genes/EH_qtl_genes.tsv) |
| Q18 | chr10 97.0 Mb | 88.8–100.5 | 6.9 | +2.40 cm ± 0.81 | 2.5 | 0.8 | Zd -1.05 (23), Zv +3.96 (26), Zx +3.13 (65) | [197 genes](qtl_genes/EH_qtl_genes.tsv) |

Bonferroni over the 18 QTL: −log10 p > 2.56 (bold). Every gene in each interval, with its GO annotation and the k-mer GWAS −log10 p of its window: [`qtl_genes/EH_qtl_genes.tsv`](qtl_genes/EH_qtl_genes.tsv).

</details>

#### 3.6.6 Ear number (EN; number of nodes with ears)

![EN: k-mer ancestry GWAS and joint QTL model](figures/kmerGWAS_joint_qtl_EN.png)

19 regions plus one QTL added by the conditional scans (Q8 chr4 120.2 Mb, open triangle); 11 of the 20 hold (Bonferroni over 20: −log10 p > 2.60). EN counts the stem nodes that bear an ear. **Every QTL that holds is positive: teosinte alleles add ear-bearing nodes**, as expected from teosinte, which bears ears at many nodes.

| QTL | Known genes | Effect (joint) | Joint −log10 p | Taxa differ | What it says |
|---|---|---|---|---|---|
| **Q3 chr2 22.4 Mb** | — | +0.17 | 16.6 | 1.1 (no) | No significant difference between taxa (Zd +0.06, Zh +0.20, Zl +0.14, Zv +0.18, Zx +0.20). |
| **Q1 chr1 22.9 Mb** | Zm00001eb007670 (ZmD14) | +0.19 | 14.9 | 2.2 (no) | No significant difference between taxa (Zd +0.19, Zh +0.15, Zl -0.12, Zv +0.26, Zx +0.20). *ZmD14* is the strigolactone receptor (shoot branching). |
| **Q4 chr2 202.3 Mb** | — | +0.14 | 12.3 | 2.7 (yes) | Differs between taxa (Zd -0.08, Zl +0.05, Zv +0.17, Zx +0.16). Not in *diploperennis*. |
| **Q13 chr7 143.6 Mb** | — | +0.16 | 11.1 | 0.1 (no) | No significant difference between taxa (Zd +0.09, Zh +0.02, Zl +0.06, Zv +0.19, Zx +0.14). |
| **Q9 chr4 166.9 Mb** | — | +0.14 | 6.4 | 0.3 (no) | No significant difference between taxa (Zd +0.19, Zh +0.07, Zl -0.13, Zv +0.08, Zx +0.19). Same region as Prolif Q13. |
| **Q8 chr4 120.2 Mb** (added) | — | +0.16 | 6.1 | 0.6 (no) | No significant difference between taxa (Zh +0.05, Zv +0.20, Zx +0.12). Second chr4 QTL besides Q9, found by the conditional scan; carriers only from Zh, Zv and Zx; wide interval (89.8–128.3 Mb) because recombination is low there. |
| **Q14 chr8 33.8 Mb** | — | +0.12 | 5.4 | 1.1 (no) | No significant difference between taxa (Zd +0.03, Zh -0.01, Zl +0.05, Zv +0.12, Zx +0.17). |
| **Q20 chr10 143.9 Mb** | — | +0.09 | 4.9 | 0.4 (no) | No significant difference between taxa (Zd +0.05, Zh +0.16, Zl +0.10, Zv +0.08, Zx +0.11). |
| **Q12 chr6 173.2 Mb** | — | +0.09 | 3.8 | 0.3 (no) | No significant difference between taxa (Zd +0.04, Zl +0.07, Zv +0.06, Zx +0.13). |
| **Q10 chr6 97.0 Mb** | — | +0.10 | 3.7 | 0.2 (no) | No significant difference between taxa (Zd +0.10, Zh +0.13, Zl +0.16, Zv +0.07, Zx +0.12). |
| **Q2 chr1 81.3 Mb** | Zm00001eb020350 (clock-coupling ZGT) | +0.09 | 3.3 | 2.0 (no) | No significant difference between taxa (Zd +0.01, Zh -0.03, Zl -0.04, Zv +0.14, Zx +0.11). |

Q5 chr3 12.4 Mb (+0.07) holds *ra2* (Zm00001eb123060) and Q15 chr8 133.6 Mb holds *Rap2.7* (Zm00001eb355240) and *lg4a* (Zm00001eb354880); neither holds in the joint model.

<details><summary>All EN QTL with per-taxon effects and all genes in the interval</summary>

| QTL | Position | Interval (Mb) | Scan −log10 p | Effect per copy (joint) | Joint −log10 p | Taxa differ −log10 p | Per-taxon effects (carrier plots) | Genes in interval |
|---|---|---|---|---|---|---|---|---|
| Q1 | chr1 22.9 Mb | 22.6–23.6 | 16.6 | +0.19 ± 0.02 | **14.9** | 2.2 | Zd +0.19 (15), Zh +0.15 (10), Zl -0.12 (13), Zv +0.26 (35), Zx +0.20 (59) | [33 genes](qtl_genes/EN_qtl_genes.tsv) |
| Q2 | chr1 81.3 Mb | 58.0–85.9 | 5.1 | +0.09 ± 0.03 | **3.3** | 2.0 | Zd +0.01 (13), Zh -0.03 (6), Zl -0.04 (9), Zv +0.14 (22), Zx +0.11 (63) | [585 genes](qtl_genes/EN_qtl_genes.tsv) |
| Q3 | chr2 22.4 Mb | 22.1–23.0 | 18.6 | +0.17 ± 0.02 | **16.6** | 1.1 | Zd +0.06 (27), Zh +0.20 (5), Zl +0.14 (20), Zv +0.18 (49), Zx +0.20 (102) | [28 genes](qtl_genes/EN_qtl_genes.tsv) |
| Q4 | chr2 202.3 Mb | 201.6–205.0 | 14.4 | +0.14 ± 0.02 | **12.3** | **2.7** | Zd -0.08 (26), Zl +0.05 (11), Zv +0.17 (53), Zx +0.16 (106) | [123 genes](qtl_genes/EN_qtl_genes.tsv) |
| Q5 | chr3 12.4 Mb | 10.9–14.0 | 6.4 | +0.07 ± 0.03 | 2.5 | 0.1 | Zd +0.08 (29), Zh +0.09 (13), Zl +0.11 (10), Zv +0.07 (32), Zx +0.04 (93) | [74 genes](qtl_genes/EN_qtl_genes.tsv) |
| Q6 | chr3 37.0 Mb | 24.4–43.4 | 5.3 | +0.04 ± 0.03 | 0.8 | 0.8 | Zd -0.02 (24), Zh -0.04 (11), Zl -0.15 (8), Zv +0.04 (27), Zx +0.09 (74) | [314 genes](qtl_genes/EN_qtl_genes.tsv) |
| Q7 | chr3 175.2 Mb | 171.4–176.8 | 4.9 | +0.07 ± 0.02 | 2.5 | 2.1 | Zd -0.07 (24), Zh +0.00 (10), Zl +0.07 (14), Zv +0.02 (39), Zx +0.14 (76) | [131 genes](qtl_genes/EN_qtl_genes.tsv) |
| Q8 (added) | chr4 120.2 Mb | 89.8–128.3 | 14.4 | +0.16 ± 0.03 | **6.1** | 0.6 | Zh +0.05 (6), Zv +0.20 (33), Zx +0.12 (51) | [272 genes](qtl_genes/EN_qtl_genes.tsv) |
| Q9 | chr4 166.9 Mb | 166.2–168.9 | 14.8 | +0.14 ± 0.03 | **6.4** | 0.3 | Zd +0.19 (7), Zh +0.07 (5), Zl -0.13 (5), Zv +0.08 (38), Zx +0.19 (68) | [79 genes](qtl_genes/EN_qtl_genes.tsv) |
| Q10 | chr6 97.0 Mb | 95.0–99.2 | 5.5 | +0.10 ± 0.03 | **3.7** | 0.2 | Zd +0.10 (10), Zh +0.13 (7), Zl +0.16 (7), Zv +0.07 (32), Zx +0.12 (45) | [100 genes](qtl_genes/EN_qtl_genes.tsv) |
| Q11 | chr6 160.8 Mb | 153.6–163.5 | 4.6 | +0.07 ± 0.03 | 2.1 | 0.2 | Zd +0.03 (17), Zh +0.14 (5), Zl -0.04 (6), Zv +0.14 (34), Zx +0.04 (57) | [314 genes](qtl_genes/EN_qtl_genes.tsv) |
| Q12 | chr6 173.2 Mb | 171.6–175.2 | 5.9 | +0.09 ± 0.03 | **3.8** | 0.3 | Zd +0.04 (16), Zl +0.07 (8), Zv +0.06 (43), Zx +0.13 (56) | [186 genes](qtl_genes/EN_qtl_genes.tsv) |
| Q13 | chr7 143.6 Mb | 142.8–144.4 | 10.8 | +0.16 ± 0.02 | **11.1** | 0.1 | Zd +0.09 (16), Zh +0.02 (8), Zl +0.06 (13), Zv +0.19 (49), Zx +0.14 (77) | [43 genes](qtl_genes/EN_qtl_genes.tsv) |
| Q14 | chr8 33.8 Mb | 29.5–37.7 | 8.4 | +0.12 ± 0.03 | **5.4** | 1.1 | Zd +0.03 (28), Zh -0.01 (11), Zl +0.05 (6), Zv +0.12 (54), Zx +0.17 (61) | [103 genes](qtl_genes/EN_qtl_genes.tsv) |
| Q15 | chr8 133.6 Mb | 127.8–139.6 | 4.6 | +0.06 ± 0.03 | 1.6 | 0.3 | Zd +0.00 (21), Zh +0.09 (11), Zl +0.07 (8), Zv +0.06 (34), Zx +0.04 (64) | [287 genes](qtl_genes/EN_qtl_genes.tsv) |
| Q16 | chr9 32.4 Mb | 31.5–32.5 | 7.5 | +0.11 ± 0.04 | 2.2 | 0.5 | Zd +0.13 (24), Zh +0.18 (8), Zl -0.03 (22), Zv +0.11 (48), Zx +0.16 (94) | [10 genes](qtl_genes/EN_qtl_genes.tsv) |
| Q17 | chr9 70.0 Mb | 68.2–70.4 | 5.6 | +0.00 ± 0.07 | 0.0 | 0.5 | Zd +0.11 (23), Zh -0.27 (7), Zl +0.05 (21), Zv -0.04 (47), Zx +0.06 (92) | [24 genes](qtl_genes/EN_qtl_genes.tsv) |
| Q18 | chr9 80.0 Mb | 75.1–86.1 | 5.4 | +0.01 ± 0.07 | 0.0 | 1.1 | Zd +0.04 (26), Zh -0.35 (7), Zl +0.23 (21), Zv -0.09 (50), Zx +0.13 (84) | [99 genes](qtl_genes/EN_qtl_genes.tsv) |
| Q19 | chr9 92.5 Mb | 90.3–101.5 | 4.9 | -0.01 ± 0.05 | 0.1 | 0.8 | Zd -0.17 (21), Zh +0.04 (6), Zl +0.14 (19), Zv +0.07 (49), Zx -0.14 (98) | [162 genes](qtl_genes/EN_qtl_genes.tsv) |
| Q20 | chr10 143.9 Mb | 143.5–146.2 | 5.8 | +0.09 ± 0.02 | **4.9** | 0.4 | Zd +0.05 (42), Zh +0.16 (8), Zl +0.10 (12), Zv +0.08 (43), Zx +0.11 (79) | [115 genes](qtl_genes/EN_qtl_genes.tsv) |

Bonferroni over the 20 QTL: −log10 p > 2.60 (bold). Every gene in each interval, with its GO annotation and the k-mer GWAS −log10 p of its window: [`qtl_genes/EN_qtl_genes.tsv`](qtl_genes/EN_qtl_genes.tsv).

</details>

#### 3.6.7 Prolificacy (Prolif; most ears on one node)

![Prolif: k-mer ancestry GWAS and joint QTL model](figures/kmerGWAS_joint_qtl_Prolif.png)

12 regions plus three QTL added by the conditional scans (Q5 chr2 9.3 Mb, Q6 chr2 34.4 Mb and Q12 chr4 120.5 Mb, open triangles); 10 of the 15 hold (Bonferroni over 15: −log10 p > 2.48). Prolificacy is the largest number of ears on a single node (teosinte clusters several ears per node). **Largely the same loci as EN, all positive** (teosinte adds ears per node): the same loci raise both the number of ear-bearing nodes and the number of ears per node.

| QTL | Known genes | Effect (joint) | Joint −log10 p | Taxa differ | What it says |
|---|---|---|---|---|---|
| **Q1 chr1 23.2 Mb** | — | +0.20 | 18.0 | 0.9 (no) | No significant difference between taxa (Zd +0.17, Zh +0.08, Zl -0.03, Zv +0.25, Zx +0.22). Same region as EN Q1. |
| **Q8 chr3 11.4 Mb** | Zm00001eb122680 (CYP734A6) | +0.16 | 13.3 | 0.5 (no) | No significant difference between taxa (Zd +0.19, Zh +0.12, Zl +0.08, Zv +0.21, Zx +0.13). |
| **Q7 chr2 191.2 Mb** | — | +0.12 | 10.3 | 1.3 (no) | No significant difference between taxa (Zd -0.02, Zh +0.20, Zl +0.02, Zv +0.09, Zx +0.16). |
| **Q13 chr4 171.3 Mb** | — | +0.14 | 7.8 | 1.4 (no) | No significant difference between taxa (Zd +0.21, Zh -0.04, Zl +0.10, Zv +0.13, Zx +0.17). Same region as EN Q9. |
| **Q15 chr9 32.1 Mb** | — | +0.09 | 5.8 | 0.1 (no) | No significant difference between taxa (Zd +0.08, Zh -0.02, Zl +0.13, Zv +0.07, Zx +0.10). |
| **Q5 chr2 9.3 Mb** (added) | Zm00001eb076990 (CKX8); Zm00001eb077080 (GA2ox2) | +0.08 | 5.2 | 0.8 (no) | No significant difference between taxa (Zd +0.13, Zh +0.14, Zl +0.11, Zv +0.02, Zx +0.07). Added by the conditional scan: a second QTL on chr2 short arm. |
| **Q12 chr4 120.5 Mb** (added) | — | +0.13 | 4.9 | 0.1 (no) | No significant difference between taxa (Zh +0.09, Zv +0.13, Zx +0.13). The same second chr4 QTL as EN Q8, found by the conditional scan. |
| **Q6 chr2 34.4 Mb** (added) | Zm00001eb076990 (CKX8); Zm00001eb077080 (GA2ox2) | +0.08 | 4.7 | 0.7 (no) | No significant difference between taxa (Zd +0.07, Zl +0.01, Zv +0.13, Zx +0.09). Added by the conditional scan: a third chr2 QTL. |
| **Q14 chr7 142.2 Mb** | — | +0.08 | 4.4 | 1.3 (no) | No significant difference between taxa (Zd +0.00, Zh +0.33, Zl +0.01, Zv +0.06, Zx +0.09). |
| **Q4 chr1 213.9 Mb** | Zm00001eb040280 (CKX10) | +0.10 | 3.1 | 0.7 (no) | No significant difference between taxa (Zd +0.18, Zh +0.11, Zl -0.03, Zv +0.10, Zx +0.10). *CKX10* is a cytokinin oxidase. |

<details><summary>All Prolif QTL with per-taxon effects and all genes in the interval</summary>

| QTL | Position | Interval (Mb) | Scan −log10 p | Effect per copy (joint) | Joint −log10 p | Taxa differ −log10 p | Per-taxon effects (carrier plots) | Genes in interval |
|---|---|---|---|---|---|---|---|---|
| Q1 | chr1 23.2 Mb | 22.9–23.3 | 20.3 | +0.20 ± 0.02 | **18.0** | 0.9 | Zd +0.17 (15), Zh +0.08 (10), Zl -0.03 (13), Zv +0.25 (33), Zx +0.22 (58) | [12 genes](qtl_genes/Prolif_qtl_genes.tsv) |
| Q2 | chr1 61.9 Mb | 56.1–64.5 | 4.6 | +0.04 ± 0.02 | 1.0 | 0.8 | Zd +0.00 (19), Zh +0.11 (10), Zl +0.01 (14), Zv +0.01 (27), Zx +0.02 (65) | [206 genes](qtl_genes/Prolif_qtl_genes.tsv) |
| Q3 | chr1 205.3 Mb | 203.8–206.6 | 4.7 | +0.02 ± 0.03 | 0.4 | 1.1 | Zd -0.18 (8), Zh -0.30 (5), Zl +0.05 (8), Zv +0.06 (33), Zx +0.05 (79) | [69 genes](qtl_genes/Prolif_qtl_genes.tsv) |
| Q4 | chr1 213.9 Mb | 211.3–218.3 | 5.8 | +0.10 ± 0.03 | **3.1** | 0.7 | Zd +0.18 (7), Zh +0.11 (7), Zl -0.03 (7), Zv +0.10 (33), Zx +0.10 (74) | [142 genes](qtl_genes/Prolif_qtl_genes.tsv) |
| Q5 (added) | chr2 9.3 Mb | 6.2–35.4 | 6.7 | +0.08 ± 0.02 | **5.2** | 0.8 | Zd +0.13 (35), Zh +0.14 (10), Zl +0.11 (17), Zv +0.02 (55), Zx +0.07 (98) | [1004 genes](qtl_genes/Prolif_qtl_genes.tsv) |
| Q6 (added) | chr2 34.4 Mb | 22.1–38.2 | 9.6 | +0.08 ± 0.02 | **4.7** | 0.7 | Zd +0.07 (33), Zl +0.01 (18), Zv +0.13 (51), Zx +0.09 (93) | [438 genes](qtl_genes/Prolif_qtl_genes.tsv) |
| Q7 | chr2 191.2 Mb | 185.7–192.0 | 12.2 | +0.12 ± 0.02 | **10.3** | 1.3 | Zd -0.02 (17), Zh +0.20 (6), Zl +0.02 (16), Zv +0.09 (51), Zx +0.16 (112) | [149 genes](qtl_genes/Prolif_qtl_genes.tsv) |
| Q8 | chr3 11.4 Mb | 10.9–11.7 | 19.0 | +0.16 ± 0.02 | **13.3** | 0.5 | Zd +0.19 (29), Zh +0.12 (13), Zl +0.08 (10), Zv +0.21 (36), Zx +0.13 (91) | [13 genes](qtl_genes/Prolif_qtl_genes.tsv) |
| Q9 | chr3 118.8 Mb | 113.5–121.7 | 6.8 | +0.03 ± 0.03 | 0.4 | 0.9 | Zd +0.13 (26), Zh -0.09 (13), Zl -0.16 (8), Zv -0.10 (28), Zx +0.07 (75) | [96 genes](qtl_genes/Prolif_qtl_genes.tsv) |
| Q10 | chr3 148.7 Mb | 146.7–154.3 | 5.0 | -0.01 ± 0.03 | 0.2 | 0.1 | Zd -0.12 (27), Zh +0.04 (13), Zl +0.12 (11), Zv +0.13 (29), Zx -0.06 (91) | [174 genes](qtl_genes/Prolif_qtl_genes.tsv) |
| Q11 | chr3 172.4 Mb | 166.7–177.1 | 5.3 | +0.07 ± 0.02 | 2.1 | 1.4 | Zd +0.04 (24), Zh +0.06 (12), Zl -0.03 (11), Zv -0.04 (36), Zx +0.14 (83) | [218 genes](qtl_genes/Prolif_qtl_genes.tsv) |
| Q12 (added) | chr4 120.5 Mb | 75.0–140.2 | 13.1 | +0.13 ± 0.03 | **4.9** | 0.1 | Zh +0.09 (6), Zv +0.13 (32), Zx +0.13 (51) | [590 genes](qtl_genes/Prolif_qtl_genes.tsv) |
| Q13 | chr4 171.3 Mb | 167.6–171.9 | 16.5 | +0.14 ± 0.03 | **7.8** | 1.4 | Zd +0.21 (11), Zh -0.04 (5), Zl +0.10 (6), Zv +0.13 (31), Zx +0.17 (66) | [111 genes](qtl_genes/Prolif_qtl_genes.tsv) |
| Q14 | chr7 142.2 Mb | 141.6–147.9 | 5.0 | +0.08 ± 0.02 | **4.4** | 1.3 | Zd +0.00 (15), Zh +0.33 (8), Zl +0.01 (10), Zv +0.06 (48), Zx +0.09 (73) | [172 genes](qtl_genes/Prolif_qtl_genes.tsv) |
| Q15 | chr9 32.1 Mb | 25.4–32.7 | 5.4 | +0.09 ± 0.02 | **5.8** | 0.1 | Zd +0.08 (24), Zh -0.02 (8), Zl +0.13 (22), Zv +0.07 (48), Zx +0.10 (95) | [156 genes](qtl_genes/Prolif_qtl_genes.tsv) |

Bonferroni over the 15 QTL: −log10 p > 2.48 (bold). Every gene in each interval, with its GO annotation and the k-mer GWAS −log10 p of its window: [`qtl_genes/Prolif_qtl_genes.tsv`](qtl_genes/Prolif_qtl_genes.tsv).

</details>

#### 3.6.8 Nodes with brace roots (NBR)

![NBR: k-mer ancestry GWAS and joint QTL model](figures/kmerGWAS_joint_qtl_NBR.png)

Number of stem nodes that carry brace roots. 8 regions, 4 hold (Bonferroni over 8: −log10 p > 2.20).

| QTL | Known genes | Effect (joint) | Joint −log10 p | Taxa differ | What it says |
|---|---|---|---|---|---|
| **Q8 chr10 98.2 Mb** | — | +0.17 | 5.7 | 3.1 (yes) | Differs between taxa (Zd +0.05, Zv +0.26, Zx +0.09). Next to the ZmCCT10 region: later-flowering plants with more nodes; largest in *parviglumis*. |
| **Q6 chr9 139.9 Mb** | Zm00001eb397240 (PRR73) | +0.09 | 4.6 | 0.1 (no) | No significant difference between taxa (Zd +0.01, Zl +0.14, Zv +0.07, Zx +0.12). Same region as PH Q15. |
| **Q5 chr8 126.2 Mb** | Zm00001eb353250 (ZCN8) | +0.12 | 4.6 | 1.1 (no) | No significant difference between taxa (Zd +0.08, Zh +0.33, Zl +0.14, Zv +0.09, Zx +0.10). *ZCN8* is the maize florigen. |
| **Q4 chr1 262.4 Mb** | — | -0.09 | 3.1 | 1.7 (no) | No significant difference between taxa (Zd -0.21, Zh +0.01, Zl -0.17, Zv -0.16, Zx -0.04). |

The other chr1 peaks (negative effects) share their signal and do not hold individually. Two QTL have no overall effect but differ between taxa: **Q7 chr10 25.4 Mb** (joint −log10 p 0.2, taxa differ 6.1: only *diploperennis* carriers add brace-root nodes, +0.40) and **Q2 chr1 23.4 Mb** (taxa differ 2.6; *ZmD14*, Zm00001eb007670, is in the interval).

<details><summary>All NBR QTL with per-taxon effects and all genes in the interval</summary>

| QTL | Position | Interval (Mb) | Scan −log10 p | Effect per copy (joint) | Joint −log10 p | Taxa differ −log10 p | Per-taxon effects (carrier plots) | Genes in interval |
|---|---|---|---|---|---|---|---|---|
| Q1 | chr1 11.1 Mb | 9.3–13.0 | 4.7 | -0.08 ± 0.03 | 2.0 | 0.5 | Zd -0.06 (17), Zh +0.02 (11), Zl -0.03 (6), Zv -0.04 (37), Zx -0.08 (64) | [125 genes](qtl_genes/NBR_qtl_genes.tsv) |
| Q2 | chr1 23.4 Mb | 21.9–25.6 | 4.4 | -0.07 ± 0.03 | 1.6 | **2.6** | Zd +0.10 (14), Zh -0.37 (10), Zl -0.08 (13), Zv -0.14 (33), Zx -0.04 (58) | [97 genes](qtl_genes/NBR_qtl_genes.tsv) |
| Q3 | chr1 205.3 Mb | 199.0–211.1 | 4.4 | -0.07 ± 0.03 | 2.2 | 1.4 | Zd -0.28 (8), Zh -0.28 (5), Zl +0.00 (8), Zv -0.00 (33), Zx -0.07 (78) | [323 genes](qtl_genes/NBR_qtl_genes.tsv) |
| Q4 | chr1 262.4 Mb | 257.2–264.5 | 5.0 | -0.09 ± 0.03 | **3.1** | 1.7 | Zd -0.21 (13), Zh +0.01 (6), Zl -0.17 (8), Zv -0.16 (34), Zx -0.04 (77) | [184 genes](qtl_genes/NBR_qtl_genes.tsv) |
| Q5 | chr8 126.2 Mb | 125.0–127.6 | 4.4 | +0.12 ± 0.03 | **4.6** | 1.1 | Zd +0.08 (23), Zh +0.33 (10), Zl +0.14 (9), Zv +0.09 (36), Zx +0.10 (61) | [66 genes](qtl_genes/NBR_qtl_genes.tsv) |
| Q6 | chr9 139.9 Mb | 131.3–145.4 | 4.3 | +0.09 ± 0.02 | **4.6** | 0.1 | Zd +0.01 (26), Zl +0.14 (25), Zv +0.07 (43), Zx +0.12 (88) | [386 genes](qtl_genes/NBR_qtl_genes.tsv) |
| Q7 | chr10 25.4 Mb | 25.0–26.9 | 5.2 | +0.02 ± 0.04 | 0.2 | **6.1** | Zd +0.40 (22), Zv +0.00 (24), Zx -0.02 (57) | [34 genes](qtl_genes/NBR_qtl_genes.tsv) |
| Q8 | chr10 98.2 Mb | 96.9–99.6 | 10.3 | +0.17 ± 0.04 | **5.7** | **3.1** | Zd +0.05 (23), Zv +0.26 (25), Zx +0.09 (64) | [50 genes](qtl_genes/NBR_qtl_genes.tsv) |

Bonferroni over the 8 QTL: −log10 p > 2.20 (bold). Every gene in each interval, with its GO annotation and the k-mer GWAS −log10 p of its window: [`qtl_genes/NBR_qtl_genes.tsv`](qtl_genes/NBR_qtl_genes.tsv).

</details>

#### 3.6.9 Sheath length (SL; 2023 only)

![SL: k-mer ancestry GWAS and joint QTL model](figures/kmerGWAS_joint_qtl_SL.png)

13 regions, 12 hold (Bonferroni over 13: −log10 p > 2.41). Sheath length shares several regions with plant height, as expected: the sheath wraps the internode, and longer internodes make both taller plants and longer sheaths.

| QTL | Known genes | Effect (joint) | Joint −log10 p | Taxa differ | What it says |
|---|---|---|---|---|---|
| **Q9 chr7 141.7 Mb** | — | -0.32 | 9.6 | 0.9 (no) | No significant difference between taxa (Zd -0.26, Zh +0.07, Zl -0.58, Zv -0.30, Zx -0.37). Same region as PH Q14. |
| **Q2 chr2 6.2 Mb** | — | -0.24 | 8.1 | 0.8 (no) | No significant difference between taxa (Zd -0.40, Zh -0.17, Zl -0.55, Zv -0.16, Zx -0.17). |
| **Q10 chr9 157.2 Mb** | — | -0.25 | 8.0 | 2.2 (no) | No significant difference between taxa (Zd -0.66, Zh -0.13, Zl -0.32, Zv -0.29, Zx -0.13). |
| **Q4 chr5 2.5 Mb** | — | +0.26 | 6.8 | 1.1 (no) | No significant difference between taxa (Zd +0.19, Zh +0.88, Zl -0.06, Zv +0.20, Zx +0.27). Same region as EH Q6. |
| **Q3 chr4 200.6 Mb** | — | +0.28 | 6.6 | 0.3 (no) | No significant difference between taxa (Zd +0.12, Zl +0.14, Zv +0.19, Zx +0.37). Same region as PH Q4 and EH Q5. |
| **Q13 chr10 151.3 Mb** | — | -0.24 | 6.3 | 0.9 (no) | No significant difference between taxa (Zd -0.05, Zl -0.23, Zv -0.14, Zx -0.36). |
| **Q1 chr1 69.7 Mb** | — | +0.28 | 5.5 | 0.1 (no) | No significant difference between taxa (Zd +0.26, Zh +0.53, Zl -0.22, Zv +0.36, Zx +0.24). |
| **Q7 chr6 116.0 Mb** | — | -0.25 | 4.8 | 0.7 (no) | No significant difference between taxa (Zd -0.31, Zh -0.53, Zl -0.09, Zv -0.19, Zx -0.12). |
| **Q8 chr6 167.2 Mb** | — | -0.22 | 4.0 | 0.2 (no) | No significant difference between taxa (Zd -0.26, Zh -0.47, Zl -0.32, Zv -0.17, Zx -0.24). |
| **Q11 chr10 25.6 Mb** | — | +0.26 | 3.8 | 1.7 (no) | No significant difference between taxa (Zd +0.57, Zv +0.18, Zx +0.15). |
| **Q6 chr5 92.0 Mb** | — | -0.26 | 3.3 | 0.4 (no) | No significant difference between taxa (Zd -0.38, Zl -0.03, Zv -0.04, Zx -0.34). |
| **Q12 chr10 129.7 Mb** | — | +0.20 | 3.2 | 0.2 (no) | No significant difference between taxa (Zd +0.23, Zh +0.43, Zl -0.10, Zv +0.18, Zx +0.29). |

<details><summary>All SL QTL with per-taxon effects and all genes in the interval</summary>

| QTL | Position | Interval (Mb) | Scan −log10 p | Effect per copy (joint) | Joint −log10 p | Taxa differ −log10 p | Per-taxon effects (carrier plots) | Genes in interval |
|---|---|---|---|---|---|---|---|---|
| Q1 | chr1 69.7 Mb | 67.6–71.1 | 5.2 | +0.28 ± 0.06 | **5.5** | 0.1 | Zd +0.26 (18), Zh +0.53 (10), Zl -0.22 (7), Zv +0.36 (19), Zx +0.24 (50) | [74 genes](qtl_genes/SL_qtl_genes.tsv) |
| Q2 | chr2 6.2 Mb | 5.8–7.4 | 7.4 | -0.24 ± 0.04 | **8.1** | 0.8 | Zd -0.40 (41), Zh -0.17 (11), Zl -0.55 (13), Zv -0.16 (56), Zx -0.17 (93) | [92 genes](qtl_genes/SL_qtl_genes.tsv) |
| Q3 | chr4 200.6 Mb | 197.7–207.2 | 6.5 | +0.28 ± 0.05 | **6.6** | 0.3 | Zd +0.12 (16), Zl +0.14 (5), Zv +0.19 (36), Zx +0.37 (66) | [234 genes](qtl_genes/SL_qtl_genes.tsv) |
| Q4 | chr5 2.5 Mb | 1.8–2.7 | 6.1 | +0.26 ± 0.05 | **6.8** | 1.1 | Zd +0.19 (31), Zh +0.88 (5), Zl -0.06 (11), Zv +0.20 (45), Zx +0.27 (70) | [76 genes](qtl_genes/SL_qtl_genes.tsv) |
| Q5 | chr5 48.6 Mb | 47.9–50.4 | 5.7 | -0.18 ± 0.08 | 1.7 | 1.7 | Zd -0.34 (10), Zl -0.20 (13), Zv -0.48 (33), Zx +0.03 (40) | [40 genes](qtl_genes/SL_qtl_genes.tsv) |
| Q6 | chr5 92.0 Mb | 91.8–92.5 | 8.0 | -0.26 ± 0.07 | **3.3** | 0.4 | Zd -0.38 (8), Zl -0.03 (9), Zv -0.04 (32), Zx -0.34 (46) | [6 genes](qtl_genes/SL_qtl_genes.tsv) |
| Q7 | chr6 116.0 Mb | 113.6–118.1 | 5.0 | -0.25 ± 0.06 | **4.8** | 0.7 | Zd -0.31 (12), Zh -0.53 (7), Zl -0.09 (9), Zv -0.19 (23), Zx -0.12 (39) | [117 genes](qtl_genes/SL_qtl_genes.tsv) |
| Q8 | chr6 167.2 Mb | 165.6–171.7 | 4.6 | -0.22 ± 0.06 | **4.0** | 0.2 | Zd -0.26 (17), Zh -0.47 (8), Zl -0.32 (8), Zv -0.17 (35), Zx -0.24 (45) | [286 genes](qtl_genes/SL_qtl_genes.tsv) |
| Q9 | chr7 141.7 Mb | 139.6–141.9 | 9.2 | -0.32 ± 0.05 | **9.6** | 0.9 | Zd -0.26 (12), Zh +0.07 (8), Zl -0.58 (9), Zv -0.30 (44), Zx -0.37 (70) | [38 genes](qtl_genes/SL_qtl_genes.tsv) |
| Q10 | chr9 157.2 Mb | 156.7–157.6 | 8.0 | -0.25 ± 0.04 | **8.0** | 2.2 | Zd -0.66 (21), Zh -0.13 (10), Zl -0.32 (14), Zv -0.29 (49), Zx -0.13 (93) | [45 genes](qtl_genes/SL_qtl_genes.tsv) |
| Q11 | chr10 25.6 Mb | 25.1–28.9 | 7.2 | +0.26 ± 0.07 | **3.8** | 1.7 | Zd +0.57 (22), Zv +0.18 (23), Zx +0.15 (53) | [52 genes](qtl_genes/SL_qtl_genes.tsv) |
| Q12 | chr10 129.7 Mb | 125.9–131.9 | 7.2 | +0.20 ± 0.06 | **3.2** | 0.2 | Zd +0.23 (32), Zh +0.43 (8), Zl -0.10 (5), Zv +0.18 (35), Zx +0.29 (61) | [182 genes](qtl_genes/SL_qtl_genes.tsv) |
| Q13 | chr10 151.3 Mb | 149.5–151.8 | 4.5 | -0.24 ± 0.05 | **6.3** | 0.9 | Zd -0.05 (31), Zl -0.23 (12), Zv -0.14 (43), Zx -0.36 (78) | [133 genes](qtl_genes/SL_qtl_genes.tsv) |

Bonferroni over the 13 QTL: −log10 p > 2.41 (bold). Every gene in each interval, with its GO annotation and the k-mer GWAS −log10 p of its window: [`qtl_genes/SL_qtl_genes.tsv`](qtl_genes/SL_qtl_genes.tsv).

</details>

#### 3.6.10 Blade length (BL; 2023 only)

![BL: k-mer ancestry GWAS and joint QTL model](figures/kmerGWAS_joint_qtl_BL.png)

5 regions plus one QTL added by the conditional scans (Q1 chr2 9.4 Mb, open triangle); 5 of the 6 hold (Bonferroni over 6: −log10 p > 2.08). **chr2 has two QTL** (Q1 and Q2), which were one region in the scan.

| QTL | Known genes | Effect (joint) | Joint −log10 p | Taxa differ | What it says |
|---|---|---|---|---|---|
| **Q2 chr2 39.8 Mb** | Zm00001eb079900 (FTIP1) | -1.58 | 8.8 | 1.8 (no) | No significant difference between taxa (Zd -0.83, Zl +1.52, Zv -1.67, Zx -1.92). Not in *luxurians*; *FT-INTERACTING PROTEIN 1* is in the interval. |
| **Q5 chr4 197.1 Mb** | — | -1.74 | 7.8 | 0.0 (no) | No significant difference between taxa (Zd -0.57, Zh -2.69, Zl -3.47, Zv -1.91, Zx -1.75). The PH/EH/SL chr4 region. |
| **Q1 chr2 9.4 Mb** (added) | — | -1.37 | 7.3 | 0.0 (no) | No significant difference between taxa (Zd -1.68, Zh -1.58, Zl -2.07, Zv -1.61, Zx -1.01). Second chr2 QTL, found by the conditional scan. |
| **Q6 chr7 185.2 Mb** | Zm00001eb331630 (PRR37-like); Zm00001eb331880 (dlf1) | -1.16 | 5.0 | 1.1 (no) | No significant difference between taxa (Zd +0.19, Zl +0.60, Zv -1.25, Zx -1.71). Same region as EH Q10. |
| **Q4 chr3 203.2 Mb** | Zm00001eb153000 (CCD8); Zm00001eb153190 (ZCN12) | +0.86 | 2.1 | 0.3 (no) | No significant difference between taxa (Zd +0.30, Zh +0.54, Zl +0.87, Zv +0.60, Zx +1.42). |

Q3 chr3 184.4 Mb (+0.91; *liguleless2*, Zm00001eb147220) is just below the cut (2.0). The chr4 238.3 Mb peak (−log10 p 5.5) lies inside the chr4 197 Mb region and the conditional scan does not add it: given chr4 197 Mb and the other QTL, it is not an independent QTL.

<details><summary>All BL QTL with per-taxon effects and all genes in the interval</summary>

| QTL | Position | Interval (Mb) | Scan −log10 p | Effect per copy (joint) | Joint −log10 p | Taxa differ −log10 p | Per-taxon effects (carrier plots) | Genes in interval |
|---|---|---|---|---|---|---|---|---|
| Q1 (added) | chr2 9.4 Mb | 2.9–11.2 | 9.9 | -1.37 ± 0.25 | **7.3** | 0.0 | Zd -1.68 (32), Zh -1.58 (11), Zl -2.07 (16), Zv -1.61 (54), Zx -1.01 (92) | [438 genes](qtl_genes/BL_qtl_genes.tsv) |
| Q2 | chr2 39.8 Mb | 37.6–41.0 | 12.0 | -1.58 ± 0.26 | **8.8** | 1.8 | Zd -0.83 (30), Zl +1.52 (14), Zv -1.67 (50), Zx -1.92 (91) | [91 genes](qtl_genes/BL_qtl_genes.tsv) |
| Q3 | chr3 184.4 Mb | 181.1–187.0 | 4.7 | +0.91 ± 0.35 | 2.0 | 0.6 | Zd +1.05 (12), Zh -1.46 (6), Zl +0.09 (15), Zv +1.97 (29), Zx +0.71 (65) | [182 genes](qtl_genes/BL_qtl_genes.tsv) |
| Q4 | chr3 203.2 Mb | 195.9–206.5 | 4.5 | +0.86 ± 0.33 | **2.1** | 0.3 | Zd +0.30 (20), Zh +0.54 (9), Zl +0.87 (15), Zv +0.60 (40), Zx +1.42 (55) | [294 genes](qtl_genes/BL_qtl_genes.tsv) |
| Q5 | chr4 197.1 Mb | 193.9–198.7 | 7.5 | -1.74 ± 0.31 | **7.8** | 0.0 | Zd -0.57 (16), Zh -2.69 (6), Zl -3.47 (5), Zv -1.91 (36), Zx -1.75 (69) | [113 genes](qtl_genes/BL_qtl_genes.tsv) |
| Q6 | chr7 185.2 Mb | 183.5–185.8 | 5.2 | -1.16 ± 0.27 | **5.0** | 1.1 | Zd +0.19 (30), Zl +0.60 (8), Zv -1.25 (44), Zx -1.71 (83) | [80 genes](qtl_genes/BL_qtl_genes.tsv) |

Bonferroni over the 6 QTL: −log10 p > 2.08 (bold). Every gene in each interval, with its GO annotation and the k-mer GWAS −log10 p of its window: [`qtl_genes/BL_qtl_genes.tsv`](qtl_genes/BL_qtl_genes.tsv).

</details>

#### 3.6.11 Blade width (BW; 2023 only)

![BW: k-mer ancestry GWAS and joint QTL model](figures/kmerGWAS_joint_qtl_BW.png)

7 regions, 6 hold (Bonferroni over 7: −log10 p > 2.15). **All effects are negative: teosinte alleles make narrower leaf blades**, as expected from teosinte's narrow leaves.

| QTL | Known genes | Effect (joint) | Joint −log10 p | Taxa differ | What it says |
|---|---|---|---|---|---|
| **Q4 chr5 39.2 Mb** | — | -0.31 | 9.4 | 2.0 (no) | No significant difference between taxa (Zd +0.09, Zl -0.16, Zv -0.18, Zx -0.45). |
| **Q5 chr6 146.4 Mb** | — | -0.29 | 8.3 | 0.1 (no) | No significant difference between taxa (Zd -0.25, Zv -0.23, Zx -0.33). |
| **Q3 chr4 237.6 Mb** | — | -0.22 | 5.8 | 0.0 (no) | No significant difference between taxa (Zd -0.22, Zl -0.23, Zv -0.28, Zx -0.20). |
| **Q7 chr9 1.8 Mb** | — | -0.17 | 5.3 | 1.1 (no) | No significant difference between taxa (Zd -0.14, Zh -0.23, Zl +0.08, Zv -0.08, Zx -0.29). The DTA chr9 region. |
| **Q6 chr7 136.2 Mb** | — | -0.17 | 4.3 | 0.3 (no) | No significant difference between taxa (Zd -0.11, Zh -0.16, Zl -0.06, Zv -0.26, Zx -0.14). |
| **Q1 chr4 166.3 Mb** | — | -0.17 | 2.2 | 0.2 (no) | No significant difference between taxa (Zd +0.02, Zv -0.12, Zx -0.18). The EN/Prolif chr4 region. |

The chr6 133.1 Mb peak (*fea4*, Zm00001eb280500; −log10 p 4.85) lies inside the chr6 146 Mb region (131.8–160.5 Mb), and the conditional scan does not add it: it is not independent of chr6 146 Mb.

<details><summary>All BW QTL with per-taxon effects and all genes in the interval</summary>

| QTL | Position | Interval (Mb) | Scan −log10 p | Effect per copy (joint) | Joint −log10 p | Taxa differ −log10 p | Per-taxon effects (carrier plots) | Genes in interval |
|---|---|---|---|---|---|---|---|---|
| Q1 | chr4 166.3 Mb | 165.1–167.9 | 5.9 | -0.17 ± 0.06 | **2.2** | 0.2 | Zd +0.02 (7), Zv -0.12 (34), Zx -0.18 (64) | [70 genes](qtl_genes/BW_qtl_genes.tsv) |
| Q2 | chr4 182.7 Mb | 178.9–183.7 | 4.7 | -0.01 ± 0.06 | 0.1 | 0.6 | Zd -0.07 (13), Zh +0.27 (8), Zl +0.17 (5), Zv +0.05 (38), Zx -0.07 (68) | [137 genes](qtl_genes/BW_qtl_genes.tsv) |
| Q3 | chr4 237.6 Mb | 237.0–238.4 | 8.7 | -0.22 ± 0.05 | **5.8** | 0.0 | Zd -0.22 (19), Zl -0.23 (6), Zv -0.28 (31), Zx -0.20 (63) | [28 genes](qtl_genes/BW_qtl_genes.tsv) |
| Q4 | chr5 39.2 Mb | 34.9–44.1 | 9.1 | -0.31 ± 0.05 | **9.4** | 2.0 | Zd +0.09 (12), Zl -0.16 (13), Zv -0.18 (30), Zx -0.45 (50) | [180 genes](qtl_genes/BW_qtl_genes.tsv) |
| Q5 | chr6 146.4 Mb | 144.6–148.7 | 9.2 | -0.29 ± 0.05 | **8.3** | 0.1 | Zd -0.25 (12), Zv -0.23 (29), Zx -0.33 (50) | [112 genes](qtl_genes/BW_qtl_genes.tsv) |
| Q6 | chr7 136.2 Mb | 130.9–138.6 | 4.5 | -0.17 ± 0.04 | **4.3** | 0.3 | Zd -0.11 (11), Zh -0.16 (8), Zl -0.06 (7), Zv -0.26 (50), Zx -0.14 (62) | [191 genes](qtl_genes/BW_qtl_genes.tsv) |
| Q7 | chr9 1.8 Mb | 0.9–7.1 | 4.5 | -0.17 ± 0.04 | **5.3** | 1.1 | Zd -0.14 (19), Zh -0.23 (9), Zl +0.08 (17), Zv -0.08 (43), Zx -0.29 (82) | [127 genes](qtl_genes/BW_qtl_genes.tsv) |

Bonferroni over the 7 QTL: −log10 p > 2.15 (bold). Every gene in each interval, with its GO annotation and the k-mer GWAS −log10 p of its window: [`qtl_genes/BW_qtl_genes.tsv`](qtl_genes/BW_qtl_genes.tsv).

</details>

#### 3.6.12 Regions shared between traits

| Region | Traits (sign of the teosinte effect) | Known genes in the intervals |
|---|---|---|
| chr10 ~92–106 Mb | DTS +, DTA +, PH +, EH + (just below the cut), NBR + | Zm00001eb418700 (*ZmCCT10*) |
| chr5 ~217–222 Mb | DTS +, DTA +, ASI + (mostly *luxurians*) | Zm00001eb257680 (HDR1-like), Zm00001eb255670 (*mads73*) |
| chr9 ~1–14 Mb | DTA −, DTS − (suggestive), BW − | Zm00001eb372960 (*MADS74*) |
| chr1 ~22–24 Mb | EN +, Prolif +, NBR − | Zm00001eb007670 (*ZmD14*) |
| chr3 ~10–14 Mb | EN +, Prolif +, PH −, EH − | Zm00001eb123060 (*ra2*), Zm00001eb122520 (BZR1-like), Zm00001eb122680 (CYP734A6) |
| chr4 ~165–172 Mb | EN +, Prolif +, BW − | |
| chr4 ~90–140 Mb (conditional QTL) | EN +, Prolif + | |
| chr2 ~6–35 Mb (conditional QTL) | Prolif +, BL − | Zm00001eb077080 (*GA2ox2*), Zm00001eb076990 (*CKX8*) |
| chr7 ~137–144 Mb | PH −, SL −, EN +, Prolif + | |
| chr4 ~194–207 Mb | PH +, EH +, SL +, BL −, ASI − | |
| chr5 ~2–3 Mb | EH +, SL +, PH + (just below the threshold) | |
| chr8 ~126 Mb | EH +, NBR + | Zm00001eb353250 (*ZCN8*) |
| chr7 ~183–186 Mb | EH +, BL − | Zm00001eb331880 (*dlf1*), Zm00001eb331630 (PRR37-like) |
| chr9 ~131–145 Mb | PH +, NBR + | Zm00001eb397240 (*PRR73*) |
| chr2 ~18–30 Mb | ASI −, EH +, EN + | Zm00001eb077080 (*GA2ox2*), Zm00001eb076990 (*CKX8*) |
| chr2 ~185–193 Mb | ASI −, Prolif + | |

All genes in every region, with their GO terms and the −log10 p of the window they lie in, are in `data/introgression_spread/GO/genes_GO_all_<TRAIT>.tsv` (one file per trait).

#### 3.6.13 *mexicana* split by race (SUBTAXON)

The same two analyses were rerun with *mexicana* split into its four races (Chalco, Mesa Central, Durango, Nobogame; section 2.1) instead of one Zx group: the fixed effect, the permutation strata and the per-taxon effects all use the race. Figures: `figures/subtaxon/kmerGWAS_joint_qtl_<TRAIT>.png`; tables: `data/introgression_spread/assoc_gridlmm_spats_subtaxon/` and `data/introgression_spread/joint_qtl/noperm_subtaxon/`.

| Trait | Threshold, taxon / race | QTL, taxon / race | Hold in joint model, taxon / race |
|---|---|---|---|
| DTS | 5.25 / 5.20 | 2 / 2 | 2 / 2 |
| DTA | 5.14 / 5.17 | 3 / 4 | 3 / 3 |
| ASI | 4.35 / 4.42 | 12 / 13 | 6 / 6 |
| PH | 4.19 / 4.26 | 16 / 20 | 8 / 8 |
| EH | 4.27 / 4.24 | 18 / 18 | 7 / 7 |
| EN | 4.57 / 4.47 | 20 / 17 | 11 / 12 |
| Prolif | 4.55 / 4.80 | 15 / 13 | 10 / 9 |
| NBR | 4.21 / 4.06 | 8 / 9 | 4 / 4 |
| SL | 4.06 / 4.14 | 13 / 13 | 12 / 12 |
| BL | 4.14 / 4.32 | 6 / 7 | 5 / 5 |
| BW | 4.22 / 4.26 | 7 / 7 | 6 / 5 |

**What changes:**
- **The QTL are essentially the same.** Most differences are a peak moving by one or two windows (e.g. PH chr7 139.1 → 141.7 Mb, ASI chr5 218.3 → 217.8 Mb) or a broad region being split differently at a slightly different threshold. So splitting *mexicana* does not create or remove signals; treating it as one taxon was not hiding race-specific QTL.
- **DTS / DTA chr5 221.8 Mb:** within *mexicana*, the effect comes from **Chalco** (+0.98 ± 0.41 d for DTS, 8 carrier plots); Mesa Central (+0.25), Durango (+0.09) and Nobogame (+0.06) are near zero. Chalco is the highland race (median 2,390 m). This is the same QTL where *luxurians* and *diploperennis* have the largest effects.
- **DTA chr9 10.9 Mb** (teosinte flowers earlier): Durango and Nobogame −0.66 d, Mesa Central −0.21, Chalco −0.03.
- **ZmCCT10 (chr10 94.5 Mb):** all four races delay flowering (DTS: Chalco +1.06, Durango +1.01, Nobogame +0.93, Mesa Central +0.50 d).
- **Caution:** a race has 6–40 carrier plots at a QTL, and Durango (3 accessions) and Nobogame (2) come from very few donor plants, so a race effect can be one accession's allele.

<details><summary>Figures with the races (DTS, DTA)</summary>

![DTS by race](figures/subtaxon/kmerGWAS_joint_qtl_DTS.png)

![DTA by race](figures/subtaxon/kmerGWAS_joint_qtl_DTA.png)

</details>

#### 3.6.14 Other traits

Scanned with the same pipeline; figures in `figures/kmerGWAS_joint_qtl_<TRAIT>.png`, tables in the same folders as above. SPAD (the first SPAD reading) is not reported; SPAD2 is used for leaf greenness.

| Trait | What it is | Plots | 5% threshold | QTL | Hold in joint model | Strongest QTL (−log10 p) |
|---|---|---|---|---|---|---|
| SPAD2 | leaf greenness (SPAD meter, a proxy for chlorophyll), 36 days after silking; 2023 | 1124 | 4.32 | 6 | 4 | chr1 148.6 Mb (6.0) |
| LAE | number of leaves above the ear; 2025 | 1274 | 4.15 | 8 | 2 | chr2 34.5 Mb (9.3) |
| Kinki | how much the stem bends, score 0–2; 2023 | 1140 | 4.33 | 6 | 4 | chr2 30.1 Mb (12.0) |
| StPu | stem pubescence (present/absent; line mean) | 1286 | 6.96 | 1 | 1 | chr4 240.8 Mb (8.5) |
| StPi | stem pigmentation (present/absent; line mean) | 1286 | 4.77 | 8 | 5 | chr2 19.9 Mb (11.8) |

StPu and StPi are scored present/absent, so their line means are mostly 0 with a few high values, and they agree poorly between years (r = 0.16 and 0.09). With such a skewed trait, a few plots can give very small p-values by chance, which is why the permutation threshold for StPu (6.96) is far above the others; the permutations account for this. A model for binary traits (e.g. a logistic mixed model) would suit them better.

### 3.7 What is new here

- **Local ancestry from k-mers in pooled, ultra-low-coverage samples.** Each sample is a pool of 6 BC2S3 plants sequenced at ~0.8×. The k-mer HMM (Section 2) turns presence and absence of B73 single-copy k-mers into B73 / SEG / DONOR calls per window, with a pool-aware emission model. k-mer GWAS methods use k-mers directly as markers in deeply sequenced individuals; here k-mers are turned into ancestry first.
- **The same ancestry calls give the genotype and the kinship.** The LOCO kinship is built from the ancestry dosage, so relatedness is measured in shared donor tracts, with no SNP calling or imputation.
- **Scan and joint model on the same QTL.** A k-mer ancestry GWAS finds the peaks; a joint QTL model on those peaks gives each QTL's effect adjusted for the others, its effect in each teosinte taxon, and a test that the taxa differ. The dominance test comes from the SEG state.

The statistics themselves (mixed models, LOCO, permutations, R/qtl2) are standard and cited above; what is new is putting k-mer ancestry from pooled skim-seq into them. To be checked against the literature before publication.

### 3.8 Limitations

- **Resolution** is set by tract length: regions are Mb wide, and genes inside them are candidates only.
- **λ_GC** is high under LOCO; the permutation threshold controls the genome-wide error, but linked windows on the same chromosome are not independent evidence, and two nearby peaks may be one QTL.
- **The pool is not the plot.** The genotype is 6 plants of a plot; the phenotype is the plot mean. For SEG plots, plants differ; the dosage of 1 is the expectation.
- **Donor reference.** DONOR levels come from MaizeGDB assemblies, not from our donor plants (Section 2, Limitations).
- **Replication:** peaks should be checked in each year separately and within taxon.
- **Small taxa:** *luxurians* (Zl) and *huehuetenangensis* (Zh) have 5–25 carrier plots at most QTL, so their per-taxon effects have wide CIs. The large *luxurians*-specific effects (DTS/DTA chr5 221.8 Mb) should be checked against how many accessions and F1 families the carriers come from.
- **How peaks are split into regions.** A region is a run of windows above the threshold, joined when less than 5 Mb apart, so one region can hold more than one QTL. The conditional scans of the joint model (section 3.4) test this. They added 6 QTL in the taxon run (BL chr2 9.4; EN chr4 120.2; Prolif chr2 9.3, chr2 34.4, chr4 120.5; StPi chr2 168.1 Mb) and 5 in the race run, and confirmed that other peaks inside regions (e.g. BL chr4 238 Mb, BW chr6 133 Mb, ASI chr3 96 Mb) are not independent QTL.
- **GO terms:** 2,300 of the 18,798 genes in the regions have no GO annotation in the lookup (`data/introgression_spread/GO/genes_missing_GO.txt`); they are listed in the gene tables without GO terms.

---

## Section 4 — References

All references cited in this README, in alphabetical order.

- Aguirre-Liguori JA, Gaut BS, Jaramillo-Correa JP, Tenaillon MI, Montes-Hernández S, García-Oliva F, Hearne SJ, Eguiarte LE (2019) Divergence with gene flow is driven by local adaptation to temperature and soil phosphorus concentration in teosinte subspecies (*Zea mays parviglumis* and *Zea mays mexicana*). *Molecular Ecology* 28:2814–2830. doi:10.1111/mec.15098
- Ahmed HI, et al. (2023) Einkorn genomics sheds light on history of the oldest domesticated wheat. *Nature* 620:830–838. (IBSpy used to find einkorn introgressions in bread wheat.)
- Baum LE, Petrie T, Soules G, Weiss N (1970) A maximization technique occurring in the statistical analysis of probabilistic functions of Markov chains. *Ann Math Stat* 41:164–171.
- Bolger AM, Lohse M, Usadel B (2014) Trimmomatic: a flexible trimmer for Illumina sequence data. *Bioinformatics* 30:2114–2120.
- Broman KW, Gatti DM, Simecek P, et al. (2019) R/qtl2: software for mapping quantitative trait loci with high-dimensional data and multiparent populations. *Genetics* 211:495–502.
- Broman KW, Sen Ś (2009) *A Guide to QTL Mapping with R/qtl*. Springer, New York.
- Brown CT, Irber L (2016) sourmash: a library for MinHash sketching of DNA. *J Open Source Softw* 1:27.
- Chen Q, et al. (2019) TeoNAM: a nested association mapping population for domestication and agronomic trait analysis in maize. *Genetics* 213:1065–1078.
- Churchill GA, Doerge RW (1994) Empirical threshold values for quantitative trait mapping. *Genetics* 138:963–971.
- Danecek P, Bonfield JK, Liddle J, et al. (2021) Twelve years of SAMtools and BCFtools. *GigaScience* 10:giab008.
- Devlin B, Roeder K (1999) Genomic control for association studies. *Biometrics* 55:997–1004.
- Haldane JBS (1919) The combination of linkage values and the calculation of distances between the loci of linked factors. *J Genet* 8:299–309.
- Hrytsenko Y, Daniels NM, Schwartz RS (2025) Determining population structure from k-mer frequencies. *PeerJ* 13:e18939. doi:10.7717/peerj.18939
- Huang C, Sun H, Xu D, et al. (2018) ZmCCT9 enhances maize adaptation to higher latitudes. *PNAS* 115:E334–E341.
- Hung H-Y, Shannon LM, Tian F, et al. (2012) ZmCCT and the genetic basis of day-length adaptation underlying the postdomestication spread of maize. *PNAS* 109:E1913–E1921.
- Irber L, Brooks PT, Reiter T, et al. (2022) Lightweight compositional analysis of metagenomes with FracMinHash and minimum metagenome covers. *bioRxiv* doi:10.1101/2022.01.11.475838
- Ito S, Song YH, Josephson-Day AR, Miller RJ, Breton G, Olmstead RG, Imaizumi T (2012) FLOWERING BHLH transcriptional activators control expression of the photoperiodic flowering regulator CONSTANS in *Arabidopsis*. *PNAS* 109:3582–3587. doi:10.1073/pnas.1118876109
- Kokot M, Długosz M, Deorowicz S (2017) KMC 3: counting and manipulating k-mer statistics. *Bioinformatics* 33:2759–2761.
- Korneliussen TS, Albrechtsen A, Nielsen R (2014) ANGSD: Analysis of Next Generation Sequencing Data. *BMC Bioinformatics* 15:356. doi:10.1186/s12859-014-0356-4
- Lander ES, Waterman MS (1988) Genomic mapping by fingerprinting random clones: a mathematical analysis. *Genomics* 2:231–239. (Coverage formula behind g.)
- Lee JH, Yoo SJ, Park SH, Hwang I, Lee JS, Ahn JH (2007) Role of SVP in the control of flowering time by ambient temperature in *Arabidopsis*. *Genes Dev* 21:397–402.
- Li H (2013) Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. *arXiv* 1303.3997.
- Lippert C, Listgarten J, Liu Y, Kadie CM, Davidson RI, Heckerman D (2011) FaST linear mixed models for genome-wide association studies. *Nat Methods* 8:833–835.
- Listgarten J, Lippert C, Kadie CM, Davidson RI, Eskin E, Heckerman D (2012) Improved linear mixed models for genome-wide association studies. *Nat Methods* 9:525–526.
- Marín-González E, Matías-Hernández L, Aguilar-Jaramillo AE, Lee JH, Ahn JH, Suárez-López P, Pelaz S (2015) SHORT VEGETATIVE PHASE up-regulates TEMPRANILLO2 floral repressor at low ambient temperatures. *Plant Physiol* 169:1214–1224. doi:10.1104/pp.15.00570
- Méndez-Vigo B, Martínez-Zapater JM, Alonso-Blanco C (2013) The flowering repressor SVP underlies a novel *Arabidopsis thaliana* QTL interacting with the genetic background. *PLoS Genet* 9:e1003289. doi:10.1371/journal.pgen.1003289
- Neupane SP, et al. (2025) Genetic basis of Fusarium ear rot resistance and productivity traits in a heterozygous multi-parent recombinant inbred intercross (RIX) maize population. *BMC Plant Biol* 25:639. (MaizeGDB annotation of *bhlh90*.)
- Oakey H, Verbyla A, Pitchford W, Cullis B, Kuchel H (2006) Joint modeling of additive and non-additive genetic line effects in single field trials. *Theor Appl Genet* 113:809–819.
- Odell SG, Hudson AI, Praud S, Dubreuil P, Tixier M-H, Ross-Ibarra J, Runcie DE (2022) Modeling allelic diversity of multiparent mapping populations affects detection of quantitative trait loci. *G3* 12:jkac011.
- Picard toolkit. Broad Institute. https://broadinstitute.github.io/picard/
- Portwood JL 2nd, Woodhouse MR, Cannon EK, et al. (2019) MaizeGDB 2018: the maize multi-genome genetics and genomics database. *Nucleic Acids Res* 47:D1146–D1154. (Source of the teosinte assemblies.)
- Quiroz-Chávez J (2023) Wheat haplotype diversity by a k-mer based approach. PhD thesis, University of East Anglia. (IBSpy.)
- Rabiner LR (1989) A tutorial on hidden Markov models and selected applications in speech recognition. *Proc IEEE* 77:257–286.
- Rahman A, Hallgrímsdóttir I, Eisen M, Pachter L (2018) Association mapping from sequencing reads using k-mers. *eLife* 7:e32920.
- Rodríguez-Álvarez MX, Boer MP, van Eeuwijk FA, Eilers PHC (2018) Correcting for spatial heterogeneity in plant breeding experiments with P-splines. *Spatial Statistics* 23:52–71.
- Rowan BA, Patel V, Weigel D, Schneeberger K (2015) Rapid and inexpensive whole-genome genotyping-by-sequencing for crossover localization and fine-scale genetic mapping. *G3* 5:385–398.
- Runcie DE, Crawford L (2019) Fast and flexible linear mixed models for genome-wide genetics. *PLoS Genet* 15:e1007978.
- sabre: barcode demultiplexing and trimming tool for FastQ files. https://github.com/najoshi/sabre
- Selvanayagam S, Quiroz-Chavez J, Ramirez-Gonzalez RH, Uauy C, Smit S, Schranz ME (2026) KCFtools: rapid alignment-free method for introgression screening and GWAS using k-mer profiles. *Bioinformatics*. doi:10.1093/bioinformatics/btag354
- Smith MW, O'Brien SJ (2005) Mapping by admixture linkage disequilibrium: advances, limitations and guidelines. *Nat Rev Genet* 6:623–632.
- Stitzer MC, Seetharam AS, et al., Hufford MB (2025) Extensive genome evolution distinguishes maize within a stable tribe of grasses. *bioRxiv* doi:10.1101/2025.01.22.633974. (PanAnd assemblies.)
- Sun X, Zhang Z, Wu J, Cui X, Feng D, Wang K, Xu M, Zhou L, Han X, Gu X, Lu T (2016) The *Oryza sativa* regulator HDR1 associates with the kinase OsK4 to control photoperiodic flowering. *PLoS Genet* 12:e1005927. doi:10.1371/journal.pgen.1005927
- Viterbi AJ (1967) Error bounds for convolutional codes and an asymptotically optimum decoding algorithm. *IEEE Trans Inf Theory* 13:260–269.
- Voichek Y, Weigel D (2020) Identifying genetic variants underlying phenotypic variation in plants without complete genomes. *Nat Genet* 52:534–540.
- Wang D, Gao Z, Dong M, et al. (2026) Integrated GWAS and Fst analysis identify flowering-time regulatory genes in maize. *Plant Genome* 19:e70185. doi:10.1002/tpg2.70185
- Yang J, Weedon MN, Purcell S, et al. (2011) Genomic inflation factors under polygenic inheritance. *Eur J Hum Genet* 19:807–812.
- Yang J, Zaitlen NA, Goddard ME, Visscher PM, Price AL (2014) Advantages and pitfalls in the application of mixed-model association methods. *Nat Genet* 46:100–106.
- Yang Q, Li Z, Li W, et al. (2013) CACTA-like transposable element in *ZmCCT* attenuated photoperiod sensitivity and accelerated the postdomestication spread of maize. *PNAS* 110:16969–16974.
- Yin Y, Yan Z, Guan J, Huo Y, Wang T, Li T, Cui Z, Ma W, Wang X, Chen W (2023) Two interacting basic helix-loop-helix transcription factors control flowering time in rice. *Plant Physiol* 192:205–221. doi:10.1093/plphys/kiad077
