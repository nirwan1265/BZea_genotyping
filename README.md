# BZea Analysis

This repository documents the end-to-end processing for **low-pass sequencing** data, starting from raw paired-end FASTQ files and proceeding through a reproducible genotype + ancestry + association workflow.

## BZea Breeding Scheme:

B73 was crossed with teosinte lines for 2 backcrosses and 3 selfs. (Need more details here)

## Expected genotype proportions after BC2 + 3 selfing generations (per locus)

We calcualte the **expected fraction of genotypes** at any SNP where the recurrent parent (B73) and the donor (teosinte) carry different alleles, assuming:
- B73 is homozygous **BB**
- the teosinte donor is homozygous **TT**
- random mating, no segregation distortion, no selection, and perfect marker informativeness.

This expectation applies **separately within each donor family** we used:
- **Zd** = *Zea diploperennis* (donor allele = T)
- **Zl** = *Zea luxurians* (donor allele = T)
- **Zx** = *Zea mexicana* (donor allele = T)
- **Zv** = *Zea parviglumis* (donor allele = T)

The Mendelian expectations shown below are the same for all families.

---

## Allelic scheme
B73 (BB) × Teo (TT) → **F1 = 100% BT**

Then:
- **BC1 = F1 × BB**
- **BC2 = BC1 × BB**
- then **three selfing generations**: S1 → S2 → S3

---

### Step 1 — After two backcrosses (BC2)

At a locus:

**BC1 (BT × BB):**
- 1/2 BB
- 1/2 BT

**BC2 (BC1 × BB):**
- P(BB) = 3/4
- P(BT) = 1/4
- P(TT) = 0

So at BC2:
- **BB = 0.75**
- **BT = 0.25**
- **TT = 0.00**

This also matches the ancestry expectation. Teosinte allele fraction is **12.5%** at BC2 (because BT loci carry one donor allele copy).

---

### Step 2 — Heterozygosity decay through selfing

After selfing a heterozygote **BT** for *s* generations, heterozygosity halves each generation as:

- **P(BT after s selfs) = (1/2)^s**
- The remaining probability becomes homozygous, split equally:
  - **P(BB) = (1 − (1/2)^s)/2**
  - **P(TT) = (1 − (1/2)^s)/2**

For **s = 3** (three selfs), starting from BT:
- BT = (1/2)^3 = **1/8 = 0.125**
- BB = TT = (1 − 1/8)/2 = **7/16 = 0.4375**


---

### Step 3 — Final genotype proportions after BC2 + 3 selfs

At BC2, only **25%** of loci are heterozygous (BT). The other **75%** are already BB and remain BB through selfing.

So final expectations:

- **P(BB) = 0.75 + 0.25×0.4375 = 0.859375 = 85.94%**
- **P(BT) = 0.25×0.125 = 0.03125 = 3.13%**
- **P(TT) = 0.25×0.4375 = 0.109375 = 10.94%**

- **Genome-wide at informative loci (BC2 → S3):**
- **~85.94% homozygous B73 (BB)**
- **~3.13% heterozygous (BT)**
- **~10.94% homozygous teosinte (TT)**

---

### Implication for low-pass genotyping
With low-pass depth, **heterozygotes are harder to call confidently**, and the expected heterozygote fraction after BC2+S3 is small (~3.1% of loci). Most donor ancestry is expected to appear as **homozygous TT** at loci where donor segments are fixed within a line, which is why imputation / probabilistic ancestry methods are helpful when coverage is low.

---

# General Overview of BZea genotyping, introgression and QTL analysis

## Genotyping

### 1) Demultiplexing by barcode (sabre)
- Split pooled FASTQs into per-sample FASTQs using barcode/index definitions using `sabre`.
- Outputs: `sample_R1.fastq.gz`, `sample_R2.fastq.gz` (+ reports/logs)

### 2) Adapter/primer trimming + quality filtering (Trimmomatic)
- Remove adapter/primer contamination and low-quality bases.
- Outputs: trimmed paired FASTQs (+ unpaired reads, trimming logs)

### 3) Alignment + BAM processing (required before GL calling)
- Align to the **B73 reference** (e.g., B73 Ref v5) using `bwa mem`.
- Post-processing:
  - sort BAM
  - mark duplicates
  - index BAM
  - basic QC (mapping %, coverage summaries)
- Outputs: `sample.sorted.markdup.bam` + `.bai` (+ QC summaries)

### 4) Genotype likelihood generation + genotype calling
Two parallel pipelines:

**A) bcftools pipeline** (for imputation + GWAS)
- Generate genotype likelihoods from low-pass BAMs (`bcftools mpileup`).
- Joint genotype calling across samples (`bcftools call`) to produce cohort VCF/BCF.
- Filtering and normalization
- Outputs: cohort VCF/BCF (raw + filtered)

**B) ANGSD pipeline** (for HMM introgression calling)
- Generate genotype likelihoods using ANGSD with de novo SNP discovery
- Per-sample GL extraction for HMM input
- Outputs: per-sample GL TSV files

### 5) Imputation
- Impute missing genotypes using `Beagle`.
- Outputs: imputed VCF + quality tables

---

## Introgression analysis (HMM-based)

We infer **local ancestry states** along the genome for each line (B73 vs teosinte ancestry) using a custom **Hidden Markov Model (HMM)** with RTIGER-style rigidity, leveraging genotype likelihoods from ANGSD.

**Inputs**
- Per-sample genotype likelihood TSV files (from ANGSD)
- Genetic maps (cM positions per chromosome)

**Outputs**
- Per-marker state calls with posterior probabilities
- Introgression tract BED files
- Summary statistics and chromosome paintings

---

## QTL / association mapping

We perform genome scans using **GridLMM** with **LOCO (leave-one-chromosome-out) kinship** to control relatedness while preserving power on the tested chromosome.

(Details in Section 10)

---

## Software and packages

### Genotyping

- **sabre** — demultiplexing by barcode (FASTQ splitting)
- **Trimmomatic** — adapter/primer trimming + quality filtering
- **bwa-mem2** *(recommended)* or **bwa** — read alignment to reference genome
- **samtools** — BAM/SAM manipulation (sort, index, stats, depth)
- **picard** — duplicate marking
- **bcftools** — genotype likelihoods + variant calling + querying
- **ANGSD** — genotype likelihood estimation + SNP calling for low-pass data
- **vcftools**  — quick VCF filters / summaries
- **bgzip / tabix** *(htslib)* — compress + index VCF/BCF outputs


### Imputation and Variant handling

- **plink2** — variant QC, LD pruning, PCA, relationship matrices, LOCO kinship (`.rel/.rel.id`)
- **Beagle** (Java) — genotype imputation

### Introgression inference (HMM)

- **Python** (>= 3.8)
- **optuna** — Bayesian optimization for parameter tuning
- **R** (>= 4.1)
- **data.table** — R package for post-processing

### QTL / association mapping

- **R**
- **GridLMM** R package
- **plink2** *(for LOCO kinship files)*:

### Gene annotation resources

- **GFF3 annotation** for the reference build (we used v5 annotation - `Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3`)

---

## Section 1 — Demultiplex pooled FASTQs

**Goal:** Split pooled lane-level paired-end FASTQs into per-sample FASTQs using a barcode map.

### Inputs
- Pooled lane FASTQs:
  - `..._R1_001.fastq.gz`
  - `..._R2_001.fastq.gz`
- Barcode map file (`.txt`)

### Barcode file format
Tab-separated **2 columns**:

1. **BARCODE** (exact barcode sequence)
2. **OUTPUT_PREFIX** (sample name)

Example:

```txt
ACGTACGT    Sample_001
TGCATGCA    Sample_002
GATCTAGA    Sample_003
```

Unmatched reads are written to the `no_bc_match_*` files.

### Outputs
- Per-sample paired FASTQs written by sabre (based on column 2 of the barcode file)
- Unmatched reads:
  - `no_bc_match_S<...>_R1.fq.gz`
  - `no_bc_match_S<...>_R2.fq.gz`

### Run
```bash
bsub < scripts/01_demultiplex_sabre.tcsh
```

Script: [`scripts/01_demultiplex_sabre.sh`](scripts/01_demultiplex_sabre.sh)

---

## Section 2 — Trim adapters/primers + quality filter (Trimmomatic PE)

**Goal:** Remove adapters/primers

### Inputs
- Demultiplexed paired FASTQs:
  - `*_R1.fq.gz`
  - `*_R2.fq.gz`

### Outputs
- Paired reads (used for alignment):
  - `*_paired_R1.fq.gz`
  - `*_paired_R2.fq.gz`

### Parameters used
- `ILLUMINACLIP:<adapters.fa>:2:30:10`
- `LEADING:3`
- `TRAILING:3`
- `SLIDINGWINDOW:4:15`
- `MINLEN:36`


### Run
```bash
bsub < scripts/02_trim_trimmomatic.tcsh
```

Script: [`scripts/02_trim_trimmomatic.sh`](scripts/02_trim_trimmomatic.sh)

---

## Section 3 — BAM preprocessing (sort → mark/remove duplicates → index)

**Goal:** Ensure each sample BAM is **coordinate-sorted**, **deduplicated**, and **indexed** before genotype likelihood (GL) generation with `bcftools mpileup`.

Ensure the files are:
- sorted
- Mark/remove duplicates
- Index BAMs files

### Inputs
- Aligned BAMs (usually from BWA/BWA-MEM2), one per sample:
  - `PN*_SID*.bam`

### Outputs
For each input `PN*_SID*.bam`
- Sorted BAM: `PN*_SID*.sorted.bam`
- Deduplicated BAM (duplicates removed): `PN*_SID*.sorted.rmdup.bam`
- Metrics: `PN*_SID*.dedup_metrics.txt`
- Index: `PN*_SID*.sorted.rmdup.bam.bai`

### Run
```bash
bash scripts/03_bam_sort_dedup_index.sh
```

Script: [`scripts/03_bam_sort_dedup_index.sh`](scripts/03_bam_sort_dedup_index.sh)

---

## Section 4 — ANGSD Genotype Likelihood Calling (for HMM Introgression)

**Goal:** Generate per-sample genotype likelihoods using ANGSD with de novo SNP discovery, optimized for downstream HMM-based introgression calling.

### Why ANGSD for introgression analysis?

ANGSD is particularly well-suited for low-pass sequencing data because it:
- Performs probabilistic SNP calling without requiring a known SNP panel
- Outputs genotype likelihoods (GL) that preserve uncertainty
- Handles low-depth data gracefully by not forcing hard genotype calls

---

### 4.1 ANGSD SNP calling + GL generation (per chromosome, parallelized)

**Goal:** Run ANGSD on each chromosome, split into 20 parts for parallelization.

**Key parameters:**
- `-GL 1`: SAMtools GL model
- `-doMajorMinor 1`: Infer major/minor alleles from GL
- `-doMaf 1`: Estimate allele frequencies
- `-doBcf 1`: Output BCF with GL fields
- `-SNP_pval 1e-8`: SNP calling p-value threshold
- `-minMapQ 30`: Minimum mapping quality
- `-minQ 20`: Minimum base quality
- `-minInd 5`: Minimum individuals with data at a site
- `-skipTriallelic 1`: Keep only biallelic sites
- `-setMinDepthInd 2`: Minimum depth per individual

**Outputs:**
- `BZea_chrN.partM.START_END.bcf` (per chromosome part)

**Run:**
```bash
bsub < scripts/HMM_introgression/0.0_angsd_genotyping.sh
```

Script: [`scripts/HMM_introgression/0.0_angsd_genotyping.sh`](scripts/HMM_introgression/0.0_angsd_genotyping.sh)

---

### 4.2 Handle failed jobs (optional)

**Goal:** Identify and re-run any failed ANGSD jobs.

**Run:**
```bash
bash scripts/HMM_introgression/0.1_find_failed_runs.sh
bsub < scripts/HMM_introgression/0.2_run_failed_parts.sh
```

Scripts:
- [`scripts/HMM_introgression/0.1_find_failed_runs.sh`](scripts/HMM_introgression/0.1_find_failed_runs.sh)
- [`scripts/HMM_introgression/0.2_run_failed_parts.sh`](scripts/HMM_introgression/0.2_run_failed_parts.sh)

---

### 4.3 Combine sub-parts (if jobs were split further)

**Goal:** Merge any sub-parts created from failed job re-runs.

**Outputs:**
- Merged BCF files with correct coordinate ordering

**Run:**
```bash
bsub < scripts/HMM_introgression/0.3_combine_sub_parts.sh
```

Script: [`scripts/HMM_introgression/0.3_combine_sub_parts.sh`](scripts/HMM_introgression/0.3_combine_sub_parts.sh)

---

### 4.4 Concatenate parts into per-chromosome BCFs

**Goal:** Merge all 20 parts for each chromosome into a single BCF.

**Outputs:**
- `per_chr/BZea_chrN.full.bcf` (+ `.csi` index)

**Run:**
```bash
bsub < scripts/HMM_introgression/0.4_combine_parts_chrom.sh
```

Script: [`scripts/HMM_introgression/0.4_combine_parts_chrom.sh`](scripts/HMM_introgression/0.4_combine_parts_chrom.sh)

---

### 4.5 Rename sample headers (BAM paths → sample IDs)

**Goal:** Replace full BAM paths in BCF headers with clean sample IDs (e.g., `PN9_SID857`).

**Outputs:**
- `per_chr_renamed/BZea_chrN.full.renamed.bcf`

**Run:**
```bash
bsub < scripts/HMM_introgression/0.5_rename_headers_files.sh
```

Script: [`scripts/HMM_introgression/0.5_rename_headers_files.sh`](scripts/HMM_introgression/0.5_rename_headers_files.sh)

---

### 4.6 Filter BCFs for HMM input

**Goal:** Apply quality filters to create a clean SNP set for introgression calling.

**Filters applied:**
- `QUAL >= 20`: Minimum variant quality
- `F_MISSING <= 0.2`: Maximum 20% missing data
- Biallelic SNPs only
- `AF >= 0.05 AND AF <= 0.95`: Remove fixed/nearly-fixed sites
- `NS >= 100`: Minimum number of samples with data

**Outputs:**
- `filtered/BZea_chrN.filtered.bcf`

**Run:**
```bash
bash scripts/HMM_introgression/0.6_filter_bcf.sh
```

Script: [`scripts/HMM_introgression/0.6_filter_bcf.sh`](scripts/HMM_introgression/0.6_filter_bcf.sh)

---

### 4.7 Extract per-sample GL for HMM input

**Goal:** Extract genotype likelihoods for each sample into a TSV format suitable for the HMM.

**Output format:**
```
CHROM   POS     REF     ALT     glRR    glRH    glHH
chr1    12345   A       G       -0.5    -2.3    -15.1
```

Where `glRR`, `glRH`, `glHH` are log10 genotype likelihoods for RR (Ref/Ref), RH (Ref/Het), HH (Alt/Alt).

**Outputs:**
- `sample.chr1-10.GL.tsv.gz` (per sample, all chromosomes)

**Run:**
```bash
bsub < scripts/HMM_introgression/1_extract_GL_emissions.sh
```

Script: [`scripts/HMM_introgression/1_extract_GL_emissions.sh`](scripts/HMM_introgression/1_extract_GL_emissions.sh)

---

## Section 5 — bcftools Genotype Calling at SNPVersity Sites

**Goal:** Generate genotype likelihood–based calls **only at known SNP sites** using SNPVersity, by running `bcftools mpileup | bcftools call` in parallel across **chromosomes × BAM-chunks**, then merging outputs.

(This section is used for imputation and GWAS, parallel to the ANGSD pipeline above)

### What you need before running Step 5
- Reference FASTA + index:
  - `Zm-B73-REFERENCE-NAM-5.0.fa`
  - `Zm-B73-REFERENCE-NAM-5.0.fa.fai`  (create with `samtools faidx ref.fa`)
- SNPVersity per-chromosome VCFs (indexed):
  - `chr1_high_coverage.vcf.gz` … `chr10_high_coverage.vcf.gz` (+ `.tbi`)
- BAMs are **sorted + (deduped recommended) + indexed** (`.bai`)

---

### 5.1 Build SNPVersity allele-target files (biallelic SNPs only)

**Goal:** Convert each SNPVersity VCF into a **tabix-indexed target file** that contains:
`CHROM  POS  REF,ALT`

**Outputs (per chromosome):**
- `SNPversity/targets_als/chrN.als.tsv.gz`
- `SNPversity/targets_als/chrN.als.tsv.gz.tbi`

**Run**
```bash
bsub < 04_SNPVersity_bialleles.sh
```
Script: `04_SNPVersity_bialleles.sh`

---

### 5.2 Create BAM list (absolute paths)

**Goal:** Generate a sorted list of all BAMs used for genotyping.

**Output:**
- `GL_work/bamlist.txt`

**Run:**
```bash
bsub < 05_BAM_list.sh
```
Script: `05_BAM_list.sh`

---

### 5.3 Split BAM list into chunks (controls open file handles)

**Goal:** Split BAM list into manageable chunks (you used `CHUNK=100`) so each mpileup job opens fewer BAMs.

**Outputs:**
- `GL_work/bam_chunks/bamlist.chunk.000`
- `GL_work/bam_chunks/bamlist.chunk.001`
- `...`

**Run:**
```bash
bsub < 05_create_chunks.sh
```
Script: `05_create_chunks.sh`

---

### 5.4 Compute GL-based calls per (chromosome × BAM-chunk)

**Goal:** For each chromosome and BAM chunk:
- run `bcftools mpileup` restricted to:
  - chromosome: `-r chrN`
  - target sites: `-T chrN.als.tsv.gz`
- then call genotypes with `bcftools call` using `-C alleles` to force SNPVersity alleles

**Outputs:**
- `GL_work/GL_by_chr_chunks/chrN.chunk###.bcf`
- `GL_work/GL_by_chr_chunks/chrN.chunk###.bcf.csi`

**Run**
```bash
bsub < 06_GL.sh
```
Script: `06_GL.sh`

---

### 5.5 Merge chunk-BCFs into one BCF per chromosome

**Goal:** Each `chrN.chunk###.bcf` contains a **subset of samples**. This step merges all chunk outputs into a single multi-sample chromosome BCF.

**Outputs (per chromosome):**
- `GL_work/GL_by_chr_merged/chrN.merged_withGT.bcf`
- `GL_work/GL_by_chr_merged/chrN.merged_withGT.bcf.csi`

**Run:**
```bash
bsub < 07_combine_bcf.sh
```
Script: `07_combine_bcf.sh`

---

### 5.6 Concatenate chr1–chr10 into one genome-wide BCF

**Goal:** Concatenate chromosome BCFs into a single file.

**Outputs:**
- `GL_work/final_genotypes/BZea.chr1_10.bcf`
- `GL_work/final_genotypes/BZea.chr1_10.bcf.csi`

**Run:**
```bash
bsub < 08_merge_one_chr.sh
```
Script: `08_merge_one_chr.sh`

---

## Section 6 — Post-calling cleanup (biallelic SNP-only VCF)

**Goal:** Convert the genome-wide BCF into a biallelic SNP-only VCF.gz for downstream **filtering + imputation**.

### 6.1 Keep only biallelic SNPs

**Output:**
- `BZea.chr1_10.biallelic_snps.vcf.gz` (+ `.tbi`)

**Run:**
```bash
bsub < 7_bialleles.sh
```
Script: `09_bialleles.sh`

---

### 6.2 Remove sites with empty ALT (recommended fix)

```bash
bsub < 10_bialleles_remove_empty_alts.sh
```
Script: `10_bialleles_remove_empty_alts.sh`

---

## Section 7 — Filtering (DP → fill-tags → MAF + missingness)

**Goal:** Starting from the biallelic SNP-only VCF, apply:
1) genotype depth (DP) filter (set low-DP genotypes to missing),
2) compute site-level tags (AN/AC/AF/MAF/NS/F_MISSING),
3) filter sites by MAF and missingness to create an imputation-ready VCF.

For PCA we used the filters described below. For RTIGER introgression analysis, we used a more relaxed filtering parameters.

### 7.1 DP filter (set low-DP genotypes to missing)

**Input**
- `BZea.chr1_10.biallelic_snps_no_missing_alts.vcf.gz`

**Parameters**
- Uses `bcftools filter -S . -e 'FMT/DP<5'`
- Any genotype with DP < 5 becomes `./.` (missing), but the variant record is retained.

**Output**
- `BZea.biallelic_snps.DP5.vcf.gz` (+ `.tbi`)

**Run**
```bash
bsub < 11_filter_DP.sh
```
Script: `11_filter_DP.sh`

---

### 7.2 Add site-level summary tags (AN/AC/AF/MAF/NS/F_MISSING)

**Input**
- `BZea.biallelic_snps.DP5.vcf.gz`

**Output**
- `BZea.biallelic_snps.DP5.tags.vcf.gz` (+ `.tbi`)

**Run**
```bash
bsub < 12_add_AF_AC_AN_MAF_NS_F_missing.sh
```
Script: `12_add_AF_AC_AN_MAF_NS_F_missing.sh`

---

### 7.3 Filter by MAF and missingness

**Goal** Keeps variants with `MAF >= 0.005` and `F_MISSING <= 0.5` for PCA and imputation

**Output**
- `BZea.DP1.MAF005.MISS50.vcf.gz` (+ `.tbi`)

**Run**
```bash
bsub < 13_filter_MAF_F_missing.sh
```
Script: `13_filter_MAF_F_missing.sh`

---

## Section 8 — Split by chromosome + imputation (Beagle)

**Goal:** Split the filtered VCF into chr-specific VCFs and run Beagle imputation per chromosome using a genetic map.


### 8.1 Split genome-wide VCF into per-chromosome files (chr1–chr10)

**Input**
- `BZea.DP2.MAF005.MISS50.vcf.gz`

**Outputs**
- `BZea.DP2.MAF005.MISS50.chr1.vcf.gz` (+ `.tbi`)
- …
- `BZea.DP2.MAF005.MISS50.chr10.vcf.gz` (+ `.tbi`)

**Run**
```bash
    bsub < 14_separate_chr.sh
```
Script: `scripts/14_separate_chr.sh`

---

### 8.2 Beagle imputation per chromosome (using NAM genetic maps)


**Inputs**
- Per-chromosome VCFs from Step 8.1
- Beagle software:
  - `beagle.27Feb25.75f.jar`
- Genetic map per chromosome:
  - `NAM_genetic_map/beagle/chrN.plink.map` or anything similar

**Outputs**
- `BZea.beagle.chr1.vcf.gz` (+ `.tbi`)
- …
- `BZea.beagle.chr10.vcf.gz` (+ `.tbi`)

**Run**
```bash
bsub < 15_beagle_impute.sh
```
Script: `scripts/15_beagle_impute.sh`

---

### 8.3 Add original FORMAT tags back to Beagle output + concatenate all chromosomes

**Goal:** After Beagle imputation, re-attach useful per-genotype fields (e.g. `AD/DP/PL`) from the **pre-imputation** chr-specific VCFs, then concatenate chr1–chr10 into one imputed genome-wide VCF.

**Inputs**
- Pre-imputation VCF:
  - `BZea.DP2.MAF005.MISS50.chr${chr}.vcf.gz`
- Beagle output VCF:
  - `BZea.beagle.chr${chr}.vcf.gz`

**Outputs**
- Beagle VCF with tags added back:
  - `BZea.beagle.chr${chr}.withPL.vcf.gz` (+ `.tbi`)
- Concatenated genome-wide imputed VCF:
  - `BZea.beagle.imputed.allchr.vcf.gz` (+ `.tbi`)

**Run**
```bash
    bsub < 16_add_tags_back_and_concat.sh
```
Script: `scripts/16_add_tags_back_and_concat.sh`

---

## Section 9 — HMM-based Introgression Analysis

We infer **local ancestry states** along the genome for each BC2S3 line using a custom Hidden Markov Model (HMM) with RTIGER-style rigidity constraints. This approach uses genotype likelihoods (GLs) rather than hard genotype calls, preserving uncertainty from low-pass sequencing.

---

### 9.1 HMM Model Overview

**Goal:** Call local ancestry states (RR, RH, HH) across the genome.

| State | Meaning | Genotype |
|-------|---------|----------|
| RR | Homozygous Reference (B73/B73) | AA |
| RH | Heterozygous (B73/Teosinte) | AB |
| HH | Homozygous Teosinte (Teo/Teo) | BB |

**Key features:**
- Uses genotype likelihoods (GL) as emissions, not hard calls
- Genetic map-based transition probabilities (Haldane function)
- RTIGER-style rigidity: requires R consecutive markers supporting a state change
- Prevents noisy single-SNP state switches

---

### 9.2 HMM Parameters

| Parameter | Description | Optimized Value |
|-----------|-------------|-----------------|
| `--rigidity` / `-R` | Consecutive markers needed to switch | **41** |
| `--rho` | Transition stickiness multiplier | **0.023** |
| `--eta_hh_from_rh` | HH rescue from RH emissions | **0.012** |
| `--prior_rr` | Prior probability of RR state | 0.85 |
| `--prior_rh` | Prior probability of RH state | 0.05 |
| `--prior_hh` | Prior probability of HH state | 0.10 |
| `--decode` | Decoding method | `posterior_hysteresis` |

---

### 9.3 Parameter Optimization (Stability-based)

Since we have **no labeled ground truth**, we optimize parameters by maximizing **stability under marker thinning**:

1. Run HMM on full data → baseline results
2. Create K thinned replicates (drop 10% markers randomly)
3. Run HMM on each replicate
4. Score = similarity between baseline and replicates
   - Jaccard similarity of donor blocks
   - State concordance at shared markers
   - Fragmentation penalty

**Run optimization:**
```bash
cd scripts/HMM_introgression

export HMM_PY="$(pwd)/hmm_viterbi_bc2s3.py"
export MAPDIR="$(pwd)/genetic_maps"
export POST_R="$(pwd)/batch_introgression_analysis.R"
export OUTROOT="$(pwd)/optuna_tuning_runs"
export MAX_GAP_KB=10000
export MIN_BLOCK_KB=5000
export DROP_FRAC=0.10
export K_REPS=3

conda activate hmm
python tune_hmm_optuna.py calib_samples.txt 100
```

**Outputs:**
- `optuna_tuning_runs/best_params.txt`
- `optuna_tuning_runs/optuna_study.db`

Scripts:
- [`scripts/HMM_introgression/tune_hmm_optuna.py`](scripts/HMM_introgression/tune_hmm_optuna.py)

---

### 9.4 Run HMM per sample

**Run:**
```bash
bsub < scripts/HMM_introgression/2_run_hmm_bySample.sh
```

**Command (single sample):**
```bash
python hmm_viterbi_bc2s3.py \
  --in_gl_tsv_gz sample.GL.tsv.gz \
  --map_dir genetic_maps/ \
  --out_statepath_gz sample.statepath.tsv.gz \
  --out_tracts_bed_gz sample.tracts.bed.gz \
  --prior_rr 0.85 --prior_rh 0.05 --prior_hh 0.10 \
  --min_morgan 1e-12 \
  --eta_hh_from_rh 0.012 \
  --rh_penalty 0 \
  --decode posterior_hysteresis \
  --rho 0.023 \
  --rigidity 41 \
  --min_run_hh 3 \
  --min_run_rh 5
```

**Outputs per sample:**
1. `sample.statepath.tsv.gz` - Per-marker state calls + posteriors
   ```
   CHROM  POS    REF  ALT  DP  STATE  P_RR   P_RH   P_HH   THETA
   chr1   12345  A    G    15  RR     0.95   0.04   0.01   0.0001
   ```

2. `sample.tracts.bed.gz` - Merged segments
   ```
   chrom  start    end      state
   chr1   0        5000000  RR
   chr1   5000001  8000000  RH
   ```

Scripts:
- [`scripts/HMM_introgression/hmm_viterbi_bc2s3.py`](scripts/HMM_introgression/hmm_viterbi_bc2s3.py)
- [`scripts/HMM_introgression/2_run_hmm_bySample.sh`](scripts/HMM_introgression/2_run_hmm_bySample.sh)

---

### 9.5 Post-analysis: Block merging + visualization

**Goal:**
- Merge adjacent same-state blocks
- Bridge small gaps between same-state blocks
- Filter out tiny blocks
- Generate summary statistics and chromosome paintings

**Parameters:**
| Parameter | Description | Value |
|-----------|-------------|-------|
| `max_gap_kb` | Max gap to bridge between same-state blocks | 10000 kb |
| `min_block_kb` | Minimum block size to keep | 5000 kb |

**Run:**
```bash
Rscript batch_introgression_analysis.R input_folder/ 10000 5000 output_folder/
```

**Outputs per sample:**
- `sample.complete_blocks.bed` - Final introgression blocks
- `sample.complete_blocks.png` - Chromosome painting visualization
- `sample.block_summary.tsv` - Block counts by genotype
- `sample.chr_summary.tsv` - Per-chromosome statistics

**Aggregate outputs:**
- `introgression_summary_total.tsv` - Ref%, Het%, Teo%, Donor% per sample
- `introgression_summary_per_chr_donor.tsv` - Donor% by chromosome
- `introgression_summary_per_chr_het.tsv` - Het% by chromosome
- `introgression_summary_per_chr_teo.tsv` - Teo% by chromosome

Scripts:
- [`scripts/HMM_introgression/batch_introgression_analysis.R`](scripts/HMM_introgression/batch_introgression_analysis.R)
- [`scripts/HMM_introgression/plot_introgression_blocks.R`](scripts/HMM_introgression/plot_introgression_blocks.R)

---

### 9.6 Chromosome painting color scheme

- **Blue (AA)** = B73/B73 (Reference)
- **Purple (AB)** = Heterozygous
- **Red (BB)** = Teo/Teo (Donor)

---

### 9.7 Complete HMM workflow

```
┌─────────────────────────────────────────────────────────────┐
│  1. ANGSD GL CALLING (Section 4)                            │
│     0.0 - 0.6 scripts → per-sample GL TSV files             │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│  2. PARAMETER OPTIMIZATION (once)                           │
│     tune_hmm_optuna.py → best_params.txt                    │
│     (R=41, rho=0.023, eta=0.012)                            │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│  3. RUN HMM (per sample)                                    │
│     2_run_hmm_bySample.sh                                   │
│     → sample.statepath.tsv.gz                               │
│     → sample.tracts.bed.gz                                  │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│  4. POST-ANALYSIS                                           │
│     batch_introgression_analysis.R                          │
│     → sample.complete_blocks.bed                            │
│     → sample.complete_blocks.png                            │
│     → introgression_summary_total.tsv                       │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│  5. DOWNSTREAM ANALYSIS                                     │
│     - GWAS on introgression bins                            │
│     - Population-level introgression patterns               │
│     - QTL mapping                                           │
└─────────────────────────────────────────────────────────────┘
```

---

### 9.8 HMM File Reference

| File | Purpose |
|------|---------|
| `hmm_viterbi_bc2s3.py` | Main HMM script |
| `batch_introgression_analysis.R` | Post-processing + plotting |
| `tune_hmm_optuna.py` | Parameter optimization |
| `split_genetic_map.py` | Split NAM map into chr*.map |
| `genetic_maps/chr*.map` | Per-chromosome genetic maps |
| `calib_samples.txt` | Sample list for optimization |
| `best_params.txt` | Optimized parameters |

For detailed documentation, see: [`scripts/HMM_introgression/README.md`](scripts/HMM_introgression/README.md)

---

## Section 10 — QTL / Association Mapping

(To be added)

---
