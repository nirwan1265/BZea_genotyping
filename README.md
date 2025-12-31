# BZea low-pass genotyping: FASTQ → BAM → GLs (bcftools)

This repository documents the end-to-end processing for **low-pass sequencing** data, starting from raw paired-end FASTQ files and proceeding through:
1) demultiplexing by barcode (sabre)
2) adapter/primer trimming and quality filtering (Trimmomatic)
3) (downstream, required before GL calling) alignment + BAM processing (sort, markdup, index)
4) genotype likelihood generation and genotype calling (bcftools mpileup/call)

> This repo stores **scripts + documentation only** (no large FASTQ/BAM/VCF files).

---

## Project assumptions

- **Compute:** LSF (`bsub`) cluster
- **Shell:** tcsh for historical consistency with existing jobs
- **Data:** paired-end lanes (e.g., `S1_L002` style)
- **Low-pass:** downstream calling uses **genotype likelihoods (GLs)** rather than “hard” genotypes early.

---

## Dependencies

- `sabre` (demultiplexing by barcode)
- `Trimmomatic` (adapter/primer trimming)
- (Downstream) `bwa-mem2` or `bwa`, `samtools`, `picard` (or `samtools markdup`)
- (Calling) `bcftools` (+ reference FASTA indexed with `.fai`)

---

## Directory conventions (recommended)

These paths are examples. On HPC we typically keep data outside the git repo:

- Raw FASTQs (not tracked):  
  `/.../BZea/raw_file/.../BZeaS1_S1_L002_R1_001.fastq.gz` etc.

- Demultiplexed output (not tracked):  
  `/.../BZea/S1_L2/`

- Trimmed output (not tracked):  
  `/.../BZea/filtered_S/filtered_S17/`

This repo contains the **reproducible scripts** to generate those outputs.

---

---

## Step 1 — Demultiplex pooled FASTQs (sabre)

**Goal:** Split pooled lane-level paired-end FASTQs into per-sample FASTQs using a barcode map.

### Inputs
- Pooled lane FASTQs:
  - `..._R1_001.fastq.gz`
  - `..._R2_001.fastq.gz`
- Barcode map file (`.tsv` / `.txt`)

### Barcode file format (required)
Tab-separated **2 columns**:

1. **BARCODE** (exact barcode sequence)
2. **OUTPUT_PREFIX** (sample name / file prefix sabre will write)

Example:

```tsv
ACGTACGT    Sample_001
TGCATGCA    Sample_002
GATCTAGA    Sample_003
```

> sabre will write outputs using the `OUTPUT_PREFIX` names (paired R1/R2 per barcode).  
> Unmatched reads are written to the `no_bc_match_*` files.

### Outputs
- Per-sample paired FASTQs written by sabre (based on column 2 of the barcode file)
- Unmatched reads:
  - `no_bc_match_S<...>_R1.fq.gz`
  - `no_bc_match_S<...>_R2.fq.gz`

### Run (LSF)
```bash
bsub < scripts/01_demultiplex_sabre.tcsh
```

📜 Script: [`scripts/01_demultiplex_sabre.tcsh`](scripts/01_demultiplex_sabre.tcsh)  
📚 Notes: [`docs/01_demultiplex.md`](docs/01_demultiplex.md)

---

## Step 2 — Trim adapters/primers + quality filter (Trimmomatic PE)

**Goal:** Remove adapters/primers, trim low-quality bases, and retain high-quality **paired** reads for alignment.

### Inputs
- Demultiplexed paired FASTQs:
  - `*_R1.fq.gz`
  - `*_R2.fq.gz`

### Outputs
- Paired reads (used for alignment):
  - `*_paired_R1.fq.gz`
  - `*_paired_R2.fq.gz`
- Optional unpaired reads (kept for debugging/QC):
  - `*_unpaired_R1.fq.gz`
  - `*_unpaired_R2.fq.gz`

### Parameters used
- `ILLUMINACLIP:<adapters.fa>:2:30:10`
- `LEADING:3`
- `TRAILING:3`
- `SLIDINGWINDOW:4:15`
- `MINLEN:36`

> Note: If you truly want to discard unpaired reads, redirect them to `/dev/null`.  
> Keeping them is often useful for QC.

### Run (LSF)
```bash
bsub < scripts/02_trim_trimmomatic.tcsh
```

📜 Script: [`scripts/02_trim_trimmomatic.tcsh`](scripts/02_trim_trimmomatic.tcsh)  
📚 Notes: [`docs/02_trimming.md`](docs/02_trimming.md)

---

## Step 3 — BAM preprocessing (sort → mark/remove duplicates → index)

**Goal:** Ensure each sample BAM is **coordinate-sorted**, optionally **deduplicated**, and **indexed** before genotype likelihood (GL) generation with `bcftools mpileup`.

### Why this matters (low-pass)
- `bcftools mpileup` expects **sorted + indexed** BAMs.
- Marking/removing duplicates is typically recommended to reduce PCR/optical duplicate inflation.
- Downstream QC becomes much easier if every BAM follows the same naming convention.

### Inputs
- Aligned BAMs (usually from BWA/BWA-MEM2), one per sample:
  - `*.bam`

> If your BAMs are already sorted, you can skip sorting and go straight to duplicate handling.

### Outputs
For each input `sample.bam`:
- Sorted BAM: `sample.sorted.bam`
- Deduplicated BAM (duplicates removed): `sample.sorted.rmdup.bam`
- Metrics: `sample.dedup_metrics.txt`
- Index: `sample.sorted.rmdup.bam.bai`

### Run (LSF)
```bash
bash scripts/03_bam_sort_dedup_index.sh
```

---

## Step 4 — Genotype likelihoods (GL) genotyping at SNPVersity sites (bcftools)

**Goal:** Use low-pass BAMs to generate genotype likelihood–based calls **only at known SNP sites** (SNPVersity), by running `bcftools mpileup | bcftools call` in parallel across **chromosomes × BAM-chunks**, then merging outputs.

### What you need before running Step 4
- ✅ Reference FASTA + index:  
  - `Zm-B73-REFERENCE-NAM-5.0.fa`  
  - `Zm-B73-REFERENCE-NAM-5.0.fa.fai`  (create with `samtools faidx ref.fa`)
- ✅ SNPVersity per-chromosome VCFs (indexed):  
  - `chr1_high_coverage.vcf.gz` … `chr10_high_coverage.vcf.gz` (+ `.tbi`)
- ✅ BAMs are **sorted + (deduped recommended) + indexed** (`.bai`)
- ✅ Tools available: `bcftools`, `bgzip`, `tabix`, `samtools`

---

### 4.1 Build SNPVersity allele-target files (biallelic SNPs only)
📜 Script: `1_SNPVersity_bialleles.sh`

**Purpose:** Convert each SNPVersity VCF into a **tabix-indexed target file** that contains:  
`CHROM  POS  REF,ALT`

**Outputs (per chromosome):**
- `SNPversity/targets_als/chrN.als.tsv.gz`
- `SNPversity/targets_als/chrN.als.tsv.gz.tbi`

**Run (job array chr1–chr10):**
```bash
bsub < 1_SNPVersity_bialleles.sh
```

---

### 4.2 Create BAM list (absolute paths)
📜 Script: `2_BAM_list.sh`

**Purpose:** Generate a stable, sorted list of all BAMs used for genotyping.

**Output:**
- `GL_work/bamlist.txt`

**Run:**
```bash
bsub < 2_BAM_list.sh
```

---

### 4.3 Split BAM list into chunks (controls open file handles)
📜 Script: `3_create_chunks.sh`

**Purpose:** Split BAM list into manageable chunks (you used `CHUNK=100`) so each mpileup job opens fewer BAMs.

**Outputs:**
- `GL_work/bam_chunks/bamlist.chunk.000`
- `GL_work/bam_chunks/bamlist.chunk.001`
- `...`

**Run:**
```bash
bsub < 3_create_chunks.sh
```

---

### 4.4 Compute GL-based calls per (chromosome × BAM-chunk)
📜 Script: `4_GL.sh`

**Purpose:** For each chromosome and BAM chunk:
- run `bcftools mpileup` restricted to:
  - chromosome: `-r chrN`
  - target sites: `-T chrN.als.tsv.gz`
- then call genotypes with `bcftools call` using `-C alleles` to force SNPVersity alleles

**Outputs:**
- `GL_work/GL_by_chr_chunks/chrN.chunk###.bcf`
- `GL_work/GL_by_chr_chunks/chrN.chunk###.bcf.csi`

**Run (job array; range depends on your total chunks):**
```bash
bsub < 4_GL.sh
```

> Note (recommended for “true GLs”): add PLs in mpileup annotation:  
> `-a FORMAT/DP,FORMAT/AD,FORMAT/PL`  
> Your current script writes DP/AD only.

---

### 4.5 Merge chunk-BCFs into one BCF per chromosome
📜 Script: `5_combine_bcf.sh`

**Purpose:** Each `chrN.chunk###.bcf` contains a **subset of samples**. This step merges all chunk outputs into a single multi-sample chromosome BCF.

**Outputs (per chromosome):**
- `GL_work/GL_by_chr_merged/chrN.merged_withGT.bcf`
- `GL_work/GL_by_chr_merged/chrN.merged_withGT.bcf.csi`

**Run:**
```bash
bsub < 5_combine_bcf.sh
```

⚠️ Important: your current script loops `{10..10}` (only chr10).  
To merge all chromosomes, use `{1..10}`.

---

### 4.6 Concatenate chr1–chr10 into one genome-wide BCF
📜 Script: `6_merge_one_chr.sh`

**Purpose:** Concatenate chromosome BCFs into a single file.

**Outputs:**
- `GL_work/final_genotypes/BZea.chr1_10.bcf`
- `GL_work/final_genotypes/BZea.chr1_10.bcf.csi`

**Run:**
```bash
bsub < 6_merge_one_chr.sh
```

---

## Unfiltered genotype QC summary

![Sample-level genotype statistics for the unfiltered post-calling dataset](figs/all_figs/Fig_genotype_statistics_unfiltered.png)

**Figure 1 | Sample-level genotype statistics for the unfiltered post-calling dataset.**  
(A) Distribution of mean read depth (DP) per sample averaged across variant sites.  
(B) Density of per-sample alternate allele frequency (fraction of called alleles that are non-reference).  
(C) Distribution of per-sample missing genotype rate.  
(D) Distribution of residual heterozygosity per sample, calculated as nHets / nCalled.  
(E) Relationship between residual heterozygosity and alternate allele burden (2×AltHom + Het), computed across called sites.  
(F) Distribution of per-sample transition/transversion (Ts/Tv) ratio.  
This panel set reflects raw calls prior to downstream QC filters; outliers in missingness, heterozygosity, ALT burden, or Ts/Tv flag samples for exclusion or closer inspection. Known controls/checks (e.g., B73 and repeated check lines) may occupy distribution extremes due to reference similarity and/or coverage differences and are assessed separately from the primary study panel.

### Results

Figure 1 summarizes per-sample QC metrics for the *unfiltered* genotype calls (immediately after variant calling, prior to any sample- or site-level filtering) and highlights the expected properties of a raw low-pass dataset along with a small set of clear outliers. Mean depth per sample is narrowly centered around ~1.4–1.7× (Panel A), consistent with uniformly low sequencing depth across the cohort at this stage. The per-sample alternate allele fraction is strongly concentrated at low values with a right-skewed tail (Panel B), indicating that most individuals contribute relatively few non-reference calls while a minority of samples show elevated ALT fractions that merit follow-up (e.g., higher divergence from the reference, contamination/mixture, or mapping artifacts). Missing genotype rate is high and broadly distributed (Panel C), with most samples in the ~60–80% missing range and a tail approaching complete missingness for a subset of individuals, consistent with incomplete site coverage before imputation and before enforcing call-rate thresholds. Residual heterozygosity (nHets / called) is generally low (Panel D) but includes outliers extending to markedly higher values; these same individuals tend to carry a larger overall alternate allele burden (2×AltHom + Het), producing the positive association between heterozygosity and non-reference load (Panel E). The Ts/Tv ratio distribution (Panel F) shows a dominant mode around ~3.4–3.6 with a broader right shoulder, suggesting that most samples share a consistent SNP spectrum while a subset display atypical spectra that often coincide with the missingness and heterozygosity outliers. As with the unimputed QC summaries, known controls/checks can disproportionately populate distribution tails due to reference similarity and/or depth differences, and are interpreted separately when assessing cohort-wide QC thresholds.

---

## Step 5 — Post-calling cleanup (biallelic SNP-only VCF)

**Goal:** Convert the genome-wide BCF into a clean, biallelic SNP-only VCF.gz for downstream **filtering + imputation**.

---

### 5.1 Keep only biallelic SNPs
📜 Script: `7_bialleles.sh`

**Output:**
- `BZea.chr1_10.biallelic_snps.vcf.gz` (+ `.tbi`)

**Run:**
```bash
bsub < 7_bialleles.sh
```

---

### 5.2 Remove sites with empty ALT (recommended fix)
📜 Script: `8_bialleles_remove_empty_alts.sh`

⚠️ Your current `8_*` is identical to `7_*`.  
If the intention is to drop empty/missing ALT alleles, use:

```bash
bcftools view --threads 16 -Oz -v snps -m2 -M2 \
  -e 'ALT="."' \
  -o BZea.chr1_10.biallelic_snps_no_missing_alts.vcf.gz \
  BZea.chr1_10.bcf

tabix -p vcf BZea.chr1_10.biallelic_snps_no_missing_alts.vcf.gz
```

---

---

## Step 6 — Filtering (DP → fill-tags → MAF + missingness)

**Goal:** Starting from the cleaned, biallelic SNP-only VCF, apply:
1) genotype depth (DP) filter (set low-DP genotypes to missing),
2) compute site-level tags (AN/AC/AF/MAF/NS/F_MISSING),
3) filter sites by MAF and missingness to create an imputation-ready VCF.

> Low-pass note: DP thresholds depend on sequencing depth. In this repo we document what was used in the scripts, but you should keep the chosen DP/MAF/missingness thresholds consistent across the run.

---

### 6.1 DP filter (set low-DP genotypes to missing)

📜 Script: `9_filter_DP.sh`

**Input**
- `BZea.chr1_10.biallelic_snps_no_missing_alts.vcf.gz`

**Output**
- `BZea.biallelic_snps.DP5.vcf.gz` (+ `.tbi`)

**Run**
    bsub < 9_filter_DP.sh

**What it does**
- Uses `bcftools filter -S . -e 'FMT/DP<5'`
- Any genotype with DP < 5 becomes `./.` (missing), but the variant record is retained.

---

### 6.2 Add site-level summary tags (AN/AC/AF/MAF/NS/F_MISSING)

📜 Script: `10_add_AF_AC_AN_MAF_NS_F_missing.sh`

**Input**
- `BZea.biallelic_snps.DP5.vcf.gz`

**Output**
- `BZea.biallelic_snps.DP5.tags.vcf.gz` (+ `.tbi`)

**Run**
    bsub < 10_add_AF_AC_AN_MAF_NS_F_missing.sh

---

### 6.3 Filter by MAF and missingness

📜 Script: `11_filter_MAF_F_missing.sh`

**What it does**
- Keeps variants with:
  - `MAF >= 0.005`
  - `F_MISSING <= 0.2`

**Output (from your script)**
- `BZea.DP1.MAF005.MISS20.vcf.gz` (+ `.tbi`)

**Run**
    bsub < 11_filter_MAF_F_missing.sh

⚠️ IMPORTANT consistency check (from your actual scripts)
- Step 6.2 writes: `BZea.biallelic_snps.DP5.tags.vcf.gz`
- Step 6.3 (your script) reads: `BZea.biallelic_snps.DP1.tags.vcf.gz`

So you must either:
- (A) change `11_filter_MAF_F_missing.sh` to read the DP5-tagged file, **or**
- (B) intentionally run a DP1 pipeline and keep naming consistent.

---

## Step 7 — Split by chromosome + imputation (Beagle)

**Goal:** Split the filtered VCF into chr-specific VCFs and run Beagle per chromosome using a genetic map.

---

### 7.1 Split genome-wide VCF into per-chromosome files (chr1–chr10)

📜 Script: `12_separate_chr.sh`

**Input (as in your script)**
- `BZea.DP2.MAF005.MISS50.vcf.gz`

**Outputs**
- `BZea.DP2.MAF005.MISS50.chr1.vcf.gz` (+ `.tbi`)
- …
- `BZea.DP2.MAF005.MISS50.chr10.vcf.gz` (+ `.tbi`)

**Run**
    bsub < 12_separate_chr.sh

⚠️ Your current loop is `{10..10}` (only chr10).
To split all chromosomes, change it to `{1..10}`.

---

### 7.2 Beagle imputation per chromosome (using NAM genetic maps)

📜 Script: `13_beagle_impute.sh`

**Inputs**
- Per-chromosome VCFs from Step 7.1 (e.g. `...chr8.vcf.gz`)
- Beagle jar:
  - `beagle.27Feb25.75f.jar`
- Genetic map per chromosome:
  - `NAM_genetic_map/beagle/chrN.plink.map`

**Outputs**
- `BZea.beagle.chr1.vcf.gz` (+ `.tbi`)
- …
- `BZea.beagle.chr10.vcf.gz` (+ `.tbi`)

**Run**
    bsub < 13_beagle_impute.sh

Notes:
- Your Beagle loop is currently commented (example shows chr8). Uncomment and set `{1..10}` (or the chr range you want).
- After Beagle finishes, index outputs with tabix (your script does this).

---

### 7.3 Add original FORMAT tags back to Beagle output + concatenate all chromosomes

📜 Script: `14_add_tags_back_and_concat.sh`  *(suggested name)*

**Goal:** After Beagle imputation, re-attach useful per-genotype fields (e.g. `AD/DP/PL`) from the **pre-imputation** chr-specific VCFs, then concatenate chr1–chr10 into one imputed genome-wide VCF.

**Inputs (per chromosome)**
- Pre-imputation VCF (source of FORMAT tags):
  - `BZea.DP2.MAF005.MISS50.chr${chr}.vcf.gz`
- Beagle output VCF:
  - `BZea.beagle.chr${chr}.vcf.gz`

**Outputs (per chromosome)**
- Beagle VCF with tags added back:
  - `BZea.beagle.chr${chr}.withPL.vcf.gz` (+ `.tbi`)

**Final output**
- Concatenated genome-wide imputed VCF:
  - `BZea.beagle.imputed.allchr.vcf.gz` (+ `.tbi`)

**Run**
    bsub < 14_add_tags_back_and_concat.sh

**Notes / fixes (based on your script)**
- Your `tabix` line and `concat` line use `withTags`, but your output filename is `withPL`.
  Pick one naming scheme and keep it consistent (recommended: use `withPL` everywhere).

---

