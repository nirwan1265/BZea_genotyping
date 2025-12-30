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



