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

## Step 1 — Demultiplex (sabre)

**Goal:** split pooled lane FASTQs into per-sample FASTQs using a barcode table.

**Inputs**
- R1 FASTQ.gz
- R2 FASTQ.gz
- barcode file: `config/barcodes/<S#_L#>_barcodes.tsv`

**Outputs**
- per-barcode paired FASTQs (sabre writes per-barcode files based on barcode table)
- unmatched reads:
  - `no_bc_match_*_R1.fq.gz`
  - `no_bc_match_*_R2.fq.gz`

Run script:
```bash
bsub < scripts/01_demultiplex_sabre.tcsh
