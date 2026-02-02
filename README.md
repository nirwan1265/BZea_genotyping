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

### 4) Genotype likelihood generation + genotype calling (ANGSD)
- Generate genotype likelihoods using ANGSD with de novo SNP discovery
- Per-sample GL extraction for HMM input
- Outputs: per-sample GL TSV files

### 5) Imputation
- Impute missing genotypes using `Beagle`.
- Outputs: imputed VCF + quality tables

## HMM introgression analysis

### 6) HMM introgression analysis
We infer **local ancestry states** along the genome for each line (B73 vs teosinte ancestry) using a custom **Hidden Markov Model (HMM)** with RTIGER-style rigidity, leveraging genotype likelihoods from ANGSD.

**Inputs**
- Per-sample genotype likelihood TSV files (from ANGSD)
- Genetic maps (cM positions per chromosome)

**Outputs**
- Per-marker state calls with posterior probabilities
- Introgression tract BED files
- Summary statistics and chromosome paintings


##  QTL / Association Mapping
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
- **bcftools** — querying
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


## Section 5 — Imputation

---

## Section 6 — HMM-based Introgression Analysis

We infer **local ancestry states** along the genome for each BC2S3 line using a custom 3-state Hidden Markov Model (HMM) with **RTIGER-style rigidity constraints**. The model operates directly on **genotype likelihoods (GLs)** rather than hard genotype calls or allele counts, preserving uncertainty typical of low-pass sequencing.

---

### 6.1 HMM Model Overview

**Goal:** Call local ancestry states (RR, RH, HH) across the genome.

| State | Meaning | Genotype label |
|-------|---------|----------------|
| RR | Homozygous Reference (B73/B73) | AA |
| RH | Heterozygous (B73/Teosinte) | AB |
| HH | Homozygous Donor (Teo/Teo) | BB |

**Key features**
- Uses **genotype likelihoods (GL)** as emissions (not allele counts, not hard calls)
- **Genetic map–based transitions** via the **Haldane mapping function**
- **Sticky transition model** controlled by a global scaling parameter (`rho`)
- **RTIGER-style rigidity (hysteresis):** requires **R consecutive markers** supporting a state change before switching which can be calculated through Bayesian optimization
- Reduces noisy single-marker switches and yields coherent introgressed tracts

---

### 6.2 HMM Model Specification (paper-ready detail)

Let markers be indexed by $t = 1, \ldots, T$, ordered by chromosome and position. The hidden state at marker $t$ is

$$
z_t \in \{\mathrm{RR},\mathrm{RH},\mathrm{HH}\}.
$$

The observed data at marker $t$ is the vector of genotype likelihoods:

$$
x_t = (GL_t(\mathrm{RR}), GL_t(\mathrm{RH}), GL_t(\mathrm{HH})),
$$

provided on the **log10** scale. These are converted to natural log scale:

$$
\ell_t(s) = GL_t(s)\cdot \ln(10).
$$

#### 6.2.1 Emission model (GL-based)
The emission log-likelihoods are:

$$
\ell_t(\mathrm{RR}) \; \ell_t(\mathrm{RH}) \; \ell_t(\mathrm{HH}).
$$

To improve robustness (especially at low depth), we apply optional emission adjustments:

**HH rescue from RH (`eta_hh_from_rh`)**  
We allow HH to “borrow” likelihood mass from RH:

$$
\ell'_t(\mathrm{HH}) = \ln\left[(1-\eta)\exp(\ell_t(\mathrm{HH}))+\eta\exp(\ell_t(\mathrm{RH}))\right],
$$

where $\eta$ = `eta_hh_from_rh`.

**Optional RR borrow from RH (`eta_rr_from_rh`, default 0)**  
Analogous mixture for RR (typically left at 0):

$$
\ell'_t(\mathrm{RR}) = \ln\left[(1-\eta_{rr})\exp(\ell_t(\mathrm{RR}))+\eta_{rr}\exp(\ell_t(\mathrm{RH}))\right].
$$

**RH penalty (`rh_penalty`)**  
We can downweight RH to reduce spurious heterozygote calls:

$$
\ell'_t(\mathrm{RH}) = \ell_t(\mathrm{RH}) - \lambda,
$$

where $\lambda$ = `rh_penalty`.

The adjusted emissions $\ell'_t(\cdot)$ are used throughout inference and decoding.

---

#### 6.2.2 Transition model (genetic map + sticky switching)
Each chromosome has a genetic map providing $(bp \rightarrow cM)$. We interpolate each marker to a genetic position $cM_t$, then compute genetic distance between consecutive markers:

$$
\Delta cM_t = |cM_t - cM_{t-1}|.
$$

Convert to Morgans:

$$
d_t = \max(\Delta cM_t / 100,\; \text{min\_morgan}).
$$

Compute recombination fraction via Haldane:

$$
\theta_t = 0.5\left(1-\exp(-2d_t)\right).
$$

We then define a **sticky** switching probability:

$$
s_t = \rho\cdot \theta_t,
$$

with $\rho$ = `rho` controlling global switching propensity. For numerical stability, $s_t$ is bounded (implementation detail):

$$
s_t \leftarrow \min(0.25,\; \max(10^{-12}, s_t)).
$$

Let $\pi = (\pi_{RR},\pi_{RH},\pi_{HH})$ be the stationary state probabilities derived from user priors:

$$
\pi_s = \frac{\text{prior\_s}}{\sum_{s'}\text{prior\_s'}}.
$$

The per-step transition matrix $A_t$ is defined as:
- Self-transition:

$$
A_t(i\rightarrow i) = 1 - s_t
$$

- Off-diagonals distributed proportional to $\pi$ (excluding self):

$$
A_t(i\rightarrow j) = s_t \cdot \frac{\pi_j}{1-\pi_i},\quad j\neq i.
$$

This favors persistence but allows switching in proportion to genetic distance and prior expectations.

---

#### 6.2.3 Forward–backward posteriors
We compute posterior state probabilities:

$$
P(z_t=s \mid x_{1:T})
$$

using the forward–backward algorithm in log space (with per-position normalization to prevent underflow). These posteriors are emitted in the **statepath** output as $P_{RR}, P_{RH}, P_{HH}$ per marker.

---

#### 6.2.4 Decoding modes (including Viterbi)
The implementation supports three decoding strategies (`--decode`):

**(A) Viterbi (`viterbi`)**  
Find the maximum a posteriori path:

$$
\hat{z}_{1:T} = \arg\max_{z_{1:T}}\left[\ln \pi_{z_1}+\sum_{t=1}^T \ell'_t(z_t)+\sum_{t=2}^T \ln A_t(z_{t-1}\rightarrow z_t)\right].
$$

Dynamic programming recursion:

$$
\delta_1(j)=\ln\pi_j+\ell'_1(j)
$$

$$
\delta_t(j)=\ell'_t(j)+\max_i\left[\delta_{t-1}(i)+\ln A_t(i\rightarrow j)\right].
$$

Backpointers store the argmax to reconstruct $\hat{z}_{1:T}$ in $O(TS^2)$, with $S=3$.

**(B) Posterior + hysteresis (`posterior_hysteresis`) — used here**  
Compute posteriors $P(z_t=s|x)$, form a per-marker best state:

$$
b_t=\arg\max_s P(z_t=s|x),
$$

then apply RTIGER-style rigidity (below). This preserves linkage information (via posteriors) and adds robustness against noisy switches.

**(C) Emission argmax + hysteresis (`emission_hysteresis`)**  
Set $b_t=\arg\max_s \ell'_t(s)$ (ignores transitions), then apply rigidity. Fast but less principled.

---

#### 6.2.5 RTIGER-style rigidity via hysteresis (`--rigidity / -R`)
After obtaining a "best state" sequence $b_t$ (from posteriors or emissions), we apply a hysteresis rule:

- Maintain current output state $c_t$.
- A change to a new state $u \neq c$ is only accepted after **R consecutive markers** support $u$.
- Chromosome boundaries reset the hysteresis counters.

This suppresses isolated single-marker flips while allowing genuine state changes supported over a stretch of markers.

---

#### 6.2.6 Post-decoding minimum-run cleanup (`--min_run_hh`, `--min_run_rh`)
After decoding (including hysteresis), we apply a run-length filter:

- HH runs shorter than `min_run_hh` are reassigned to RH (HH→RH),
- RH runs shorter than `min_run_rh` are reassigned to RR (RH→RR),

applied within each chromosome. This removes extremely short segments that are unlikely to represent true introgression at the marker density used.

---

### 6.3 HMM Parameters

Below are the key parameters commonly tuned or reported.

| Parameter | Meaning | Notes |
|-----------|---------|------|
| `--prior_rr`, `--prior_rh`, `--prior_hh` | Stationary/initial state priors | Normalized to π; also used to distribute off-diagonal transition mass |
| `--rho` | Sticky transition multiplier | Scales switching $s_t=\rho\theta_t$; higher → more switching |
| `--min_morgan` | Lower bound on Morgan distance | Prevents $\theta_t$ from becoming exactly zero when ΔcM≈0 |
| `--eta_hh_from_rh` | HH rescue from RH emissions | GL-mixture: HH borrows mass from RH; stabilizes HH in low-information regions |
| `--eta_rr_from_rh` | RR borrow from RH emissions | Usually 0; included for completeness |
| `--rh_penalty` | Penalize RH emission | Downweights RH likelihood by a constant |
| `--decode` | Decoding strategy | `posterior_hysteresis` recommended for robustness |
| `--rigidity` / `-R` | Hysteresis rigidity | Require R consecutive markers supporting a switch |
| `--min_run_hh`, `--min_run_rh` | Post-decoding cleanup | Removes short HH/RH runs by collapsing them toward less “donor-like” states |

**Example optimized values:**

| Parameter | Description | Example optimized value |
|-----------|-------------|--------------------------|
| `--rigidity` / `-R` | Consecutive markers needed to switch | **41** |
| `--rho` | Transition stickiness multiplier | **0.023** |
| `--eta_hh_from_rh` | HH rescue from RH emissions | **0.012** |
| `--prior_rr` | Prior probability of RR | 0.85 |
| `--prior_rh` | Prior probability of RH | 0.05 |
| `--prior_hh` | Prior probability of HH | 0.10 |
| `--decode` | Decoding method | `posterior_hysteresis` |
| `--min_run_hh` | Minimum HH run length (markers) | 3 |
| `--min_run_rh` | Minimum RH run length (markers) | 5 |

---

### 6.4 Parameter Optimization (Stability-based; no ground truth)

Since we do not have labeled local ancestry truth at each marker, we optimize parameters by maximizing **stability under marker thinning**:

1. Run the full pipeline on the full marker set → **baseline**
2. Create $K$ thinned replicates (randomly drop ~10% markers)
3. Run the pipeline on each replicate
4. Compute stability score = similarity between baseline and replicate outputs

**Stability metrics (objective components)**
- **Donor-block Jaccard similarity** after post-analysis (donor = AB+BB):

  $$
  J = \frac{bp(\text{baseline} \cap \text{replicate})}{bp(\text{baseline} \cup \text{replicate})}
  $$
  
- **State concordance** at shared markers:
  fraction of shared (CHROM,POS) with identical state calls
- **Fragmentation penalty**:
  penalize excessive donor breakpoints or overly fragmented donor tracts per Mb

We search parameter space using Bayesian optimization (Optuna), typically tuning:
- `R` (rigidity, integer),
- `rho` (transition scaling),
- `eta_hh_from_rh` (emission rescue strength),
while keeping priors and post-run filters fixed unless there is strong evidence they require tuning.

**Run optimization**
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

---

### 6.4 Run HMM per sample

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

### 6.5 Post-analysis: Block merging + visualization

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

### 6.6 Chromosome painting color scheme

- **Blue (AA)** = B73/B73 (Reference)
- **Purple (AB)** = Heterozygous
- **Red (BB)** = Teo/Teo (Donor)

---

### 6.7 Complete HMM workflow

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

### 6.8 HMM File Reference

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

## Section 7 — QTL / Association Mapping

(To be added)

---
