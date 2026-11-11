# HMM Introgression Calling Pipeline

A Hidden Markov Model (HMM) pipeline for calling local ancestry / introgression states in BC2S3 maize × teosinte samples using genotype likelihoods.

---

## 1. HMM Model Overview

**Goal:** Call local ancestry states (RR, RH, HH) across the genome in BC2S3 maize × teosinte samples.

### States

| State | Meaning | Genotype |
|-------|---------|----------|
| RR | Homozygous Reference (B73/B73) | AA |
| RH | Heterozygous (B73/Teosinte) | AB |
| HH | Homozygous Teosinte (Teo/Teo) | BB |

### Input Format

Genotype Likelihoods (GL) per SNP in TSV format (gzipped):

```
CHROM  POS     REF  ALT  DP  glRR    glRH    glHH
chr1   12345   A    G    15  -0.5    -2.3    -15.1
```

- `glRR`, `glRH`, `glHH` are log10 genotype likelihoods (relative or absolute)

---

## 2. HMM Components

### A. Emissions (Observation Probabilities)

- Convert GL from log10 → natural log scale
- **Adjustments:**
  - `eta_hh_from_rh`: Allows HH to "borrow" likelihood from RH (helps when low depth causes RH overcalling)
  - `rh_penalty`: Subtract penalty from RH log-likelihood to reduce spurious het calls

### B. Transitions

- Use **genetic distance** (cM) between adjacent markers
- Convert cM → Morgan → Haldane's θ (recombination probability)
- **Sticky transition matrix:**
  ```
  switch_prob = rho × theta
  P(stay in state) = 1 - switch_prob
  P(switch to other) = switch_prob × prior[other] / (1 - prior[current])
  ```
- Lower `rho` = stickier (fewer state switches)

### C. Decoding Methods

| Method | Description |
|--------|-------------|
| `viterbi` | Standard MAP path, no rigidity |
| `posterior_hysteresis` | Forward-backward posteriors → argmax → hysteresis rigidity **(recommended)** |
| `emission_hysteresis` | Argmax emission only → hysteresis (fast, ignores transitions) |

### D. Rigidity (R) - RTIGER-style

- Requires **R consecutive markers** supporting a new state before switching
- Prevents noisy single-SNP state flips
- Higher R = smoother calls, fewer false switches

---

## 3. HMM Parameters

| Parameter | Description | Optimized Value |
|-----------|-------------|-----------------|
| `--rigidity` / `-R` | Consecutive markers needed to switch | **41** |
| `--rho` | Transition stickiness multiplier | **0.023** |
| `--eta_hh_from_rh` | HH rescue from RH emissions | **0.012** |
| `--prior_rr` | Prior probability of RR state | 0.85 |
| `--prior_rh` | Prior probability of RH state | 0.05 |
| `--prior_hh` | Prior probability of HH state | 0.10 |
| `--rh_penalty` | Penalty on RH log-likelihood | 0 |
| `--min_run_hh` | Min consecutive HH markers to keep | 3 |
| `--min_run_rh` | Min consecutive RH markers to keep | 5 |
| `--decode` | Decoding method | `posterior_hysteresis` |

---

## 4. Parameter Optimization (No Ground Truth)

Since we have **no labeled truth data**, we use **stability-based optimization**.

### Concept

```
If parameters are good → results should be STABLE under small perturbations
If parameters are bad → results will be NOISY and unstable
```

### Process

```
For each trial (R, rho, eta combination):
    1. Run HMM on FULL data → baseline results
    2. Create K thinned replicates (drop 10% markers randomly)
    3. Run HMM on each replicate → replicate results
    4. Compare baseline vs replicates:
       - Jaccard similarity of donor blocks (AB + BB)
       - State concordance at shared markers
       - Breakpoints per Mb (fragmentation penalty)
       - Het% guardrail (penalize if het collapses to ~0)
    5. Score = Jaccard + 0.5×Concordance - 0.02×Fragmentation - HetPenalty
```

### Optimization

- **Optuna** (Bayesian optimization) searches parameter space
- 50-100 trials typically sufficient
- Maximizes stability score

### Files Needed

```
tune_hmm_optuna.py              # Optimizer
hmm_viterbi_bc2s3.py            # HMM
batch_introgression_analysis.R  # Post-processing
genetic_maps/chr*.map           # Genetic maps
calib_samples.txt               # Sample list
```

### Run Command

```bash
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

### Output

```
optuna_tuning_runs/
├── best_params.txt      # Best parameters found
└── optuna_study.db      # SQLite database with all trials
```

---

## 5. Running the HMM

### Single Sample

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

### Outputs

1. **`sample.statepath.tsv.gz`** - Per-marker state calls + posteriors
   ```
   CHROM  POS    REF  ALT  DP  STATE  P_RR   P_RH   P_HH   THETA
   chr1   12345  A    G    15  RR     0.95   0.04   0.01   0.0001
   ```

2. **`sample.tracts.bed.gz`** - Merged segments
   ```
   chrom  start    end      state
   chr1   0        5000000  RR
   chr1   5000001  8000000  RH
   ```

---

## 6. Post-Analysis (R Script)

### Purpose

- Merge adjacent same-state blocks
- Bridge small gaps between same-state blocks
- Filter out tiny blocks
- Generate summary statistics and plots

### Parameters

| Parameter | Description | Value |
|-----------|-------------|-------|
| `max_gap_kb` | Max gap to bridge between same-state blocks | 10000 kb |
| `min_block_kb` | Minimum block size to keep | 5000 kb |

### Run Command

```bash
Rscript batch_introgression_analysis.R input_folder/ 10000 5000 output_folder/
```

### Outputs Per Sample

| File | Description |
|------|-------------|
| `sample.complete_blocks.bed` | Final introgression blocks |
| `sample.complete_blocks.png` | Chromosome painting visualization |
| `sample.block_summary.tsv` | Block counts by genotype |
| `sample.chr_summary.tsv` | Per-chromosome statistics |

### Aggregate Outputs

| File | Description |
|------|-------------|
| `introgression_summary_total.tsv` | Ref%, Het%, Teo%, Donor% per sample |
| `introgression_summary_per_chr_donor.tsv` | Donor% by chromosome |
| `introgression_summary_per_chr_het.tsv` | Het% by chromosome |
| `introgression_summary_per_chr_teo.tsv` | Teo% by chromosome |

---

## 7. Plotting

### Automatic Plot (from R script)

Each sample gets a chromosome painting:

```
sample.complete_blocks.png
```

**Color scheme:**
- **Blue (AA)** = B73/B73 (Reference)
- **Purple (AB)** = Heterozygous
- **Red (BB)** = Teo/Teo (Donor)

### Custom Plotting

Use the BED files for custom visualization:

```r
library(data.table)
blocks <- fread("sample.complete_blocks.bed",
                col.names = c("chrom", "start", "end", "genotype"))
```

---

## 8. Complete Workflow

```
┌─────────────────────────────────────────────────────────────┐
│  1. PREPARE INPUTS                                          │
│     - GL files (from ANGSD)                                 │
│     - Genetic maps (chr1.map ... chr10.map)                 │
│     - Sample list                                           │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│  2. PARAMETER OPTIMIZATION (once)                           │
│     python tune_hmm_optuna.py calib_samples.txt 100         │
│     → best_params.txt (R=41, rho=0.023, eta=0.012)          │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│  3. RUN HMM (per sample)                                    │
│     python hmm_viterbi_bc2s3.py --rigidity 41 --rho 0.023   │
│     → sample.statepath.tsv.gz                               │
│     → sample.tracts.bed.gz                                  │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│  4. POST-ANALYSIS (batch)                                   │
│     Rscript batch_introgression_analysis.R folder/ 10000 5000│
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

## 9. File Reference

| File | Purpose |
|------|---------|
| `hmm_viterbi_bc2s3.py` | Main HMM script |
| `batch_introgression_analysis.R` | Post-processing + plotting |
| `tune_hmm_optuna.py` | Parameter optimization |
| `split_genetic_map.py` | Split NAM map into chr*.map |
| `genetic_maps/chr*.map` | Per-chromosome genetic maps |
| `calib_samples.txt` | Sample list for optimization |
| `best_params.txt` | Optimized parameters |

---

## 10. Dependencies

### Python
- Python 3.8+
- optuna
- numpy (optional, for faster computation)

### R
- data.table

### Installation

```bash
# Create conda environment
conda create -n hmm python=3.10
conda activate hmm
pip install optuna

# R packages
Rscript -e "install.packages('data.table')"
```

---

## 11. Citation

If you use this pipeline, please cite:

- RTIGER method: Campos-Martin et al. (2023) for the rigidity/hysteresis concept
- Haldane mapping function for recombination probability calculation
- Optuna: Akiba et al. (2019) for Bayesian optimization

---

## 12. Troubleshooting

### Common Issues

1. **"Map too small" error**
   - Ensure genetic map files have at least 2 markers
   - Check file format: `bp_position\tcM_position`

2. **0% heterozygosity**
   - This is biologically valid for some BC2S3 samples
   - The optimizer applies a mild penalty but doesn't exclude these samples

3. **Missing rigidity parameter**
   - Add `--rigidity 41` to your HMM command
   - This is the most important tuned parameter

4. **R script fails with NA errors**
   - Ensure STATE column contains only RR, RH, HH
   - Check for empty statepath files

### Contact

For issues, open a GitHub issue or contact the maintainers.
