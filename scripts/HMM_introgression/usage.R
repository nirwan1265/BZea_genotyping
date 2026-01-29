# ============================================================================
#                    HMM INTROGRESSION ANALYSIS - USAGE GUIDE
# ============================================================================
#
# This file documents the complete workflow for introgression calling using
# rigidity-based smoothing (marker-count approach, inspired by RTIGER).
#
# Author: Nirwan Tandukar
# Date: 2025-01-28
#
# ============================================================================

# ============================================================================
# OVERVIEW
# ============================================================================
#
# The pipeline has these main steps:
#
#   1. HMM calling (scripts 1-4):
#      - Extracts genotype likelihoods from BCF
#      - Masks unreliable B73 positions
#      - Runs HMM (Viterbi + forward-backward) to call states per SNP
#      - Output: statepath.tsv.gz files with per-SNP states + posteriors
#
#   2. Introgression tract calling (this workflow):
#      - Applies rigidity-based smoothing to remove noisy short segments
#      - Merges adjacent donor segments into continuous tracts
#      - Output: tract files with introgression regions
#
# Key difference from original cM-based bridging:
#   - OLD: Bridge gaps if <= 0.10 cM AND <= 25 Mb
#   - NEW: Smooth segments with < min_markers SNPs (RTIGER-style rigidity)
#
# ============================================================================


# ============================================================================
# STEP 1: ANALYZE SINGLE SAMPLE
# ============================================================================
#
# Use analyze_sample.R to process one sample:
#
#   Rscript analyze_sample.R <input.statepath.tsv.gz> [min_markers]
#
# Examples:
#   Rscript analyze_sample.R PN12_SID1120.statepath.tsv.gz 5
#   Rscript analyze_sample.R /path/to/sample.statepath.tsv.gz 10
#
# Arguments:
#   - input.statepath.tsv.gz: Output from HMM (script 4)
#   - min_markers: Rigidity value (default: 5)
#     - Higher = more smoothing (fewer small segments)
#     - Lower = less smoothing (keep more detail)
#     - Recommended: 3-10
#
# Output:
#   - Console: Summary statistics
#   - File: <sample>.tracts.rigidity<N>.tsv
#
# ============================================================================


# ============================================================================
# STEP 2: OPTIMIZE RIGIDITY VALUE (OPTIONAL)
# ============================================================================
#
# Use optimize_rigidity.R to find the best min_markers value:
#
#   Rscript optimize_rigidity.R
#
# This script:
#   1. Tests multiple rigidity values (1, 2, 3, 5, 7, 10, 15, 20, 30, 50)
#   2. Measures for each value:
#      - Posterior improvement (higher = removing low-confidence noise)
#      - Removal quality ratio (good removals / bad removals)
#      - Donor signal loss (how much real introgression we lose)
#   3. Computes composite score to find optimal value
#   4. Outputs results table + PDF plot
#
# IMPORTANT: Edit the script first to set:
#   - state_dir: Path to your statepath files
#   - n_samples_to_test: Number of samples to use (20 for quick test)
#
# Output:
#   - rigidity_optimization_results.tsv
#   - rigidity_optimization_plot.pdf
#
# ============================================================================


# ============================================================================
# STEP 3: BATCH PROCESS ALL SAMPLES
# ============================================================================
#
# Use introgression_stats_rigidity.R for batch processing:
#
#   Rscript introgression_stats_rigidity.R
#
# This script processes all samples in state_dir and:
#   1. Applies rigidity-based smoothing
#   2. Refines boundaries using posterior probabilities
#   3. Extracts donor tracts
#   4. Computes per-sample and genome-wide statistics
#
# IMPORTANT: Edit the script first to set:
#   - state_dir: Path to your statepath files
#   - out_dir: Where to save results
#   - min_markers: Rigidity value (use optimal from Step 2)
#   - smoothing_method: "probability" (recommended) or "simple"
#
# Output:
#   - Per-sample tract files: <sample>.alt_het.tracts.rigidity.tsv.gz
#   - Summary: all_samples.ref_alt_het_perc.rigidity.tsv
#   - Combined tracts: all_samples.alt_het.tracts.rigidity.tsv.gz
#
# ============================================================================


# ============================================================================
# EXAMPLE: COMPLETE WORKFLOW
# ============================================================================

# --- Option A: Quick single-sample analysis ---

# In terminal:
# cd /path/to/HMM_introgression
# Rscript analyze_sample.R /path/to/sample.statepath.tsv.gz 5


# --- Option B: Full batch analysis ---

# Step 1: Optimize rigidity (optional but recommended)
# Edit optimize_rigidity.R to set paths, then:
# Rscript optimize_rigidity.R

# Step 2: Run batch processing with optimal rigidity
# Edit introgression_stats_rigidity.R:
#   min_markers <- 5  # Use value from optimization
# Then:
# Rscript introgression_stats_rigidity.R


# ============================================================================
# UNDERSTANDING THE OUTPUT
# ============================================================================
#
# Tract file columns:
#   - Chr: Chromosome
#   - Start: Start position (bp)
#   - End: End position (bp)
#   - State: "HH" (homozygous teosinte) or "RH" (heterozygous)
#   - n_SNPs: Number of SNPs in tract
#   - len_bp: Length in base pairs
#   - len_Mb: Length in megabases
#   - mean_posterior: Average posterior probability for called state
#
# Summary file columns:
#   - Sample: Sample name
#   - Ref_perc: % genome as RR (B73 reference)
#   - Alt_perc: % genome as HH (teosinte)
#   - Het_perc: % genome as RH (heterozygous)
#   - *_bp: Absolute values in base pairs
#
# ============================================================================


# ============================================================================
# RTIGER vs YOUR HMM: KEY DIFFERENCES
# ============================================================================
#
# 1. MULTI-SAMPLE USAGE:
#
#    RTIGER:
#    - Uses ALL samples together for EM training
#    - Learns emission parameters (Beta-Binomial) from population
#    - Then applies Viterbi PER SAMPLE for state calling
#
#    Your HMM:
#    - Uses fixed emission parameters (genotype likelihoods)
#    - Runs Viterbi PER SAMPLE independently
#    - Multi-sample benefit comes from B73 masking (script 2)
#
#    Bottom line: Both do per-sample Viterbi. RTIGER's multi-sample
#    advantage is in parameter learning, not in final calling.
#
# 2. RIGIDITY IMPLEMENTATION:
#
#    RTIGER:
#    - Rigidity is BUILT INTO the Viterbi algorithm
#    - Requires N consecutive markers to switch states
#    - Done DURING HMM inference
#
#    Your approach (this workflow):
#    - Rigidity applied POST-HOC after Viterbi
#    - Short segments smoothed using posterior probabilities
#    - Done AFTER HMM inference
#
#    Both achieve similar results, but RTIGER's is more principled.
#    Your approach is simpler and uses the posteriors you already have.
#
# 3. DATA INPUT:
#
#    RTIGER: Allele counts -> Beta-Binomial emissions
#    Your HMM: Genotype likelihoods -> Direct emissions
#
# ============================================================================


# ============================================================================
# RECOMMENDED PARAMETERS
# ============================================================================
#
# For BC2S3 populations:
#
#   min_markers = 5
#     - Removes noise while preserving signal
#     - Good balance for most datasets
#
#   smoothing_method = "probability"
#     - Uses posterior probabilities for decisions
#     - More accurate than simple neighbor-based smoothing
#
#   donor_set = c("HH", "RH")
#     - Includes both homozygous and heterozygous teosinte
#     - Use "HH" only if you want strict homozygous introgression
#
# ============================================================================


# ============================================================================
# QUICK R CODE TO RUN ANALYSIS
# ============================================================================

# --- Load libraries ---
library(data.table)

# --- Set paths ---
input_file <- "/path/to/sample.statepath.tsv.gz"
min_markers <- 5

# --- Read data ---
dt <- fread(cmd = paste("gunzip -c", shQuote(input_file)))

# --- Check original states ---
table(dt$STATE)

# --- Source the analyze script for functions (or run directly) ---
# source("analyze_sample.R")  # If you want to use as library

# --- Or just run from command line ---
# system("Rscript analyze_sample.R /path/to/sample.statepath.tsv.gz 5")


# ============================================================================
# TROUBLESHOOTING
# ============================================================================
#
# 1. "gunzip: command not found"
#    - On Mac: Should work by default
#    - On cluster: module load gzip
#
# 2. "No statepath files found"
#    - Check state_dir path in script
#    - Make sure files end with .statepath.tsv.gz
#
# 3. "Memory error"
#    - Process fewer samples at once
#    - Increase memory allocation: srun --mem=32G ...
#
# 4. "Posterior columns missing"
#    - Your statepath files need P_RR, P_RH, P_HH columns
#    - These come from HMM forward-backward (script 4)
#    - If missing, use smoothing_method = "simple"
#
# ============================================================================


# ============================================================================
# FILE SUMMARY
# ============================================================================
#
# Scripts in this folder:
#
#   analyze_sample.R
#     - Single sample analysis
#     - Usage: Rscript analyze_sample.R <input> [min_markers]
#
#   optimize_rigidity.R
#     - Find optimal min_markers value
#     - Usage: Rscript optimize_rigidity.R
#
#   introgression_stats_rigidity.R
#     - Batch process all samples
#     - Usage: Rscript introgression_stats_rigidity.R
#
#   introgression_stats.R (ORIGINAL)
#     - Old cM-based bridging approach
#     - Kept for comparison
#
#   usage.R (THIS FILE)
#     - Documentation and examples
#
# ============================================================================
