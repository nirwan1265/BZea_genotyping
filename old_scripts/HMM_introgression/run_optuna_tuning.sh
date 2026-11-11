#!/bin/bash
# ============================================================================
# Run HMM parameter tuning with Optuna
# ============================================================================
#
# This script sets up the environment and runs the Optuna-based parameter
# tuning for HMM introgression calling.
#
# Prerequisites:
#   1. conda environment "hmm" with optuna installed
#   2. Genetic map files split into chr1.map, chr2.map, ..., chr10.map
#   3. Calibration samples file (calib_samples.txt) with paths to GL files
#
# Usage:
#   ./run_optuna_tuning.sh [n_trials]
#
# ============================================================================

set -e

# Directory containing this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Activate conda environment
echo "Activating conda environment 'hmm'..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate hmm

# ============================================================================
# CONFIGURATION - Edit these paths as needed
# ============================================================================

# Path to HMM script
export HMM_PY="${SCRIPT_DIR}/hmm_viterbi_bc2s3.py"

# Directory containing genetic maps (chr1.map, chr2.map, ...)
# If maps don't exist, split them from NAM_genetic_map.txt first
export MAPDIR="${SCRIPT_DIR}/genetic_maps"

# Path to R post-analysis script
export POST_R="${SCRIPT_DIR}/batch_introgression_analysis.R"

# Output directory for tuning results
export OUTROOT="${SCRIPT_DIR}/optuna_tuning_runs"

# Post-analysis parameters (frozen during tuning)
export MAX_GAP_KB=10000    # Max gap to bridge between same-state blocks (kb)
export MIN_BLOCK_KB=5000   # Minimum block size to keep (kb)

# Stability experiment parameters
export DROP_FRAC=0.10      # Fraction of markers to drop in thinning (10 percent)
export K_REPS=3            # Number of thinned replicates per sample

# Calibration samples file
CALIB_FILE="${SCRIPT_DIR}/calib_samples.txt"

# Number of optimization trials
N_TRIALS=${1:-50}

# ============================================================================
# PRE-FLIGHT CHECKS
# ============================================================================

echo ""
echo "=== Pre-flight checks ==="

# Check HMM script
if [[ ! -f "$HMM_PY" ]]; then
    echo "ERROR: HMM script not found: $HMM_PY"
    exit 1
fi
echo "HMM script: $HMM_PY [OK]"

# Check R script
if [[ ! -f "$POST_R" ]]; then
    echo "ERROR: R script not found: $POST_R"
    exit 1
fi
echo "R script: $POST_R [OK]"

# Check/create genetic maps
if [[ ! -d "$MAPDIR" ]]; then
    echo "Genetic map directory not found. Creating from NAM_genetic_map.txt..."
    mkdir -p "$MAPDIR"

    if [[ -f "${SCRIPT_DIR}/NAM_genetic_map.txt" ]]; then
        python "${SCRIPT_DIR}/split_genetic_map.py" "${SCRIPT_DIR}/NAM_genetic_map.txt" "$MAPDIR"
    else
        echo "ERROR: NAM_genetic_map.txt not found in ${SCRIPT_DIR}"
        exit 1
    fi
fi

# Check that all chr maps exist
for i in {1..10}; do
    if [[ ! -f "${MAPDIR}/chr${i}.map" ]]; then
        echo "ERROR: Missing ${MAPDIR}/chr${i}.map"
        echo "Run: python ${SCRIPT_DIR}/split_genetic_map.py NAM_genetic_map.txt $MAPDIR"
        exit 1
    fi
done
echo "Genetic maps: $MAPDIR [OK]"

# Check calibration file
if [[ ! -f "$CALIB_FILE" ]]; then
    echo "ERROR: Calibration file not found: $CALIB_FILE"
    exit 1
fi
echo "Calibration file: $CALIB_FILE [OK]"

# Check sample files exist
echo "Checking sample files..."
n_samples=0
n_missing=0
while IFS= read -r line; do
    # Skip empty lines and comments
    [[ -z "$line" || "$line" == \#* ]] && continue
    if [[ ! -f "$line" ]]; then
        echo "  WARNING: Sample not found: $line"
        ((n_missing++))
    else
        ((n_samples++))
    fi
done < "$CALIB_FILE"
echo "Found $n_samples samples, $n_missing missing"

if [[ $n_samples -eq 0 ]]; then
    echo "ERROR: No valid samples found!"
    exit 1
fi

# ============================================================================
# RUN OPTIMIZATION
# ============================================================================

echo ""
echo "=== Starting Optuna optimization ==="
echo "Trials: $N_TRIALS"
echo "Output: $OUTROOT"
echo ""

mkdir -p "$OUTROOT"

python "${SCRIPT_DIR}/tune_hmm_optuna.py" "$CALIB_FILE" "$N_TRIALS"

echo ""
echo "=== Done ==="
echo "Results saved to: $OUTROOT"
echo "Best parameters: $OUTROOT/best_params.txt"
