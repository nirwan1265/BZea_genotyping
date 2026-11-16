#!/bin/bash
#SBATCH --job-name=joint_qtl
#SBATCH --output=logs/29_joint_qtl.%A_%a.out
#SBATCH --error=logs/29_joint_qtl.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=04:00:00
#SBATCH --mem=12G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# 29 — Joint multi-QTL model with per-taxon QTL effects (29_joint_qtl.R). One array task per trait in $TRAITS.
#   Needs the step-25 inputs (${D25}) and R with data.table only (no GridLMM).
#   Per trait: 1,000 permutations (~1 s each) for the threshold, forward/backward selection, final model.
#   Submit:  N=$(wc -w <<< "$TRAITS"); sbatch --array=1-$N 29_joint_qtl.sh
#   Output:  ${D29}/joint_qtl_<TRAIT>{,_by_taxon,_genes,_steps,_perm}.tsv
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"
need Rscript
T=$(echo $TRAITS | tr ' ' '\n' | sed -n "${SLURM_ARRAY_TASK_ID:-1}p")
[[ -n "$T" ]] || { echo "no trait for task ${SLURM_ARRAY_TASK_ID:-1}"; exit 0; }
GFF1=$(ls $GFF_GLOB 2>/dev/null | head -1 || true)
if (( $(ls $GFF_GLOB 2>/dev/null | wc -l) > 1 )); then        # several per-chromosome GFFs: join them once
  GFF1=${D29}/b73_v5_genes.gff3; mkdir -p "$D29"
  if [[ ! -s "$GFF1" ]]; then cat $GFF_GLOB | awk -F'\t' '$3=="gene"' > "$GFF1.$$" && mv "$GFF1.$$" "$GFF1"; fi
fi
echo "trait $T | inputs $D25 | out $D29 | GFF ${GFF1:-none}"
GFF="$GFF1" NPERM=${NPERM:-1000} Rscript "${FINAL_DIR}/29_joint_qtl.R" "$D25" "$T" "$D29"
