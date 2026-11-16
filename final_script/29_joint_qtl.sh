#!/bin/bash
#SBATCH --job-name=qtl2_joint
#SBATCH --output=logs/29_joint_qtl.%A_%a.out
#SBATCH --error=logs/29_joint_qtl.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --mem=64G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# 29 — R/qtl2: scan, permutation thresholds, QTL, joint model and per-taxon effects (29_joint_qtl.R).
#   One array task per trait in $TRAITS. Needs the step-25 inputs (${D25}) and R with qtl2 + data.table:
#     Rscript -e 'install.packages(c("qtl2","data.table"), repos = "https://cloud.r-project.org")'
#   Time per trait on 16 cores: ~5 min for the scans, ~1-2 h for 2 x 1,000 permutations (~22 s each on 1 core).
#   Submit:  N=$(wc -w <<< "$TRAITS"); sbatch --array=1-$N 29_joint_qtl.sh
#   Output:  ${D29}/qtl2_<TRAIT>_{qtl,by_taxon,genes,scan,perm}.tsv, qtl2_<TRAIT>_manhattan.png
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"
need Rscript
T=$(echo $TRAITS | tr ' ' '\n' | sed -n "${SLURM_ARRAY_TASK_ID:-1}p")
[[ -n "$T" ]] || { echo "no trait for task ${SLURM_ARRAY_TASK_ID:-1}"; exit 0; }
mkdir -p "$D29"
GFF1=$(ls $GFF_GLOB 2>/dev/null | head -1 || true)
if (( $(ls $GFF_GLOB 2>/dev/null | wc -l) > 1 )); then        # several per-chromosome GFFs: join the genes once
  GFF1=${D29}/b73_v5_genes.gff3
  if [[ ! -s "$GFF1" ]]; then cat $GFF_GLOB | awk -F'\t' '$3=="gene"' > "$GFF1.$$" && mv "$GFF1.$$" "$GFF1"; fi
fi
echo "trait $T | inputs $D25 | out $D29 | GFF ${GFF1:-none} | cores ${SLURM_CPUS_PER_TASK:-1}"
GFF="$GFF1" NPERM=${NPERM:-1000} CORES=${SLURM_CPUS_PER_TASK:-1} Rscript "${FINAL_DIR}/29_joint_qtl.R" "$D25" "$T" "$D29"
