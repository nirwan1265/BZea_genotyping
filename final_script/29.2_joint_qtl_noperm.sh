#!/bin/bash
#SBATCH --job-name=qtl2_noperm
#SBATCH --output=logs/29.2_joint_qtl_noperm.%A_%a.out
#SBATCH --error=logs/29.2_joint_qtl_noperm.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --mem=24G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# 29.2 — Joint QTL model and per-taxon effects (R/qtl2) WITHOUT permutations, on the SAME QTL as the k-mer
#   GridLMM scan: the top windows and support intervals in peaks_<TRAIT>_ALL.tsv (step 27, laptop run).
#   R/qtl2 rescans the genome (no threshold) and, at each k-mer-scan peak, reports its own LOD / -log10 p next to
#   GridLMM's, the joint (conditional) effect, the per-taxon effects and the taxa-differ test.
#   Copy the laptop peaks first (one folder, all traits):
#     scp ~/Documents/Research/data/BZea/introgression_spread/assoc_gridlmm_spats/genes/peaks_*_ALL.tsv \
#         <server>:${D21}/assoc_gridlmm_spats/genes/
#   Submit:  N=$(wc -w <<< "$TRAITS"); sbatch --array=1-$N 29.2_joint_qtl_noperm.sh
#   Output:  ${D29}/noperm/qtl2_<TRAIT>_{qtl,by_taxon,genes,scan}.tsv, qtl2_<TRAIT>_manhattan.png  (~5 min per trait)
#   Later, with permutations: 29_joint_qtl.sh
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"
need Rscript
T=$(echo $TRAITS | tr ' ' '\n' | sed -n "${SLURM_ARRAY_TASK_ID:-1}p")
[[ -n "$T" ]] || { echo "no trait for task ${SLURM_ARRAY_TASK_ID:-1}"; exit 0; }
PEAKS_DIR=${PEAKS_DIR:-${D21}/assoc_gridlmm_spats/genes}
PEAKS=${PEAKS_DIR}/peaks_${T}_ALL.tsv
[[ -s "$PEAKS" ]] || { echo "$T: no $PEAKS (no significant k-mer peaks, or not copied yet), skipped"; exit 0; }
OUT=${D29}/noperm; mkdir -p "$OUT"
GFF1=$(ls $GFF_GLOB 2>/dev/null | head -1 || true)
if (( $(ls $GFF_GLOB 2>/dev/null | wc -l) > 1 )); then        # several per-chromosome GFFs: join the genes once
  GFF1=${D29}/b73_v5_genes.gff3
  if [[ ! -s "$GFF1" ]]; then cat $GFF_GLOB | awk -F'\t' '$3=="gene"' > "$GFF1.$$" && mv "$GFF1.$$" "$GFF1"; fi
fi
echo "trait $T | inputs $D25 | peaks $PEAKS | out $OUT"
NPERM=0 PEAKS="$PEAKS" GFF="$GFF1" CORES=${SLURM_CPUS_PER_TASK:-1} Rscript "${FINAL_DIR}/29_joint_qtl.R" "$D25" "$T" "$OUT"
