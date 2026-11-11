#!/bin/bash
#SBATCH --job-name=rename
#SBATCH --array=1-10
#SBATCH --output=logs/06_rename.%A_%a.out
#SBATCH --error=logs/06_rename.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mem=16G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# 06 — Replace BAM paths in the BCF header with PN*_SID* sample IDs.
#      (replaces 0.5_rename_headers_files.sh; builds the name list itself and
#       checks it against the BAM list instead of relying on a hand-made file)
#
#   Output: ${D06}/BZea_<chr>.full.renamed.bcf (+ .csi)
#   Submit: sbatch 06_rename_samples.sh
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"

IDX=${SLURM_ARRAY_TASK_ID:-${1:?usage: bash $0 <chrom index 1-10>}}
CHR=${CHROMS[$(( IDX - 1 ))]}
mkdir -p "$D06"
IN=${D05}/BZea_${CHR}.full.bcf
OUT=${D06}/BZea_${CHR}.full.renamed.bcf
[[ -s $IN ]] || die "missing ${IN} (run 05)"

NEW=${D06}/BZea_${CHR}.newnames.txt
OLD=${D06}/BZea_${CHR}.oldnames.txt
bcftools query -l "$IN" > "$OLD"
sample_names_from_bamlist > "$NEW"

# ANGSD names samples by BAM path, in BAM-list order. Check the order before renaming.
sed -E 's#^.*/##; s/\.bam$//' "$OLD" | sed -E 's/^.*(PN[0-9]+_SID[0-9]+).*$/\1/' | cmp -s - "$NEW" \
  || die "BCF sample order differs from ${BAMLIST}"
[[ $(sort "$NEW" | uniq -d | wc -l) -eq 0 ]] || die "duplicated sample IDs in ${NEW}"

bcftools reheader -s "$NEW" -o "$OUT" "$IN"
bcftools index -f "$OUT"
echo "First samples: $(bcftools query -l "$OUT" | head -n 3 | tr '\n' ' ')"
echo "WROTE: ${OUT}"
