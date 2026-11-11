#!/bin/bash
#BSUB -J "rename[1-10]%10"
#BSUB -n 1
#BSUB -W 12:00
#BSUB -q sara
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=16GB]"
#BSUB -o logs/06_rename.%J.%I.out
#BSUB -e logs/06_rename.%J.%I.err
# =============================================================================
# 06 — Replace BAM paths in the BCF header with PN*_SID* sample IDs.
#      (replaces 0.5_rename_headers_files.sh; builds the name list itself and
#       checks it against the BAM list instead of relying on a hand-made file)
#
#   Output: ${D06}/BZea_<chr>.full.renamed.bcf (+ .csi)
#   Submit: bsub < 06_rename_samples.sh
# =============================================================================
set -euo pipefail
source "${LS_SUBCWD:-$PWD}/config.sh"

IDX=${LSB_JOBINDEX:-${1:?usage: bash $0 <chrom index 1-10>}}
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
