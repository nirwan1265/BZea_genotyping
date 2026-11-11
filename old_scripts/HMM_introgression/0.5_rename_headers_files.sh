#!/bin/bash
#BSUB -J "rehead_full[1-10]%10"
#BSUB -n 1
#BSUB -W 12:00
#BSUB -q sara
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=16GB]"
#BSUB -o logs/rehead_full.%J.%I.out
#BSUB -e logs/rehead_full.%J.%I.err




#DO THIS ONCE

#cd /rsstu/users/r/rrellan/BZea/angsd_genotyping/angsd_full_snps

#BCF=per_chr/BZea_chr1.full.bcf

# old names exactly as stored in the BCF (full BAM paths in your case)
#bcftools query -l "$BCF" > rename.old.list

# new names: extract PN*_SID### from anywhere in the string
# (works even if basename has extra suffixes, as long as PN..._SID... exists)
#awk '{
#  if (match($0, /(PN[0-9]+_SID[0-9]+)/, m)) print m[1];
#  else {
#    # fallback: basename without .bam
#    n=$0; sub(/^.*\//,"",n); sub(/\.bam$/,"",n); print n
#  }
#}' rename.old.list > rename.new.list

# sanity check
#paste rename.old.list rename.new.list | head


set -euo pipefail

ROOT=/rsstu/users/r/rrellan/BZea/angsd_genotyping/angsd_full_snps
INDIR=${ROOT}/per_chr
OUTDIR=${ROOT}/per_chr_renamed

OLDLIST=${INDIR}/rename.old.list
NEWLIST=${INDIR}/rename.new.list

mkdir -p "${OUTDIR}" "${ROOT}/logs"

CHR="chr${LSB_JOBINDEX}"

INBCF="${INDIR}/BZea_${CHR}.full.bcf"
OUTBCF="${OUTDIR}/BZea_${CHR}.full.renamed.bcf"

echo "CHR=$CHR"
echo "INBCF=$INBCF"
echo "OUTBCF=$OUTBCF"

[[ -s "$INBCF" ]] || { echo "Missing input: $INBCF" >&2; exit 2; }
[[ -s "$NEWLIST" ]] || { echo "Missing rename list: $NEWLIST" >&2; exit 3; }

# reheader (order-based mapping: line i in NEWLIST replaces sample i in header)
bcftools reheader -s "$NEWLIST" -o "$OUTBCF" "$INBCF"

# index
bcftools index -f "$OUTBCF"

# sanity
echo "First 5 sample names:"
bcftools query -l "$OUTBCF" | head -n 5

echo "DONE: $OUTBCF"
