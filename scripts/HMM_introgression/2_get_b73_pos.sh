#!/usr/bin/env bash
#BSUB -J "b73MaskStrict[1-10]%10"
#BSUB -n 1
#BSUB -W 01:00
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=4GB]"
#BSUB -o /rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf/b73_mask/logs/b73MaskStrict.%J.%I.out
#BSUB -e /rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf/b73_mask/logs/b73MaskStrict.%J.%I.err



#cat > b73.list <<'EOF'
#PN3_SID203
#PN3_SID236
#PN3_SID255
#PN4_SID326
#PN5_SID429
#PN5_SID442
#PN5_SID468
#PN6_SID576
#PN8_SID693
#PN8_SID751
#PN9_SID838
#PN10_SID893
#EOF

set -euo pipefail

BCFDIR="/rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf"
OUTDIR="/rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf/b73_mask_strict"
B73LIST="/share/maize/ntanduk/angsd_genotyping/HMM_introgression/b73.list"

mkdir -p "$OUTDIR"

CHR="chr${LSB_JOBINDEX}"
BCF="${BCFDIR}/${CHR}.teosnp.clean.rehead.bcf"

# --------- TUNE THESE ----------
MIN_CALLED_B73=12              # require >= this many B73 samples called
MIN_AN_B73=$((2*MIN_CALLED_B73))
# --------------------------------

echo "CHR=$CHR"
echo "BCF=$BCF"
echo "B73LIST=$B73LIST"
echo "MIN_AN_B73=$MIN_AN_B73"

# BAD if:
#   (AN < MIN_AN_B73)  OR  (AC > 0)
# i.e. B73 not confidently RR & well-called
bcftools view -Ou -S "$B73LIST" "$BCF" \
  | bcftools +fill-tags -Ou -- -t AC,AN \
  | bcftools query -f'%CHROM\t%POS\t%AC\t%AN\n' \
  | awk -v min_an="$MIN_AN_B73" '
      ($4 < min_an) || ($3 > 0) { print $1 "\t" $2 }
    ' \
  | sort -k1,1 -k2,2n \
  > "${OUTDIR}/${CHR}.b73_notclean.pos"

wc -l "${OUTDIR}/${CHR}.b73_notclean.pos"
echo "WROTE: ${OUTDIR}/${CHR}.b73_notclean.pos"
