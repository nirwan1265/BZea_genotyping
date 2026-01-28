#!/usr/bin/env bash
#BSUB -J "maskGL[1-1600]%200"
#BSUB -n 1
#BSUB -W 04:00
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=2GB]"
#BSUB -o /rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf/HMM_outputs/logs/maskGL.%J.%I.out
#BSUB -e /rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf/HMM_outputs/logs/maskGL.%J.%I.err


#GLDIR="/rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf/HMM_inputs_GL"
#LIST="${GLDIR}/gl_files.list"

#ls -1 "${GLDIR}"/*.chr1-10.GL.tsv.gz | sort > "$LIST"
#wc -l "$LIST"
#head "$LIST"




set -euo pipefail

GLDIR="/rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf/HMM_inputs_GL"
INLIST="${GLDIR}/gl_files.list"

MASKPOS="/rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf/b73_mask/allchr.b73_bad.pos"  # chr<TAB>pos

OUTDIR="/rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf/HMM_inputs_GL_maskedB73"
LOGDIR="/rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf/HMM_outputs/logs"
mkdir -p "$OUTDIR" "$LOGDIR"

INFILE="$(sed -n "${LSB_JOBINDEX}p" "$INLIST")"
if [[ -z "${INFILE}" ]]; then
  echo "ERROR: No input file for JOBINDEX=${LSB_JOBINDEX} in $INLIST" >&2
  echo "List length = $(wc -l < "$INLIST")" >&2
  exit 2
fi

base="$(basename "$INFILE")"
sample="${base%%.chr1-10.GL.tsv.gz}"

OUT="${OUTDIR}/${sample}.chr1-10.GL.maskedB73.tsv.gz"

echo "JOBINDEX=${LSB_JOBINDEX}"
echo "INFILE=${INFILE}"
echo "MASKPOS=${MASKPOS}"
echo "OUT=${OUT}"

# Remove masked positions (exact match on CHROM+POS)
zcat "$INFILE" \
| awk -v OFS="\t" '
    NR==FNR { key=$1 "\t" $2; mask[key]=1; next }
    {
      key=$1 "\t" $2
      if(!(key in mask)) print
    }
  ' "$MASKPOS" - \
| bgzip -c > "$OUT"

echo "DONE: $OUT"
