#!/usr/bin/env bash
#BSUB -J "stateStats[1-1600]%200"
#BSUB -n 1
#BSUB -W 01:00
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=2GB]"
#BSUB -o /rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf/HMM_outputs/logs/stateStats.%J.%I.out
#BSUB -e /rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf/HMM_outputs/logs/stateStats.%J.%I.err

set -euo pipefail

STATEDIR="/rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf/HMM_outputs/statepaths"
OUTDIR="/rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf/HMM_outputs/state_stats"
LOGDIR="/rsstu/users/r/rrellan/BZea/angsd_genotyping/clean/final_bcf/HMM_outputs/logs"
CODE="/share/maize/ntanduk/angsd_genotyping/HMM_introgression/compute_state_sites.py"

mkdir -p "$OUTDIR" "$LOGDIR"

LIST="${OUTDIR}/statepath_files.list"
if [[ ! -f "$LIST" ]]; then
  ls -1 "${STATEDIR}"/*.statepath.tsv.gz | sort > "$LIST"
fi

IN="$(sed -n "${LSB_JOBINDEX}p" "$LIST")"
if [[ -z "$IN" ]]; then
  echo "ERROR: No input file for JOBINDEX=${LSB_JOBINDEX} in $LIST" >&2
  exit 2
fi

base="$(basename "$IN")"
sample="${base%%.statepath.tsv.gz}"
OUT="${OUTDIR}/${sample}.stats.tsv"

CHR_LEN="chr1=308452471,chr2=243675191,chr3=238017767,chr4=250330460,chr5=226353449,chr6=181357234,chr7=185808916,chr8=182411202,chr9=163004744,chr10=152435371"

python3 "$CODE" --statepath_gz "$IN" --chr_len "$CHR_LEN" --out_tsv "$OUT"
echo "DONE: $OUT"
