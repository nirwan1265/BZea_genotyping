#!/bin/bash
#SBATCH --job-name=targseq_aln
#SBATCH --output=logs/34_targseq_align.%A_%a.out
#SBATCH --error=logs/34_targseq_align.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --mem=24G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# 34 — Validation data: targeted sequencing (MolBreeding panel, high depth) of 96 of our plots.
#   One array task per sample: adapter trimming (fastp), alignment to B73 v5 (bwa mem), sorted + indexed BAM.
#   No duplicate removal: the panel may be amplicon-based, where duplicates are expected.
#   Input:  ${TARGSEQ_DIR}/<NAME>_S<n>_L001_R{1,2}_001.fastq.gz, NAME e.g. Zd0010_P3_P3_P3211
#           (= plot Zd.0010_P3_P3_P3.2.1.1; step 34c links it to the pool)
#   Output: ${D34}/bam/<NAME>.bam(.bai), ${D34}/fastp/<NAME>.json
#   Needs fastp, bwa, samtools on the PATH, and a bwa index of ${REF} (built once by task 1 if missing: ~1-2 h, 24G).
#   Submit:  N=$(ls ${TARGSEQ_DIR}/*_R1_001.fastq.gz | wc -l); sbatch --array=1-$N 34_targseq_align.sh
#   (if the index is missing, run task 1 alone first: sbatch --array=1 34_targseq_align.sh)
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"
need fastp bwa samtools
TARGSEQ_DIR=${TARGSEQ_DIR:-/rsstu/users/r/rrellan/BZea/targseq/data}
D34=${D34:-${D21}/targseq}
mkdir -p "$D34/bam" "$D34/fastp"
T=${SLURM_CPUS_PER_TASK:-8}
R1=$(ls "$TARGSEQ_DIR"/*_R1_001.fastq.gz | sed -n "${SLURM_ARRAY_TASK_ID:-1}p")
[[ -n "$R1" ]] || { echo "no sample for task ${SLURM_ARRAY_TASK_ID:-1}"; exit 0; }
R2=${R1/_R1_001/_R2_001}
NAME=$(basename "$R1" | sed -E 's/_S[0-9]+_L00[0-9]_R1_001\.fastq\.gz$//')
if [[ ! -s "${REF}.bwt" ]]; then
  [[ "${SLURM_ARRAY_TASK_ID:-1}" == 1 ]] || { echo "no bwa index for ${REF}: run task 1 first"; exit 1; }
  echo "building bwa index of ${REF}"; bwa index "$REF"
fi
echo "$NAME | $R1"
fastp -i "$R1" -I "$R2" -o "$D34/fastp/$NAME.R1.fq.gz" -O "$D34/fastp/$NAME.R2.fq.gz" \
      --detect_adapter_for_pe -w "$T" -j "$D34/fastp/$NAME.json" -h /dev/null
bwa mem -t "$T" -R "@RG\tID:$NAME\tSM:$NAME\tPL:ILLUMINA" "$REF" "$D34/fastp/$NAME.R1.fq.gz" "$D34/fastp/$NAME.R2.fq.gz" |
  samtools sort -@ 4 -m 2G -o "$D34/bam/$NAME.bam" -
samtools index "$D34/bam/$NAME.bam"
rm -f "$D34/fastp/$NAME.R1.fq.gz" "$D34/fastp/$NAME.R2.fq.gz"
samtools flagstat "$D34/bam/$NAME.bam" | head -8
