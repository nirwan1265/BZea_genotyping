#!/bin/bash
#SBATCH --job-name=sourmash_sketch
#SBATCH --output=logs/32a_sourmash_sketch.%A_%a.out
#SBATCH --error=logs/32a_sourmash_sketch.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH --mem=8G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# 32a — k-mer sketches (sourmash FracMinHash, k = 31, scaled = SOURMASH_SCALED) for population structure.
#   Task 1          : B73 (${KMER_DIR}/B73_chr1-10.fa)
#   Tasks 2..9      : the teosinte assemblies in DONOR_ASSEMBLIES (NAME:TAXON:FASTA)
#   Tasks 10..      : pools (control + donor in pools.tsv), SOURMASH_POOLS_PER_TASK per task, from the trimmed reads
#   A sketch keeps every k-mer whose hash falls in the lowest 1/scaled of the hash space, the same k-mers in every
#   sample, so sketches can be compared directly (Brown & Irber 2016; Irber et al. 2022).
#   Needs: sourmash (pip install sourmash).  Submit (all at once):
#     N=$(awk -F'\t' 'NR>1 && ($2=="control"||$2=="donor")' ${D11}/pools.tsv | wc -l)
#     sbatch --array=1-$(( 9 + (N + 9) / 10 ))%60 32a_sourmash_sketch.sh
#   Output: ${D32}/refs/<NAME>.sig.gz, ${D32}/pools/<SAMPLE>.sig.gz. Done sketches are skipped.
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"
need sourmash
IDX=${SLURM_ARRAY_TASK_ID:-${1:?usage: bash $0 <task index>}}
PER=${SOURMASH_POOLS_PER_TASK:-10}
P="k=31,scaled=${SOURMASH_SCALED},noabund"
mkdir -p "${D32}/refs" "${D32}/pools"
NREF=$(( 1 + ${#DONOR_ASSEMBLIES[@]} ))

if (( IDX <= NREF )); then
  if (( IDX == 1 )); then NAME=B73; FA=${KMER_DIR}/B73_chr1-10.fa
  else IFS=: read -r NAME _ FA <<< "${DONOR_ASSEMBLIES[$(( IDX - 2 ))]}"; fi
  out=${D32}/refs/${NAME}.sig.gz
  [[ -s $out ]] && { echo "${NAME}: done already"; exit 0; }
  sourmash sketch dna -p "$P" --name "$NAME" -o "$out.tmp.sig.gz" "$FA" && mv "$out.tmp.sig.gz" "$out"
  echo "${NAME}: sketched"; exit 0
fi

mapfile -t POOLS < <(awk -F'\t' 'NR > 1 && ($2 == "control" || $2 == "donor") {print $1}' "${D11}/pools.tsv")
FIRST=$(( (IDX - NREF - 1) * PER ))
(( FIRST < ${#POOLS[@]} )) || { echo "task ${IDX}: no pools"; exit 0; }
for S in "${POOLS[@]:$FIRST:$PER}"; do
  out=${D32}/pools/${S}.sig.gz
  [[ -s $out ]] && { echo "${S}: done already"; continue; }
  plate=${S%%_*}; plate=${plate#PN}
  R1=${FASTQ_DIR}/filtered_S${plate}/${S}/${S}_paired_R1.fq.gz; R2=${R1/_R1/_R2}
  [[ -s $R1 && -s $R2 ]] || { echo "${S}: MISSING FASTQ ${R1}"; continue; }
  sourmash sketch dna -p "$P" --name "$S" -o "$out.tmp.sig.gz" "$R1" "$R2" && mv "$out.tmp.sig.gz" "$out"
  echo "${S}: sketched"
done
