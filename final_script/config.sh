#!/bin/bash
# =============================================================================
# config.sh — shared settings for every numbered script in final_script/.
# Edit paths HERE only. All scripts are submitted from inside final_script/:
#     cd final_script && mkdir -p logs && sbatch 07_make_site_tables.sh
# =============================================================================

# ---------- software ----------
# Batch jobs need angsd, bcftools, bgzip, python3 and Rscript on the PATH.
# Put the line(s) that load them here. Find them on the login node with:
#     which angsd bcftools bgzip Rscript
#     conda env list            # an env that has bcftools/angsd?
#     module avail 2>&1 | grep -i -E "bcftools|angsd|htslib|samtools"
# Examples (uncomment / edit ONE):
# source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate hmm
# module load bcftools angsd htslib R
# export PATH=/path/to/angsd:/path/to/bcftools/bin:$PATH

# stop early with a clear message if a tool is missing
need() {
  local t
  for t in "$@"; do
    command -v "$t" > /dev/null 2>&1 || {
      echo "ERROR: '$t' not found in PATH inside the job. Load it in the software section of config.sh." >&2
      exit 1
    }
  done
}

# ---------- inputs ----------
BAMLIST=/rsstu/users/r/rrellan/DOE_CAREER/BZea/dedup_filtered/bam_list.txt
REF=/rsstu/users/r/rrellan/sara/ref/Zm-B73-REFERENCE-NAM-5.0.fa
SAMPLE_META=/share/maize/ntanduk/BZea_genotyping/data/name_change/name_change_final.txt

# ---------- genome split ----------
CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10)
NPARTS=20          # parts per chromosome. Slurm arrays in 01 and 08 are sized 10 x 20 = 200.

# ---------- outputs ----------
WORK=/rsstu/users/r/rrellan/BZea/angsd_genotyping/final
D01=${WORK}/01_angsd_parts          # discovery BCF per part
D03=${WORK}/03_angsd_rerun          # discovery BCF per sub-part (killed parts)
D05=${WORK}/05_per_chr              # discovery BCF per chromosome
D06=${WORK}/06_per_chr_renamed      # same, samples renamed to PN*_SID*
D07=${WORK}/07_sites                # B73-polarized site tables + ANGSD -sites files
D08=${WORK}/08_counts_parts         # read counts per part
D09=${WORK}/09_counts_by_chr        # read counts per chromosome (final)
D10=${WORK}/10_qc                   # QC tables and plots

# Site discovery (01-06) has ALREADY been run. Steps 07+ read the renamed
# per-chromosome BCFs from here. If you ever rerun 01-06, set this to ${D06}.
DISCOVERY_BCF_DIR=/rsstu/users/r/rrellan/BZea/angsd_genotyping/angsd_full_snps/4_per_chr_renamed
discovery_bcf() { echo "${DISCOVERY_BCF_DIR}/BZea_${1}.full.renamed.bcf"; }

# ---------- ANGSD ----------
THREADS=${SLURM_CPUS_PER_TASK:-4}   # follows --cpus-per-task
# Read filters shared by discovery (01/03) and counting (08), so counts come from the same reads.
READ_FILTERS=(-minMapQ 30 -minQ 20 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -baq 1 -C 50)

# 01: discovery settings, as run
SNP_PVAL=1e-8
MININD=5
SETMINDEPTHIND=2

# 03: settings used for the parts that were killed and rerun in 2 halves.
# These are what was ACTUALLY run (stricter than 01). Set them equal to the 01
# values if you redo discovery from scratch. See README.
RERUN_SNP_PVAL=1e-10
RERUN_MININD=20
RERUN_SETMINDEPTHIND=4
RERUN_NSUB=2

# 08: keep ANGSD's raw A/C/G/T table (~4x larger than the AD BCF)? 1 = keep
KEEP_RAW_COUNTS=0

# =============================================================================
# helpers (no edits needed below)
# =============================================================================
FINAL_DIR="${SLURM_SUBMIT_DIR:-$PWD}"

die() { echo "ERROR: $*" >&2; exit 1; }

chr_len() {
  local len
  len=$(awk -v c="$1" '$1==c{print $2; exit}' "${REF}.fai")
  [[ -n "$len" ]] || { echo "ERROR: $1 not found in ${REF}.fai" >&2; return 1; }
  echo "$len"
}

# Slurm array index (1 .. nchrom*NPARTS) -> "CHR PART"
index_to_chr_part() {
  local i=$1 n=${#CHROMS[@]}
  if (( i < 1 || i > n * NPARTS )); then
    echo "ERROR: array index $i outside 1..$((n * NPARTS))" >&2; return 1
  fi
  echo "${CHROMS[$(( (i - 1) / NPARTS ))]} $(( (i - 1) % NPARTS + 1 ))"
}

# "CHR PART" -> "START END" (1-based, inclusive). Same split as the original 0.0 script,
# so part file names match the existing runs.
part_interval() {
  local len chunk start end
  len=$(chr_len "$1") || return 1
  chunk=$(( (len + NPARTS - 1) / NPARTS ))
  start=$(( ($2 - 1) * chunk + 1 ))
  end=$(( $2 * chunk ))
  if (( end > len )); then end=$len; fi
  echo "$start $end"
}

# PN*_SID* sample names in BAM-list order (= column order in every ANGSD output)
sample_names_from_bamlist() {
  tr -d '\r' < "$BAMLIST" | grep -v '^[[:space:]]*$' \
    | sed -E 's#^.*/##; s/\.bam$//' \
    | sed -E 's/^.*(PN[0-9]+_SID[0-9]+).*$/\1/'
}

# ANGSD writes "[ALL done]" to <prefix>.arg only when it finishes
angsd_finished() { grep -q "ALL done" "$1.arg" 2>/dev/null; }

# Parts of one chromosome in coordinate order: prints "PART START END" per line
parts_of_chr() {
  local p
  for (( p = 1; p <= NPARTS; p++ )); do echo "$p $(part_interval "$1" "$p")"; done
}
