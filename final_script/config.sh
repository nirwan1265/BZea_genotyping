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
export PATH=/rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk/seqanal/bin:$PATH   # conda env "seqanal"
# Alternatives:
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
D11=${WORK}/11_sites_filtered       # site filters, teosinte panel, donor-group read sums
D12=${WORK}/12_signal_check         # donor signal per taxon: panel vs non-panel sites
D13=${WORK}/13_pool_gl              # per-pool GL files for the HMM
D14=${WORK}/14_hmm                  # HMM state paths, tracts, summaries per pool
D15=${WORK}/15_introgression        # all pools together: tracts + summary
# family-informed version (BC1 sister donor haplotypes), run alongside 13-15:
D13B=${WORK}/13b_pool_gl_all        # per-pool GLs at ALL PASS sites (panel + non-panel)
D13C=${WORK}/13c_family_donor       # donor allele probability (PALT) per BC1 family
D14B=${WORK}/14b_hmm_family         # HMM results, family-informed
D15B=${WORK}/15b_introgression      # all pools together, family-informed

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

# 08: max BAMs per ANGSD run. Empty = as many as the node's open-file limit allows
#     (all 1,535 in one run if the limit is high enough; otherwise split + merged).
BAMS_PER_RUN=

# 08: keep ANGSD's raw A/C/G/T table (~4x larger than the AD BCF)? 1 = keep
KEEP_RAW_COUNTS=0

# ---------- 11: site filters + teosinte panel ----------
# Teosinte reference panel: Schnable 2023 set, 227 teosinte samples (see QC_NOTES.md)
TEO_PANEL_DIR=/rsstu/users/r/rrellan/sara/ref/Schnable_SNPset/teosinte_only
teo_panel() { echo "${TEO_PANEL_DIR}/schnable2023_${1}.teosinte.vcf.gz"; }
MAX_CTRL_ALT=0      # remove a site if the clean B73 controls have more than this many ALT reads
MAX_DP_MULT=2.0     # remove a site if total depth > 2 x chromosome median (repeats)
MIN_DP_MULT=0.25    # remove a site if total depth < 0.25 x median (hard to map)
MAX_OTH_FRAC=0.05   # remove a site if > 5% of its reads show a third base

# ---------- 14: HMM (Nirwan's BC2S3 HMM, pooled-read emissions) ----------
MAP_DIR=${SLURM_SUBMIT_DIR:-$PWD}/genetic_maps   # NAM genetic map per chromosome (<chr>.map)
HMM_RHO=0.2                 # scales genetic-map recombination into switch probability
HMM_DECODE=posterior_hysteresis   # posterior_hysteresis | posterior | viterbi
HMM_RIGIDITY=50             # markers in a row needed before switching (posterior_hysteresis)
SAMPLES_PER_TASK=10         # pools per array task in step 14 (array size = ceil(pools / this))

# ---------- 13c-14b: family-informed HMM ----------
# Donor haplotype per BC1 family from the sister sequencing: one column per family
# (e.g. Zd.0010_P2), FORMAT GT:LLR:PRIOR:PALT:SRC. A pool's family = first two fields
# of its pedigree (Zd.0010_P2_P3_P4.5.1.1 -> Zd.0010_P2).
DONOR_ALLELES_DIR=/rsstu/users/r/rrellan/BZea/ZEAL/variant_calling/variant_calling/union
donor_alleles() { echo "${DONOR_ALLELES_DIR}/${1}.donor_alleles.vcf.gz"; }

# =============================================================================
# helpers (no edits needed below)
# =============================================================================
FINAL_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
# ---------- 16-18: k-mer windows (alignment-free introgression signal) ----------
# For each pool: count 31-mers in its trimmed reads (KMC), keep the ones that occur once
# in the whole B73 genome, and count how many of each 100-kb window's B73 k-mers the pool
# has (KCFtools). Donor DNA lacks many B73 k-mers, so the count drops in introgressions.
# Run with the k-mer conda env active (kmc, kmc_tools, kcftools, samtools); jobs inherit it.
KMER_DIR=/share/maize/ntanduk/kmer                        # B73 k-mer databases, per-pool outputs
FASTQ_DIR=/rsstu/users/r/rrellan/sara/BZea/filtered_S     # filtered_S<plate>/<SAMPLE>/<SAMPLE>_paired_R1/R2.fq.gz
D16=${KMER_DIR}/windows                                   # per-pool window tables: <SAMPLE>.tsv.gz
KMER_K=31
KMER_WINDOW=100000          # window size (bp)
KMER_POOLS_PER_TASK=4       # pools per array task in step 17
KMER_KEEP_KMC=0             # 1 = keep each pool's KMC read database (~2-3 GB per pool) in ${KMER_DIR}/kmc

# Teosinte genome assemblies (PanAnd) for the expected DONOR level per window (step 19).
# Format NAME:TAXON:FASTA. Each is scored like a pool: which B73 single-copy k-mers it has.
ASM=/rsstu/users/r/rrellan/sara/ref/combined_genome
DONOR_ASSEMBLIES=(
  "Zd_Gigi:Zd:${ASM}/diplo/Zd-Gigi-REFERENCE-PanAnd-1.0_rename_chr.fa"
  "Zd_Momo:Zd:${ASM}/diplo/Zd-Momo-REFERENCE-PanAnd-1.0_rename_chr.fa"
  "Zh_RIMHU001:Zh:${ASM}/huehuetenagensis/Zh-RIMHU001-REFERENCE-PanAnd-1.0_rename.fa"
  "Zl_RIL003:Zl:${ASM}/luxurians/Zl-RIL003-REFERENCE-PanAnd-1.0_rename.fa"
  "Zv_TIL01:Zv:${ASM}/parviglumis/Zv-TIL01-REFERENCE-PanAnd-1.0_rename.fa"
  "Zv_TIL11:Zv:${ASM}/parviglumis/Zv-TIL11-REFERENCE-PanAnd-1.0_rename.fa"
  "Zx_TIL18:Zx:${ASM}/mexicana/Zx-TIL18-REFERENCE-PanAnd-1.0_rename.fa"
  "Zx_TIL25:Zx:${ASM}/mexicana/Zx-TIL25-REFERENCE-PanAnd-1.0_rename.fa"
)
D19=${KMER_DIR}/donor_ref                                 # <NAME>.tsv.gz per assembly + B73.tsv.gz

# ---------- 20-21: HMM on k-mer windows ----------
# ROUND=0: DONOR level from the teosinte assemblies. ROUND=k>0: from the family's own pools,
# measured on round k-1 (step 20b). Set it on the sbatch line: sbatch --export=ALL,ROUND=1 ...
ROUND=${ROUND:-0}
if (( ROUND == 0 )); then
  D20=${KMER_DIR}/hmm                 # per-pool results (<SAMPLE>.npz), all rho values
  D21=${KMER_DIR}/introgression       # final calls at the best rho
else
  D20=${KMER_DIR}/hmm_fam${ROUND}
  D21=${KMER_DIR}/introgression_fam${ROUND}
fi
# SPREAD=1 (the default, current method since 15 Nov 2026): widen the DONOR/SEG levels by how much the
# taxon's assemblies disagree per window (20_kmer_hmm.py --d_spread), because our donor plants differ
# from the assembly plants. Output: hmm_spread / introgression_spread.
# SPREAD=0 reproduces the first method (round 0, hmm / introgression): sbatch --export=ALL,SPREAD=0 ...
if (( ROUND == 0 )); then SPREAD=${SPREAD:-1}; else SPREAD=${SPREAD:-0}; fi
RUN_TAG=round${ROUND}
if (( SPREAD == 1 )); then
  (( ROUND == 0 )) || { echo "SPREAD=1 is for ROUND=0 (assembly levels) only" >&2; exit 1; }
  D20=${KMER_DIR}/hmm_spread; D21=${KMER_DIR}/introgression_spread; RUN_TAG=spread
fi
D20B=${KMER_DIR}/family_levels_fam${ROUND}   # family DONOR levels used in round ROUND (made by 20b)
KMER_HMM_POOLS_PER_TASK=30          # array size = ceil(pools / this)
KMER_PLANTS=6                       # plants per pool (SEG pools: donor copies ~ Binomial(2n, 1/2)); comma list if it varies
KMER_SKIP_GROUPS="Zm"               # taxa not run in steps 20-21 (space-separated); Zm has no assembly yet
KMER_NOISY_SIGMA=0.15               # pools with sigma above this are flagged NOISY, not used to fit rho

SAMPLE_LISTS=${FINAL_DIR}/sample_lists   # clean controls + exclusions (see QC_NOTES.md)

# ---- Association mapping (steps 25-26) ----
PHENO=${PHENO:-${FINAL_DIR}/../data/phenotypes/BZea_spats_BLUEs_combined.csv}   # step 28: plot name, traits, Sequencing_ID (2023+2025)
# old flowering-time file: ${FINAL_DIR}/../data/phenotypes/BZea_FloweringTime_spats_BLUE_weighted.csv
TRAITS=${TRAITS:-"DTS DTA ASI PH EH EN BW BL SL SPAD SPAD2 Prolif Kinki NBR LAE StPu StPi"}   # trait columns in PHENO (ST left out: ~all 0)
ASSOC_TAG=${ASSOC_TAG:-}             # optional run name, e.g. y2023 (one year's BLUEs): own folders assoc_<tag>, ..._<tag>
TAGSUF=${ASSOC_TAG:+_$ASSOC_TAG}
D25=${D21}/assoc${TAGSUF}           # step-25 inputs (dosage, SEG indicator, LOCO K, pheno)
# D26 is set below, after ASSOC_FIXED_COL (one folder per fixed effect)
D29=${D21}/joint_qtl                # step-29 joint multi-QTL model with per-taxon effects
D32=${KMER_DIR}/sourmash           # step-32 k-mer sketches (sourmash) for population structure
SOURMASH_SCALED=1000                # keep 1 in 1000 k-mer hashes (FracMinHash)
ASSOC_MIN_CARRIERS=20               # test only windows where at least this many plots carry donor (all plots)
ASSOC_GROUP_COL=none                # extra scans per group: TAXON (Zd/Zl/Zv/Zx), F1FAMILY, or none
ASSOC_FIXED_COL=${ASSOC_FIXED_COL:-TAXON}   # fixed effect in the all-plots scan: TAXON (default), SUBTAXON (mexicana split
                                    # by race), or F1FAMILY (check only; over-corrects). Override: sbatch --export=ALL,ASSOC_FIXED_COL=SUBTAXON
ASSOC_SUF=$( [[ "$ASSOC_FIXED_COL" == TAXON ]] && echo "" || echo "_${ASSOC_FIXED_COL,,}" )
D26=${D21}/assoc_gridlmm_spats${ASSOC_SUF}${TAGSUF}   # step-26 scans, permutations, Manhattan plots (two-year BLUEs; _subtaxon etc.)
ASSOC_MIN_N=20                      # groups with fewer plots are not scanned
ASSOC_GROUP_MIN_CARRIERS=10         # per-group scans: at least this many carrier plots in the group
ASSOC_PERM_PER_TASK=10              # permutations per array task: sbatch --array=1-100 -> 1,000
GFF_GLOB="/rsstu/users/r/rrellan/sara/ref/gff3/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.57.chromosome.*.gff3"   # B73 v5 genes (steps 25, 27)
METADATA=${FINAL_DIR}/../data/metadata_BZea.csv   # accession -> race (taxa_code), elevation (SUBTAXON column)
MARK_GENES=${MARK_GENES:-0}         # 1 = label the genes in GENES_TSV on the Manhattan plots (default: no labels)
GENES_TSV=${FINAL_DIR}/../data/flowering_genes.tsv   # genes to label when MARK_GENES=1 (NAME, GENE_ID)
# R needs GridLMM (github.com/deruncie/GridLMM), data.table and ggplot2. If R lives in another env, put it on PATH:
# export PATH=/path/to/r_env/bin:$PATH

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
