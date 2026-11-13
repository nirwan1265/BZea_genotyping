#!/bin/bash
#SBATCH --job-name=pool_counts
#SBATCH --array=1-10
#SBATCH --output=logs/13_pool_counts.%A_%a.out
#SBATCH --error=logs/13_pool_counts.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mem=4G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# 13 — Per-pool marker files for the HMM (replaces old 1_extract_GL_emissions.sh).
#      One job per chromosome; one file per pool, only sites where the pool has reads.
#
#   Input : ${D11}/BZea_<chr>.PASS.panel.bcf  (PASS sites that are also teosinte-panel SNPs)
#   Output: ${D13}/<chr>/<SAMPLE>.tsv.gz
#     CHROM POS B73 ALT REF_READS ALT_READS PANEL_AC PANEL_AN
#   (PANEL_AC / PANEL_AN = how many of the panel's teosinte chromosomes carry ALT;
#    the HMM uses it as the chance that a donor carries ALT at that site)
#
#   Submit after 11:  sbatch 13_extract_pool_counts.sh
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"
need bcftools

IDX=${SLURM_ARRAY_TASK_ID:-${1:?usage: bash $0 <chrom index 1-10>}}
CHR=${CHROMS[$(( IDX - 1 ))]}
BCF=${D11}/BZea_${CHR}.PASS.panel.bcf
[[ -s $BCF ]] || die "missing ${BCF} (run 11)"
OUTDIR=${D13}/${CHR}
rm -rf "$OUTDIR"; mkdir -p "$OUTDIR"
NAMES=${OUTDIR}/.names
bcftools query -l "$BCF" > "$NAMES"
ulimit -n "$(ulimit -Hn)" 2> /dev/null || true      # one open output file per pool
echo "${CHR}: $(wc -l < "$NAMES") pools, $(bcftools index -n "$BCF") sites"

bcftools query -f '%POS\t%REF\t%ALT\t%INFO/PANEL_AC\t%INFO/PANEL_AN[\t%AD{0}\t%AD{1}]\n' "$BCF" \
  | awk -F'\t' -v OFS='\t' -v chr="$CHR" -v dir="$OUTDIR" -v names="$NAMES" '
      BEGIN { while ((getline s < names) > 0) { n++; f[n] = dir "/" s ".tsv" } }
      {
        for (i = 1; i <= n; i++) {
          r = $(4 + 2 * i); a = $(5 + 2 * i)
          if (r + a > 0) print chr, $1, $2, $3, r, a, $4, $5 > f[i]
        }
      }
      END { for (i = 1; i <= n; i++) close(f[i]) }'
rm -f "$NAMES"
gzip -f "${OUTDIR}"/*.tsv
echo "pool files: $(ls "${OUTDIR}" | wc -l)"
echo "WROTE: ${OUTDIR}/<SAMPLE>.tsv.gz"
