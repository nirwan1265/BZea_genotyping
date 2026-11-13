#!/bin/bash
#SBATCH --job-name=pool_gl
#SBATCH --array=1-10
#SBATCH --output=logs/13_pool_gl.%A_%a.out
#SBATCH --error=logs/13_pool_gl.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mem=8G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# 13 — Per-pool GL files for the HMM (replaces old 1_extract_GL_emissions.sh).
#      One job per chromosome; one file per pool, only sites where the pool has data.
#
#   Sites : the PASS teosinte-panel sites from step 11 (${D11}/BZea_<chr>.PASS.panel.bcf)
#   GLs   : ANGSD genotype likelihoods from the discovery BCF (DISCOVERY_BCF_DIR), no rerun.
#           ANGSD stores alleles as major,minor. Where the major allele is not the B73 base,
#           the three GLs are flipped, so GL_0 is always "no ALT", GL_2 "all ALT".
#           GL = 0,0,0 means no reads in that pool there -> left out.
#   Output: ${D13}/<chr>/<SAMPLE>.tsv.gz
#     CHROM POS B73 ALT GL_0 GL_1 GL_2 PANEL_AC PANEL_AN
#     GL_0/1/2 = log10 likelihood of the pool's reads if the ALT allele fraction is 0, 1/2, 1
#     PANEL_AC/AN = teosinte panel ALT count / chromosomes (chance a donor carries ALT)
#
#   Submit after 11:  sbatch 13_extract_pool_gl.sh
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"
need bcftools

IDX=${SLURM_ARRAY_TASK_ID:-${1:?usage: bash $0 <chrom index 1-10>}}
CHR=${CHROMS[$(( IDX - 1 ))]}
PASS=${D11}/BZea_${CHR}.PASS.panel.bcf
DISC=$(discovery_bcf "$CHR")
[[ -s $PASS ]] || die "missing ${PASS} (run 11)"
[[ -s $DISC ]] || die "missing discovery BCF ${DISC}"
[[ -s ${D11}/pools.tsv ]] || die "missing ${D11}/pools.tsv (run 11)"

OUTDIR=${D13}/${CHR}
rm -rf "$OUTDIR"; mkdir -p "$OUTDIR"
TMP=${OUTDIR}/.tmp
mkdir -p "$TMP"
trap 'rm -rf "$TMP"' EXIT
ulimit -n "$(ulimit -Hn)" 2> /dev/null || true      # one open output file per pool

# sites to use (B73 allele, ALT allele, panel AC/AN) and the pools to keep
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/PANEL_AC\t%INFO/PANEL_AN\n' "$PASS" > "${TMP}/sites.tsv"
cut -f1,2 "${TMP}/sites.tsv" > "${TMP}/targets.tsv"
awk -F'\t' 'NR > 1 && ($2 == "control" || $2 == "donor") {print $1}' "${D11}/pools.tsv" > "${TMP}/keep.txt"
bcftools view -h "$DISC" | tail -n 1 | tr '\t' '\n' | tail -n +10 | grep -Fxf "${TMP}/keep.txt" > "${TMP}/names.txt" || true
echo "${CHR}: $(wc -l < "${TMP}/sites.tsv") panel PASS sites | pools in discovery BCF and kept: $(wc -l < "${TMP}/names.txt") of $(wc -l < "${TMP}/keep.txt")"

bcftools view -T "${TMP}/targets.tsv" -S "${TMP}/names.txt" -Ou "$DISC" \
  | bcftools query -f '%POS\t%REF\t%ALT[\t%GL]\n' \
  | awk -F'\t' -v OFS='\t' -v chr="$CHR" -v dir="$OUTDIR" -v sites="${TMP}/sites.tsv" -v names="${TMP}/names.txt" '
      BEGIN {
        while ((getline l < sites) > 0) { split(l, s, "\t"); b[s[2]] = s[3]; al[s[2]] = s[4]; ac[s[2]] = s[5]; an[s[2]] = s[6] }
        while ((getline l < names) > 0) { n++; f[n] = dir "/" l ".tsv" }
      }
      {
        p = $1
        if (!(p in b)) { nsk++; next }
        if ($2 == b[p] && $3 == al[p]) flip = 0
        else if ($3 == b[p] && $2 == al[p]) { flip = 1; nflip++ }
        else { nmis++; next }
        nsite++
        for (i = 1; i <= n; i++) {
          if (split($(3 + i), g, ",") != 3) continue
          if (g[1] == 0 && g[2] == 0 && g[3] == 0) continue          # no reads in this pool
          if (flip) print chr, p, b[p], al[p], g[3], g[2], g[1], ac[p], an[p] > f[i]
          else      print chr, p, b[p], al[p], g[1], g[2], g[3], ac[p], an[p] > f[i]
        }
      }
      END {
        for (i = 1; i <= n; i++) close(f[i])
        printf "sites used: %d (GL flipped to B73 orientation at %d) | allele mismatch: %d | not a PASS site: %d\n", nsite, nflip, nmis, nsk > "/dev/stderr"
      }'
gzip -f "${OUTDIR}"/*.tsv
echo "pool files: $(find "${OUTDIR}" -maxdepth 1 -name "*.tsv.gz" | wc -l)"
echo "WROTE: ${OUTDIR}/<SAMPLE>.tsv.gz"
