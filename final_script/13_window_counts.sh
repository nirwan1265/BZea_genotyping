#!/bin/bash
#SBATCH --job-name=win_counts
#SBATCH --array=1-10
#SBATCH --output=logs/13_win_counts.%A_%a.out
#SBATCH --error=logs/13_win_counts.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mem=4G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# 13 — B73 / ALT reads per pool in fixed windows (WIN_SIZE, default 1 Mb),
#      on the PASS teosinte-panel sites from step 11. One job per chromosome.
#
#   Output: ${D13}/BZea_<chr>.win_counts.tsv.gz
#     CHROM START END NSITES SAMPLE REF ALT
#     (one row per window x pool; windows without any site are left out)
#
#   Submit after 11:  sbatch 13_window_counts.sh
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"
need bcftools

IDX=${SLURM_ARRAY_TASK_ID:-${1:?usage: bash $0 <chrom index 1-10>}}
CHR=${CHROMS[$(( IDX - 1 ))]}
BCF=${D11}/BZea_${CHR}.PASS.panel.bcf
[[ -s $BCF ]] || die "missing ${BCF} (run 11)"
LEN=$(chr_len "$CHR")
mkdir -p "$D13"
OUT=${D13}/BZea_${CHR}.win_counts.tsv.gz
NAMES=${D13}/.names_${CHR}.txt
bcftools query -l "$BCF" > "$NAMES"
echo "${CHR}: $(wc -l < "$NAMES") pools, window ${WIN_SIZE} bp, $(bcftools index -n "$BCF") sites"

# positions are sorted, so each window is summed and printed as soon as the next one starts
bcftools query -f '%POS[\t%AD{0}\t%AD{1}]\n' "$BCF" \
  | awk -F'\t' -v OFS='\t' -v W="$WIN_SIZE" -v len="$LEN" -v chr="$CHR" -v names="$NAMES" '
      BEGIN {
        while ((getline s < names) > 0) name[++n] = s
        print "CHROM", "START", "END", "NSITES", "SAMPLE", "REF", "ALT"
        cur = -1
      }
      function flush(   i, e) {
        if (cur < 0) return
        e = (cur + 1) * W; if (e > len) e = len
        for (i = 1; i <= n; i++) print chr, cur * W + 1, e, ns, name[i], r[i] + 0, a[i] + 0
        delete r; delete a; ns = 0
      }
      {
        w = int(($1 - 1) / W)
        if (w != cur) { flush(); cur = w }
        ns++
        for (i = 1; i <= n; i++) { r[i] += $(2 * i); a[i] += $(2 * i + 1) }
      }
      END { flush() }' \
  | gzip -c > "$OUT"
rm -f "$NAMES"

echo "windows with sites: $(gzip -dc "$OUT" | awk -F'\t' 'NR > 1 {w[$2]} END {print length(w)}')"
echo "WROTE: ${OUT}"
