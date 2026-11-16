#!/bin/bash
# =============================================================================
# Association mapping (steps 25-27) on a laptop, no Slurm. One or more traits per call; run each trait
# in its own terminal to go in parallel. From the repo root:
#   bash final_script/run_assoc_local.sh DTS
#   bash final_script/run_assoc_local.sh PH EH
# Per trait: scan (~3 min) -> NPERM permutations (~1.5 min each) -> threshold -> significant regions + genes.
# The inputs (step 25) are built once; other terminals wait for them.
#
# Options (environment):
#   IN        folder with tracts_all.tsv.gz and pool_summary.tsv from the cluster
#             (default ~/Documents/Research/data/BZea/introgression_spread)
#   GFF       B73 v5 GFF3 (default ~/Documents/Research/data/GFF3/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.62.chr.gff3)
#   PHENO     phenotype CSV (default data/phenotypes/BZea_spats_BLUEs_combined.csv, step 28)
#   NPERM     permutations per trait (default 100)
#   THRESHOLD -log10 p used when a trait has no permutations (NPERM=0)
#   CORES     cores for GridLMM (default 4)
#   PY        python with pandas (default: python, else python3)
#   FIXED_COL TAXON (default) or F1FAMILY (check only)
#   MARK_GENES=1  label the known flowering genes on the Manhattan plots (default 0: no labels)
# Output in $IN/assoc_gridlmm_spats/:
#   manhattan_<TRAIT>_ALL_{add,dom}.png, summary_<TRAIT>.tsv, scan_<TRAIT>.tsv
#   genes/regions_genes_<TRAIT>_ALL.tsv   significant regions and every B73 v5 gene in them (one row per gene)
#   genes/peaks_<TRAIT>_ALL.tsv           one row per region
# =============================================================================
set -euo pipefail
(( $# )) || { echo "usage: bash final_script/run_assoc_local.sh TRAIT [TRAIT ...]"; exit 1; }
REPO=$(cd "$(dirname "$0")/.." && pwd)
IN=${IN:-$HOME/Documents/Research/data/BZea/introgression_spread}
GFF=${GFF:-$HOME/Documents/Research/data/GFF3/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.62.chr.gff3}
PHENO=${PHENO:-$REPO/data/phenotypes/BZea_spats_BLUEs_combined.csv}
NPERM=${NPERM:-100}; CORES=${CORES:-4}; MARK_GENES=${MARK_GENES:-0}
PY=${PY:-$(command -v python || command -v python3)}   # conda python first (has pandas)
export GROUP_COL=none FIXED_COL=${FIXED_COL:-TAXON} MIN_N=20 MIN_CARRIERS=10
for f in "$IN/tracts_all.tsv.gz" "$IN/pool_summary.tsv" "$GFF" "$PHENO"; do [[ -f "$f" ]] || { echo "missing: $f"; exit 1; }; done
D25=$IN/assoc_spats; D26=$IN/assoc_gridlmm_spats; mkdir -p "$D26"

# ---- step 25 once (a lock so parallel terminals don't build it at the same time) ----
LOCK=$IN/.assoc_spats.lock
until mkdir "$LOCK" 2>/dev/null; do echo "waiting for another terminal to finish the inputs ..."; sleep 15; done
trap 'rmdir "$LOCK" 2>/dev/null || true' EXIT
if [[ ! -f "$D25/pheno.csv" || "$PHENO" -nt "$D25/pheno.csv" || "$IN/tracts_all.tsv.gz" -nt "$D25/pheno.csv" ]]; then
  echo "== 25: inputs ($PHENO)"
  "$PY" "$REPO/final_script/25_assoc_inputs.py" --tracts "$IN/tracts_all.tsv.gz" --pool_summary "$IN/pool_summary.tsv" \
    --pheno "$PHENO" --out_dir "$D25" --min_carriers 20 --genes "$REPO/data/flowering_genes.tsv" --gff "$GFF"
fi
if [[ "$MARK_GENES" == 1 ]]; then
  if [[ -f "$D25/gene_positions.tsv.off" ]]; then mv "$D25/gene_positions.tsv.off" "$D25/gene_positions.tsv"; fi
else
  if [[ -f "$D25/gene_positions.tsv" ]]; then mv "$D25/gene_positions.tsv" "$D25/gene_positions.tsv.off"; fi
fi
rmdir "$LOCK"; trap - EXIT

for T in "$@"; do
  head -1 "$D25/pheno.csv" | tr ',' '\n' | grep -qx "$T" || { echo "$T: not a column in $PHENO, skipped"; continue; }
  echo "== $T: scan ($FIXED_COL fixed + LOCO K)"
  Rscript "$REPO/final_script/26_assoc_gridlmm.R" scan "$D25" "$T" "$D26" "$CORES"
  if (( NPERM > 0 )); then
    echo "== $T: $NPERM permutations"
    Rscript "$REPO/final_script/26_assoc_gridlmm.R" perm "$D25" "$T" "$D26" "$NPERM" 1 "$CORES"
  fi
  Rscript "$REPO/final_script/26_assoc_gridlmm.R" summary "$D25" "$T" "$D26"
  thr=$(awk -F'\t' 'NR>1 && $2=="ALL"{print $6}' "$D26/summary_${T}.tsv")
  if [[ -z "$thr" || "$thr" == "NA" ]]; then thr=${THRESHOLD:-}; fi
  if [[ -n "$thr" ]]; then
    echo "== $T: significant regions and genes (-log10 p > $thr)"
    "$PY" "$REPO/final_script/27_assoc_genes.py" --scan "$D26/scan_${T}.tsv" --threshold "$thr" --gff "$GFF" \
      --out_dir "$D26/genes" --tag "${T}_ALL"
  else
    echo "$T: no threshold (NPERM=0 and no THRESHOLD), regions not written"
  fi
done
echo "done: $D26"
