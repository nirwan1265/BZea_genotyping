#!/bin/bash
# =============================================================================
# Association mapping (steps 25-27) on a laptop, no Slurm. Run from the repo root:
#   bash final_script/run_assoc_local.sh <dir with tracts_all.tsv.gz + pool_summary.tsv> <B73 v5 GFF3> [traits]
# Example:
#   bash final_script/run_assoc_local.sh ~/Documents/Research/data/BZea/introgression_spread \
#        ~/Documents/Research/data/GFF3/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.62.chr.gff3 "DTS PH"
# Options (environment):
#   PHENO   phenotype CSV (default: data/phenotypes/BZea_spats_BLUEs_combined.csv, step 28)
#   NPERM   permutations per trait (default 0 = scans only; one permutation ~1.5 min, so 100 = ~2.5 h per trait)
#   THRESHOLD  -log10 p used for peaks/genes when a trait has no permutations (e.g. from DTS and PH permutations)
#   CORES   cores for GridLMM (default 4)
#   FIXED_COL  F1FAMILY (default) or TAXON (faster)
# Output: <dir>/assoc_spats (inputs) and <dir>/assoc_gridlmm_spats (scans, summary_<TRAIT>.tsv, Manhattan plots, genes/)
# =============================================================================
set -euo pipefail
IN=${1:?dir with tracts_all.tsv.gz and pool_summary.tsv}
GFF=${2:?B73 v5 GFF3}
REPO=$(cd "$(dirname "$0")/.." && pwd)
source <(grep '^TRAITS=' "$REPO/final_script/config.sh")
TRAITS=${3:-$TRAITS}
PHENO=${PHENO:-$REPO/data/phenotypes/BZea_spats_BLUEs_combined.csv}
NPERM=${NPERM:-0}; CORES=${CORES:-4}
export GROUP_COL=none FIXED_COL=${FIXED_COL:-F1FAMILY} MIN_N=20 MIN_CARRIERS=10
D25=$IN/assoc_spats; D26=$IN/assoc_gridlmm_spats; mkdir -p "$D26"   # new dirs: old DTS/DTA permutations are not mixed in

echo "== 25: inputs ($PHENO)"
python3 "$REPO/final_script/25_assoc_inputs.py" --tracts "$IN/tracts_all.tsv.gz" --pool_summary "$IN/pool_summary.tsv" \
  --pheno "$PHENO" --out_dir "$D25" --min_carriers 20 --genes "$REPO/data/flowering_genes.tsv" --gff "$GFF"

for T in $TRAITS; do
  echo "== $T: scan"
  Rscript "$REPO/final_script/26_assoc_gridlmm.R" scan "$D25" "$T" "$D26" "$CORES"
  if (( NPERM > 0 )); then
    echo "== $T: $NPERM permutations"
    Rscript "$REPO/final_script/26_assoc_gridlmm.R" perm "$D25" "$T" "$D26" "$NPERM" 1 "$CORES"
  fi
  Rscript "$REPO/final_script/26_assoc_gridlmm.R" summary "$D25" "$T" "$D26"
  thr=$(awk -F'\t' 'NR>1 && $2=="ALL"{print $6}' "$D26/summary_${T}.tsv")
  [[ -z "$thr" || "$thr" == "NA" ]] && thr=${THRESHOLD:-}
  if [[ -n "$thr" && "$thr" != "NA" ]]; then
    python3 "$REPO/final_script/27_assoc_genes.py" --scan "$D26/scan_${T}.tsv" --threshold "$thr" --gff "$GFF" \
      --out_dir "$D26/genes" --tag "${T}_ALL" --find 'CCT'
  fi
done
echo "done: $D26"
