#!/bin/bash
# =============================================================================
# More permutations for the k-mer ancestry GWAS (step 26, GridLMM) on a laptop, then the new threshold,
# peaks and genes (steps 26 summary + 27). Permutations are added in chunks of CHUNK (100), one file per seed
# (perm_<TRAIT>_seed<N>.tsv); the summary pools every seed file. The first 100 (seed 1, from run_assoc_local.sh)
# are kept, so going to 1,000 means adding seeds 2..10. Seeds already done are skipped, so a stopped run can
# simply be restarted.
#   bash final_script/run_perm_local.sh DTS            # one trait
#   bash final_script/run_perm_local.sh DTS DTA PH     # several, one after another
# Options (environment):
#   TOTAL (1000)  permutations wanted in total      CHUNK (100)  permutations per seed file
#   JOBS (4)      seed chunks run at the same time   CORES (2)    cores per chunk   (JOBS x CORES <= your cores)
#   IN            default ~/Documents/Research/data/BZea/introgression_spread
#   GFF           default ~/Documents/Research/data/GFF3/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.62.chr.gff3
# =============================================================================
set -euo pipefail
(( $# )) || { echo "usage: bash final_script/run_perm_local.sh TRAIT [TRAIT ...]"; exit 1; }
REPO=$(cd "$(dirname "$0")/.." && pwd)
IN=${IN:-$HOME/Documents/Research/data/BZea/introgression_spread}
GFF=${GFF:-$HOME/Documents/Research/data/GFF3/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.62.chr.gff3}
TOTAL=${TOTAL:-1000}; CHUNK=${CHUNK:-100}; JOBS=${JOBS:-4}; CORES=${CORES:-2}
PY=${PY:-$(command -v python || command -v python3)}
export GROUP_COL=none FIXED_COL=${FIXED_COL:-TAXON} MIN_N=20 MIN_CARRIERS=10
D25=$IN/assoc_spats; D26=$IN/assoc_gridlmm_spats
[[ -f "$D25/pheno.csv" ]] || { echo "missing $D25/pheno.csv: run run_assoc_local.sh first"; exit 1; }
NSEED=$(( (TOTAL + CHUNK - 1) / CHUNK ))

for T in "$@"; do
  todo=()
  for s in $(seq 1 $NSEED); do [[ -s "$D26/perm_${T}_seed${s}.tsv" ]] || todo+=("$s"); done
  echo "== $T: ${#todo[@]} of $NSEED chunks of $CHUNK permutations to run (JOBS=$JOBS, CORES=$CORES each)"
  if (( ${#todo[@]} )); then
    printf '%s\n' "${todo[@]}" | xargs -P "$JOBS" -I{} sh -c \
      "Rscript '$REPO/final_script/26_assoc_gridlmm.R' perm '$D25' '$T' '$D26' $CHUNK {} $CORES > '$D26/perm_${T}_seed{}.log' 2>&1 && echo '  $T seed {} done'"
  fi
  Rscript "$REPO/final_script/26_assoc_gridlmm.R" summary "$D25" "$T" "$D26"
  thr=$(awk -F'\t' 'NR>1 && $2=="ALL"{print $6}' "$D26/summary_${T}.tsv")
  echo "== $T: threshold $thr from $(awk -F'\t' 'NR>1 && $2=="ALL"{print $5}' "$D26/summary_${T}.tsv") permutations"
  if [[ -z "$thr" || "$thr" == "NA" ]]; then echo "$T: no threshold (fewer than 20 permutations)"; continue; fi
  "$PY" "$REPO/final_script/27_assoc_genes.py" --scan "$D26/scan_${T}.tsv" --threshold "$thr" --gff "$GFF" \
    --out_dir "$D26/genes" --tag "${T}_ALL"
done
