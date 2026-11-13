#!/bin/bash
#SBATCH --job-name=signal_check
#SBATCH --output=logs/12_signal_check.%j.out
#SBATCH --error=logs/12_signal_check.%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem=8G
##SBATCH --partition=compute   # uncomment and set if your account needs a partition
# =============================================================================
# 12 — Donor signal per taxon on teosinte-panel sites vs other sites.
#
#   Answers: are the panel sites alone enough to see introgressions for each taxon
#   (Zd, Zh, Zl, Zm, Zv, Zx)? For a typical pool it reports ALT reads per 1 Mb in a
#   B73 region vs a fully-donor region, and how clearly they separate (z_fixed).
#   z_fixed >= 10 on panel sites -> panel sites are enough for that taxon.
#
#   Output in ${D12}: taxon_signal.tsv, pool_signal.tsv, site_filter_summary.tsv
#   Submit after 11:  sbatch 12_signal_check.sh
# =============================================================================
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/config.sh"
need python3
python3 -c "import numpy, pandas" 2> /dev/null || die "python3 needs numpy and pandas"

for c in "${CHROMS[@]}"; do
  [[ -s ${D11}/BZea_${c}.pool_class.tsv ]] || die "missing ${D11}/BZea_${c}.pool_class.tsv (run 11)"
done
GENOME_MB=$(awk -v list=" ${CHROMS[*]} " 'index(list, " " $1 " ") {s += $2} END {printf "%.1f", s / 1e6}' "${REF}.fai")
echo "genome size (Mb, ${#CHROMS[@]} chromosomes): ${GENOME_MB}"
python3 "${FINAL_DIR}/12_signal_check.py" --d11 "$D11" --genome-mb "$GENOME_MB" --out "$D12"
