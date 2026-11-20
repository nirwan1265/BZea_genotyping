#!/bin/bash
# new_model/fix1 — settings shared by the fx1_* scripts (sourced AFTER final_script/config.sh).
# Nothing here changes final_script/ or the current run; all outputs go to ${FX}.
FX=${FX:-${KMER_DIR}/new_model/fix1}                     # all outputs of the full build
PILOT_SETS=${KMER_DIR}/new_model/pilot_fix1/sets         # Zv / Zx informative sets of the pilot (RULE=all), reused
KMER_ENV=/rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk/kmer_env            # kmc, kmc_tools, kcftools
PY_ENV=/rsstu/users/r/rrellan/sara/nirwan_backup/ntanduk/machine_learning      # numpy, scipy, pandas, matplotlib
FX_TAXA=(Zd Zh Zl Zv Zx)
# optional second run of steps 4-6 (fx1_submit_levels.sh): FX_RUN names its subfolder of ${FX}, FX_LEVEL_DIR gives
# DONOR levels for the windows that keep the old counts (e.g. new_model/20c's taxon levels)
FX_RUN=${FX_RUN:-}; FX_LEVEL_DIR=${FX_LEVEL_DIR:-}
FXR=${FX}${FX_RUN:+/${FX_RUN}}                           # outputs of steps 4-6 (hmm/, introgression/, checks)
FX_POOLS_PER_TASK=2                                      # pools per array task in fx1_2 (array = ceil(pools / 2))
# fx1_3: a window of a taxon switches to the informative k-mers only if it has >= FX_MIN_NEW of them in B73,
# both contrasts can be measured (>= FX_MIN_POOLS pools called DONOR there, >= 5 B73 pools/controls),
# and the new contrast is >= FX_GAIN x the old one.
FX_MIN_NEW=200
FX_MIN_POOLS=5
FX_GAIN=1.2
use_kmer_env() { export PATH=${KMER_ENV}/bin:$PATH; }
use_py_env()   { export PATH=${PY_ENV}/bin:$PATH; }
# pools of fx1_2, in a fixed order: the 10 B73 controls first (they need every taxon's set), then the donor pools
fx_pools() {
  awk -F'\t' 'NR > 1 && $2 == "control" {print $1 "\tcontrol"}' "${D11}/pools.tsv"
  awk -F'\t' -v skip=" ${KMER_SKIP_GROUPS} " 'NR > 1 && $2 == "donor" && index(skip, " " $3 " ") == 0 {print $1 "\t" $3}' "${D11}/pools.tsv"
}
