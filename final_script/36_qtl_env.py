#!/usr/bin/env python3
"""
36 — Does a trait, or a QTL's teosinte effect, depend on the environment the donor accession comes from?

Each plot's donor accession (pedigree name, e.g. Zl.0010_P1_P2_P2.3.1.1 -> Zl.0010) is joined to the environment at
the accession's collection site (data/env/BZea_accession_env.csv: WorldClim BIO1-19, aridity index = precipitation /
potential evapotranspiration, so HIGHER = MORE HUMID; SoilGrids topsoil N, pH, ...). Each layer is standardised over
the plots used (1 SD).

Three tests, ordinary least squares with standard errors clustered by accession (the environment varies between
accessions, not plots, so the accession is the unit; kinship is not modelled, so treat p as approximate):
  A  per QTL (step 29 lead window): y ~ TAXON + dos + dos:env
     -> change of the teosinte effect per copy per SD of env, across all taxa
  A' the same within taxa: y ~ TAXON + TAXON:dos + dos:env
     -> the same, but only from differences between accessions of the same taxon (taxon climate is not used)
  B  whole genome: y ~ TAXON + donor_share + donor_share:env
     -> does more teosinte DNA from wetter (or drier, ...) sites raise the trait more?
  C  per accession: mean of y over the accession's plots vs env (Spearman), all accessions and within taxon
dos = teosinte copies at the QTL window (0/1/2, step 25); donor_share = the pool's DONOR + SEG/2 share (step 21).

Usage (from the repo root; step-25 inputs are not in git, so build them first, outside the repo):
  python final_script/25_assoc_inputs.py --tracts data/introgression_fix1/tracts_all.tsv.gz \
     --pool_summary data/introgression_fix1/pool_summary.tsv --pheno data/phenotypes/BZea_spats_BLUEs_combined.csv \
     --out_dir $HOME/bzea_assoc_fix1 --min_carriers 20 --metadata data/metadata_BZea.csv
  python final_script/36_qtl_env.py --assoc_in $HOME/bzea_assoc_fix1 \
     --qtl data/introgression_fix1/joint_qtl/noperm/qtl2_NBR_qtl.tsv --pool_summary data/introgression_fix1/pool_summary.tsv \
     --env data/env/BZea_accession_env.csv --trait NBR --out data/introgression_fix1/env/NBR_env.tsv
Needs: numpy, pandas, scipy, statsmodels.
"""
import argparse
import os
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import spearmanr

ap = argparse.ArgumentParser()
ap.add_argument("--assoc_in", required=True)
ap.add_argument("--qtl", required=True)
ap.add_argument("--pool_summary", required=True)
ap.add_argument("--env", required=True)
ap.add_argument("--trait", required=True)
ap.add_argument("--layers", default="aridity BIO12 BIO15 BIO18 nitrogen_0_5cm elevation")
ap.add_argument("--out", required=True)
a = ap.parse_args()
WHAT = {"aridity": "aridity index (higher = more humid)", "BIO12": "annual precipitation", "BIO15": "precipitation seasonality",
        "BIO18": "precipitation of warmest quarter", "nitrogen_0_5cm": "topsoil nitrogen", "elevation": "elevation"}

P = pd.read_csv(f"{a.assoc_in}/pheno.csv")
P = P.rename(columns={P.columns[0]: "line"}).dropna(subset=[a.trait])
P["acc"] = P.line.str.split("_").str[0]
E = pd.read_csv(a.env).rename(columns={"accession-id": "acc"})
P = P.merge(E[["acc"] + a.layers.split()], on="acc", how="left")
ps = pd.read_csv(a.pool_summary, sep="\t").set_index("SAMPLE")
P["donor_share"] = P.SAMPLE.map(ps.FRAC_DONOR + ps.FRAC_SEG / 2)
G = pd.read_csv(f"{a.assoc_in}/geno_dosage.csv", index_col=0)
Q = pd.read_csv(a.qtl, sep="\t")
y = a.trait
rows = []


def fit(f, d, term):
    m = smf.ols(f, d).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(d.acc)[0]})
    return m.params[term], m.bse[term], m.pvalues[term]


for L in a.layers.split():
    d0 = P.dropna(subset=[L]).copy()
    d0["env"] = (d0[L] - d0[L].mean()) / d0[L].std()
    for i, q in enumerate(Q.itertuples(), 1):
        if q.qtl not in G.columns:
            continue
        d = d0.assign(dos=G.loc[d0.line, q.qtl].to_numpy())
        car = d[d.dos > 0]
        b, se, p = fit(f"{y} ~ C(TAXON) + dos + dos:env", d, "dos:env")
        bw, sew, pw = fit(f"{y} ~ C(TAXON) + C(TAXON):dos + dos:env", d, "dos:env")
        rows.append({"test": "A per QTL", "qtl": f"Q{i} {q.qtl}", "layer": L, "what": WHAT.get(L, L), "plots": len(d),
                     "carrier_plots": len(car), "carrier_accessions": car.acc.nunique(),
                     "effect_per_copy_joint": q.effect_per_copy_joint, "change_per_copy_per_SD": round(b, 3), "se": round(se, 3),
                     "p": p, "change_within_taxa": round(bw, 3), "se_within": round(sew, 3), "p_within": pw})
    b, se, p = fit(f"{y} ~ C(TAXON) + donor_share + donor_share:env", d0, "donor_share:env")
    rows.append({"test": "B whole genome", "qtl": "donor share", "layer": L, "what": WHAT.get(L, L), "plots": len(d0),
                 "carrier_accessions": d0.acc.nunique(), "change_per_copy_per_SD": round(b, 3), "se": round(se, 3), "p": p})
    acc = d0.groupby(["acc", "TAXON"]).agg(y=(y, "mean"), env=(L, "first")).reset_index()
    r, p = spearmanr(acc.y, acc.env)
    res = acc.y - acc.groupby("TAXON").y.transform("mean")
    rw, pw = spearmanr(res, acc.env - acc.groupby("TAXON").env.transform("mean"))
    rows.append({"test": "C per accession (Spearman rho)", "qtl": "accession mean", "layer": L, "what": WHAT.get(L, L),
                 "carrier_accessions": len(acc), "change_per_copy_per_SD": round(r, 3), "p": p,
                 "change_within_taxa": round(rw, 3), "p_within": pw})
R = pd.DataFrame(rows)
os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
R.to_csv(a.out, sep="\t", index=False)
pd.set_option("display.width", 250)
print(R.drop(columns=["what"]).to_string(index=False, float_format=lambda v: f"{v:.3g}"))
