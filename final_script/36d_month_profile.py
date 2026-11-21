#!/usr/bin/env python3
"""
36d — Which months of rain go with a QTL's teosinte effect? Month-by-month version of steps 36 (test A) and 36b.

For each month m (prec_m = WorldClim 2.1 monthly precipitation at the donor's home site, step 36c):
  per accession   the accession's effect (carrier minus non-carrier plots, as in 36b) vs prec_m: Spearman rho,
                  across all accessions
  within taxa     y ~ TAXON + TAXON:dos + dos:z(prec_m), errors clustered by accession (as 36 test A')
  within races    the same with SUBTAXON (mexicana split into Chalco, Durango, Mesa Central, Nobogame): only
                  differences between accessions of the same taxon or race, so lineage cannot make the pattern
z(prec_m) = prec_m standardised over the plots (1 SD). 12 months are tested: judge them against a Bonferroni
cut of p < 0.05 / 12 = 0.004, or look at the shape of the whole season, not the best single month.

Outputs: <out_prefix>.tsv (one row per month) and <out_prefix>.png (top: median monthly rain of the carrier
accessions, by taxon/race; middle: per-accession rho; bottom: within-races change per SD with 95% interval).
Usage (from the repo root; step-25 inputs and 36c output first):
  python final_script/36d_month_profile.py --assoc_in $HOME/bzea_assoc_fix1 --monthly data/env/BZea_accession_monthly.csv \
     --trait NBR --window chr10:98250000 --out_prefix data/introgression_fix1/env/NBR_Q8_months
Needs: numpy, pandas, scipy, statsmodels, matplotlib.
"""
import argparse
import os
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--assoc_in", required=True)
ap.add_argument("--monthly", required=True)
ap.add_argument("--trait", required=True)
ap.add_argument("--window", required=True)
ap.add_argument("--var", default="prec", help="prec (rain, mm) or tavg (mean temperature)")
ap.add_argument("--min_carriers", type=int, default=1)
ap.add_argument("--min_noncarriers", type=int, default=2)
ap.add_argument("--out_prefix", required=True)
a = ap.parse_args()
MON = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
cols = [f"{a.var}_{m:02d}" for m in range(1, 13)]

P = pd.read_csv(f"{a.assoc_in}/pheno.csv")
P = P.rename(columns={P.columns[0]: "line"}).dropna(subset=[a.trait])
P["acc"] = P.line.str.split("_").str[0]
G = pd.read_csv(f"{a.assoc_in}/geno_dosage.csv", index_col=0, usecols=lambda c: c in ("Line", a.window))
P["dos"] = G.loc[P.line, a.window].to_numpy()
M = pd.read_csv(a.monthly).rename(columns={"accession-id": "acc"})[["acc"] + cols]
P = P.merge(M, on="acc", how="left")

acc = []
for k, d in P.groupby("acc"):
    c, n = d[d.dos > 0], d[d.dos == 0]
    if len(c) >= a.min_carriers and len(n) >= a.min_noncarriers:
        acc.append({"acc": k, "group": d.SUBTAXON.iloc[0], "effect": c[a.trait].mean() - n[a.trait].mean(),
                    **{x: d[x].iloc[0] for x in cols}})
A = pd.DataFrame(acc)

rows = []
for i, x in enumerate(cols):
    d = P.dropna(subset=[x]).copy()
    d["env"] = (d[x] - d[x].mean()) / d[x].std()
    cl = {"groups": pd.factorize(d.acc)[0]}
    mt = smf.ols(f"{a.trait} ~ C(TAXON) + C(TAXON):dos + dos:env", d).fit(cov_type="cluster", cov_kwds=cl)
    mr = smf.ols(f"{a.trait} ~ C(SUBTAXON) + C(SUBTAXON):dos + dos:env", d).fit(cov_type="cluster", cov_kwds=cl)
    aa = A.dropna(subset=[x])
    r, p = spearmanr(aa[x], aa.effect)
    rows.append({"month": MON[i], "layer": x, "median_value": round(d.drop_duplicates("acc")[x].median(), 1),
                 "accessions": len(aa), "rho_per_accession": round(r, 3), "p_rho": p,
                 "within_taxa": round(mt.params["dos:env"], 3), "se_taxa": round(mt.bse["dos:env"], 3), "p_taxa": mt.pvalues["dos:env"],
                 "within_races": round(mr.params["dos:env"], 3), "se_races": round(mr.bse["dos:env"], 3), "p_races": mr.pvalues["dos:env"]})
R = pd.DataFrame(rows)
os.makedirs(os.path.dirname(a.out_prefix) or ".", exist_ok=True)
R.to_csv(f"{a.out_prefix}.tsv", sep="\t", index=False)
pd.set_option("display.width", 200)
print(f"{a.trait} at {a.window}, {a.var}: {len(A)} accessions; Bonferroni over 12 months: p < 0.0042")
print(R.to_string(index=False, float_format=lambda v: f"{v:.3g}"))

x = np.arange(12)
fig, ax = plt.subplots(3, 1, figsize=(9.5, 8.2), sharex=True, gridspec_kw={"height_ratios": [1.1, 1, 1]})
COL = {"Zd": "#E69F00", "Zh": "#999999", "Zl": "#009E73", "Zv": "#CC79A7", "Zx_Chal": "#D55E00", "Zx_Dura": "#0072B2",
       "Zx_Mesa": "#56B4E9", "Zx_Nobo": "#8C564B"}
NAME = {"Zd": "diploperennis", "Zh": "huehuetenangensis", "Zl": "luxurians", "Zv": "parviglumis", "Zx_Chal": "mexicana, Chalco",
        "Zx_Dura": "mexicana, Durango", "Zx_Mesa": "mexicana, Mesa Central", "Zx_Nobo": "mexicana, Nobogame"}
for g, d in A.groupby("group"):
    ax[0].plot(x, d[cols].median().to_numpy(), color=COL.get(g, "#888"), lw=2, label=f"{NAME.get(g, g)} ({len(d)})")
ax[0].set_ylabel("median monthly rain (mm)" if a.var == "prec" else "median monthly mean temp. (°C)", fontsize=9)
ax[0].legend(frameon=False, fontsize=7.5, loc="upper left", bbox_to_anchor=(1.0, 1.0))
ax[1].bar(x, R.rho_per_accession, color=["#2a78d6" if p < 0.0042 else "#b8c9e0" for p in R.p_rho])
ax[1].axhline(0, color="#999", lw=0.8)
ax[1].set_ylabel("per-accession\nSpearman ρ", fontsize=9)
ax[2].errorbar(x, R.within_races, yerr=1.96 * R.se_races, fmt="o", color="#2a78d6", ecolor="#9bbbe3", capsize=3)
for i, p in enumerate(R.p_races):
    if p < 0.0042:
        ax[2].annotate("*", (x[i], R.within_races[i] + 1.96 * R.se_races[i]), ha="center", va="bottom", fontsize=12)
ax[2].axhline(0, color="#999", lw=0.8)
ax[2].set_ylabel(f"within races: extra {a.trait}\nper copy per SD", fontsize=9)
ax[2].set_xticks(x, MON, fontsize=9)
ax[1].set_title("dark bars / *: p < 0.05 / 12", fontsize=8, loc="right")
for z in ax:
    z.spines[["top", "right"]].set_visible(False)
    z.tick_params(labelsize=8)
fig.suptitle(f"{a.trait}, {a.window}: teosinte effect vs monthly {'rain' if a.var == 'prec' else 'temperature'} at the donor's home site",
             fontsize=10, x=0.01, ha="left")
fig.tight_layout()
fig.savefig(f"{a.out_prefix}.png", dpi=200)
print("wrote", f"{a.out_prefix}.tsv", f"{a.out_prefix}.png")
