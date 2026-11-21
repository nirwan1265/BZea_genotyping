#!/usr/bin/env python3
"""
36b — One QTL, one accession at a time: does the teosinte effect follow the environment of the donor's home site?

For every donor accession with at least --min_carriers plots carrying teosinte DNA at the QTL window (dosage > 0) and
at least --min_noncarriers plots without it, the effect = mean trait of its carrier plots - mean trait of its
non-carrier plots (same accession, so the accession's background and site are the same on both sides).
This effect is then correlated with each environment layer of the accession's site (Spearman), across all accessions
and within taxa, or within taxa and mexicana races with --group_col SUBTAXON (effect and layer centred on their group means). Environment: data/env/BZea_accession_env.csv
(aridity = precipitation / potential evapotranspiration, so HIGHER = MORE HUMID; nitrogen_0_5cm = SoilGrids topsoil
total nitrogen, cg/kg). Complements step 36 (test A: the same question as a dosage x environment regression).

Outputs: <out_prefix>.tsv (one row per accession) and <out_prefix>.png (one panel per layer; dot size = carrier plots).
Usage (from the repo root; --assoc_in = step-25 inputs, built as in the header of 36_qtl_env.py):
  python final_script/36b_qtl_env_accession.py --assoc_in $HOME/bzea_assoc_fix1 --env data/env/BZea_accession_env.csv \
     --trait NBR --window chr10:98250000 --layers aridity --group_col SUBTAXON --out_prefix data/introgression_fix1/env/NBR_Q8_accessions
Needs: numpy, pandas, scipy, matplotlib.
"""
import argparse
import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--assoc_in", required=True)
ap.add_argument("--env", required=True)
ap.add_argument("--trait", required=True)
ap.add_argument("--window", required=True, help="QTL lead window as in geno_dosage.csv, e.g. chr10:98250000")
ap.add_argument("--layers", default="aridity")
ap.add_argument("--group_col", default="TAXON", help="TAXON, or SUBTAXON to split mexicana into its four races")
ap.add_argument("--min_carriers", type=int, default=1)
ap.add_argument("--min_noncarriers", type=int, default=2)
ap.add_argument("--out_prefix", required=True)
a = ap.parse_args()
LAB = {"aridity": "aridity index of the home site (higher = more humid)", "nitrogen_0_5cm": "topsoil nitrogen at the home site (cg/kg)",
       "BIO12": "annual precipitation (mm)", "BIO18": "precipitation of warmest quarter (mm)", "elevation": "elevation (m)"}
COL = {"Zd": "#E69F00", "Zh": "#56B4E9", "Zl": "#009E73", "Zv": "#CC79A7", "Zx": "#0072B2",
       "Zx_Chal": "#D55E00", "Zx_Dura": "#0072B2", "Zx_Mesa": "#56B4E9", "Zx_Nobo": "#8C564B"}
NAME = {"Zd": "diploperennis", "Zh": "huehuetenangensis", "Zl": "luxurians", "Zv": "parviglumis", "Zx": "mexicana",
        "Zx_Chal": "mexicana, Chalco", "Zx_Dura": "mexicana, Durango", "Zx_Mesa": "mexicana, Mesa Central",
        "Zx_Nobo": "mexicana, Nobogame"}
if a.group_col == "SUBTAXON":
    COL["Zh"] = "#999999"
L = a.layers.split()

P = pd.read_csv(f"{a.assoc_in}/pheno.csv")
P = P.rename(columns={P.columns[0]: "line"}).dropna(subset=[a.trait])
P["acc"] = P.line.str.split("_").str[0]
G = pd.read_csv(f"{a.assoc_in}/geno_dosage.csv", index_col=0, usecols=lambda c: c in ("Line", a.window) or c == "")
P["dos"] = G.loc[P.line, a.window].to_numpy()
E = pd.read_csv(a.env).rename(columns={"accession-id": "acc"})
rows = []
for acc, d in P.groupby("acc"):
    c, n = d[d.dos > 0], d[d.dos == 0]
    if len(c) >= a.min_carriers and len(n) >= a.min_noncarriers:
        rows.append({"accession": acc, "taxon": d[a.group_col].iloc[0], "carrier_plots": len(c), "noncarrier_plots": len(n),
                     "effect": c[a.trait].mean() - n[a.trait].mean()})
A = pd.DataFrame(rows).merge(E[["acc"] + L].rename(columns={"acc": "accession"}), on="accession", how="left")
os.makedirs(os.path.dirname(a.out_prefix) or ".", exist_ok=True)
A.to_csv(f"{a.out_prefix}.tsv", sep="\t", index=False)

fig, axes = plt.subplots(1, len(L), figsize=(6.4 * len(L), 4.6), squeeze=False)
print(f"{a.trait} at {a.window}: {len(A)} accessions (>= {a.min_carriers} carrier and >= {a.min_noncarriers} non-carrier plots)")
for ax, lay in zip(axes[0], L):
    d = A.dropna(subset=[lay])
    r, p = spearmanr(d[lay], d.effect)
    cen = d.groupby("taxon")[[lay, "effect"]].transform(lambda s: s - s.mean())
    rw, pw = spearmanr(cen[lay], cen.effect)
    print(f"  {lay}: {len(d)} accessions, Spearman rho {r:.2f} (p {p:.2g}); within groups rho {rw:.2f} (p {pw:.2g})")
    print(d.groupby("taxon").agg(accessions=("effect", "size"), carrier_plots=("carrier_plots", "sum"),
                                 mean_effect=("effect", "mean"), median_layer=(lay, "median")).round(2).to_string())
    ax.axhline(0, color="#d6d4cc", lw=1, zorder=0)
    for t, x in d.groupby("taxon"):
        ax.scatter(x[lay], x.effect, s=18 + 14 * x.carrier_plots, color=COL.get(t, "#888888"), edgecolor="white", lw=0.8,
                   label=f"{NAME.get(t, t)} ({len(x)})", zorder=2)
    ax.set_xlabel(LAB.get(lay, lay), fontsize=9)
    ax.set_title(f"Spearman ρ = {r:.2f} (p = {p:.2g}); within {'races and taxa' if a.group_col == 'SUBTAXON' else 'taxa'} ρ = {rw:.2f} (p = {pw:.2g})", fontsize=8.5, loc="left")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=8)
axes[0, 0].set_ylabel(f"{a.trait} effect of the teosinte segment\n(carrier minus non-carrier plots, same accession)", fontsize=9)
axes[0, -1].legend(frameon=False, fontsize=8, title=("taxon or race" if a.group_col == "SUBTAXON" else "taxon") + " (accessions)", title_fontsize=8, loc="upper left", bbox_to_anchor=(1.0, 1.0))
fig.suptitle(f"{a.trait}, {a.window}: one dot per donor accession (size = carrier plots)", fontsize=10, x=0.01, ha="left")
fig.tight_layout()
fig.savefig(f"{a.out_prefix}.png", dpi=200)
print("wrote", f"{a.out_prefix}.tsv", f"{a.out_prefix}.png")
