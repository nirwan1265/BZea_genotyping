#!/usr/bin/env python3
"""
36e — README figure for one QTL: does the teosinte segment add more to the trait when it comes from a wetter rainy
season? Two panels, mexicana split into its four races (SUBTAXON):
  a  per taxon/race, ordered by the median rainy-season rain at the donors' home sites: mean trait of plots WITHOUT
     the segment (open) and WITH it (filled), from the same taxon/race; in grey, the same difference within BC2
     families (sibs from the same BC2 plant, carriers minus non-carriers; families with at least one of each)
  b  one dot per donor accession: its effect (carrier minus non-carrier plots, same accession; as step 36b) against
     its rainy-season rain (sum of the --months, WorldClim 2.1, step 36c); Spearman rho across accessions and within
     taxa/races (both centred on their group means)
Usage (from the repo root; step-25 inputs and 36c output first):
  python final_script/36e_qtl_rain_figure.py --assoc_in $HOME/bzea_assoc_fix1 --monthly data/env/BZea_accession_monthly.csv \
     --trait NBR --window chr10:98250000 --out figures/fix1/NBR_Q8_rain.png
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
ap.add_argument("--monthly", required=True)
ap.add_argument("--trait", required=True)
ap.add_argument("--window", required=True)
ap.add_argument("--months", default="7 8 9 10", help="months summed as the rainy season (default July-October)")
ap.add_argument("--label", default="nodes with brace roots", help="trait name on the axes")
ap.add_argument("--out", required=True)
a = ap.parse_args()
MON = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
ms = [int(m) for m in a.months.split()]
season = f"{MON[ms[0] - 1]}–{MON[ms[-1] - 1]}"
COL = {"Zd": "#E69F00", "Zh": "#999999", "Zl": "#009E73", "Zv": "#CC79A7", "Zx_Chal": "#D55E00", "Zx_Dura": "#0072B2",
       "Zx_Mesa": "#56B4E9", "Zx_Nobo": "#8C564B"}
NAME = {"Zd": "diploperennis", "Zh": "huehuetenangensis", "Zl": "luxurians", "Zv": "parviglumis", "Zx_Chal": "mexicana Chalco",
        "Zx_Dura": "mexicana Durango", "Zx_Mesa": "mexicana Mesa Central", "Zx_Nobo": "mexicana Nobogame"}

P = pd.read_csv(f"{a.assoc_in}/pheno.csv")
P = P.rename(columns={P.columns[0]: "line"}).dropna(subset=[a.trait])
P["acc"] = P.line.str.split("_").str[0]
G = pd.read_csv(f"{a.assoc_in}/geno_dosage.csv", index_col=0, usecols=lambda c: c in ("Line", a.window))
P["dos"] = G.loc[P.line, a.window].to_numpy()
M = pd.read_csv(a.monthly).rename(columns={"accession-id": "acc"})
M["rain"] = M[[f"prec_{m:02d}" for m in ms]].sum(axis=1)
P = P.merge(M[["acc", "rain"]], on="acc", how="left")
P["car"] = P.dos > 0

# a: per group
grp = []
for g, d in P.groupby("SUBTAXON"):
    grp.append({"g": g, "rain": d.drop_duplicates("acc").rain.median(), "non": d.loc[~d.car, a.trait].mean(),
                "wi": d.loc[d.car, a.trait].mean(), "n_with": int(d.car.sum()), "n_non": int((~d.car).sum())})
Gr = pd.DataFrame(grp).sort_values("rain").reset_index(drop=True)
# sibs: within each BC2 family (lines from the same BC2 plant, e.g. Zv.0200_P3_P2_P1 for Zv.0200_P3_P2_P1.1.1.1),
# carrier minus non-carrier sibs; families with at least one of each; mean over the families of the taxon/race
P["bc2"] = P.line.str.split(r"[._]").str[:5].str.join("_")
fam = []
for b, d in P.groupby("bc2"):
    c, n = d[d.car], d[~d.car]
    if len(c) and len(n):
        fam.append({"g": d.SUBTAXON.iloc[0], "diff": c[a.trait].mean() - n[a.trait].mean()})
Fm = pd.DataFrame(fam)
sib = Fm.groupby("g")["diff"].agg(["mean", "size"])
Gr["sib"] = Gr.g.map(sib["mean"]); Gr["n_fam"] = Gr.g.map(sib["size"]).fillna(0).astype(int)
from scipy.stats import wilcoxon
p_sib = wilcoxon(Fm["diff"]).pvalue
# b: per accession
acc = []
for k, d in P.groupby("acc"):
    c, n = d[d.car], d[~d.car]
    if len(c) >= 1 and len(n) >= 2:
        acc.append({"acc": k, "g": d.SUBTAXON.iloc[0], "rain": d.rain.iloc[0], "n": len(c), "effect": c[a.trait].mean() - n[a.trait].mean()})
A = pd.DataFrame(acc)
r, p = spearmanr(A.rain, A.effect)
cen = A.groupby("g")[["rain", "effect"]].transform(lambda s: s - s.mean())
rw, pw = spearmanr(cen.rain, cen.effect)

fig, ax = plt.subplots(1, 2, figsize=(14, 6.0), gridspec_kw={"width_ratios": [1, 1.15]})
y = np.arange(len(Gr))
for i, q in Gr.iterrows():
    c = COL.get(q.g, "#888")
    ax[0].plot([q.non, q.wi], [i, i], color=c, lw=2.2, alpha=0.6, zorder=1)
    ax[0].scatter(q.non, i, s=55, facecolor="white", edgecolor=c, lw=2, zorder=2)
    ax[0].scatter(q.wi, i, s=70, color=c, zorder=3)
    sibtxt = f"   sibs {q.sib:+.2f} ({q.n_fam} fam.)" if q.n_fam else "   sibs –"
    ax[0].text(max(q.non, q.wi) + 0.03, i, f"{q.wi - q.non:+.2f}  ({q.n_with} plots)", va="center", fontsize=8, color="#3d3d3a")
    ax[0].text(max(q.non, q.wi) + 0.03, i - 0.33, sibtxt.strip(), va="center", fontsize=7.5, color="#73726c", style="italic")
ax[0].set_yticks(y, [f"{NAME.get(g, g)}  ·  {rr:.0f} mm" for g, rr in zip(Gr.g, Gr.rain)], fontsize=8.5)
ax[0].set_xlabel(f"mean {a.label} (plot BLUE)", fontsize=9)
ax[0].set_title(f"a  Plots without (○) and with (●) the segment,\n    ordered by median {season} rain at the donors' sites\n"
                f"    grey: within BC2 families (carrier minus non-carrier sibs);\n"
                f"    {len(Fm)} families, mean {Fm['diff'].mean():+.2f}, Wilcoxon p = {p_sib:.2g}",
                fontsize=9, loc="left")
ax[0].set_xlim(right=ax[0].get_xlim()[1] + 0.25)
ax[0].set_ylim(-0.75, len(Gr) - 0.5)
for g, d in A.groupby("g"):
    ax[1].scatter(d.rain, d.effect, s=18 + 14 * d.n, color=COL.get(g, "#888"), edgecolor="white", lw=0.8,
                  label=f"{NAME.get(g, g)} ({len(d)})", zorder=2)
ax[1].axhline(0, color="#d6d4cc", lw=1, zorder=0)
ax[1].set_xlabel(f"{season} rain at the donor accession's site (mm, WorldClim 2.1)", fontsize=9)
ax[1].set_ylabel(f"effect of the segment on {a.label}\n(carrier minus non-carrier plots, same accession)", fontsize=9)
ax[1].set_title(f"b  One dot per accession (size = carrier plots)\n    Spearman ρ = {r:.2f} (p = {p:.1g}); "
                f"within taxa and races ρ = {rw:.2f} (p = {pw:.2g})", fontsize=9, loc="left")
ax[1].legend(frameon=False, fontsize=8, title="taxon / race (accessions)", title_fontsize=8, loc="upper left", bbox_to_anchor=(1.0, 1.0))
for z in ax:
    z.spines[["top", "right"]].set_visible(False)
    z.tick_params(labelsize=8)
fig.tight_layout()
os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
fig.savefig(a.out, dpi=200)
print(Gr.round(2).to_string(index=False))
print(f"sibs: {len(Fm)} BC2 families, mean {Fm['diff'].mean():+.3f}, Wilcoxon p {p_sib:.2g}")
print(f"{len(A)} accessions: Spearman rho {r:.3f} (p {p:.2g}); within taxa and races rho {rw:.3f} (p {pw:.2g})")
print("wrote", a.out)
