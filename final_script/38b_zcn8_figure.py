#!/usr/bin/env python
"""
38b — README figure for NBR Q6 / ZCN8 (step 38). Top row: dose response; bottom row: the format of 36e
(NBR_Q8_rain.png, DTS_CCT10_elevation.png).
  a  NBR by teosinte copies at the window (0 = B73, 1 = segregating pool, 2 = donor): taxon-adjusted mean +- SE over
     plots; effect per copy and p from step 38 (mixed model y ~ TAXON + dos + LOCO kinship)
  b  DTA and DTS, the same
  c  NBR per taxon/race (mexicana split into its four races): mean of plots WITHOUT the segment (open) and WITH it
     (filled), same taxon/race; in grey, the same difference within BC2 families (carrier minus non-carrier sibs from
     the same BC2 plant; families with at least one of each)
  d  the same for flowering time (--trait2, default DTA), rows in the same order
  e  one dot per donor accession: its effect on flowering (x) against its effect on NBR (y), carrier minus non-carrier
     plots of the same accession (size = carrier plots); Spearman rho across accessions and within taxa/races
Usage (from the repo root; step-25 inputs as in step 38):
  python final_script/38b_zcn8_figure.py --assoc_in $HOME/bzea_assoc_fix1 --window chr8:126650000 \
     --tables data/introgression_fix1/zcn8 --out figures/fix1/NBR_Q6_ZCN8.png
Needs: numpy, pandas, scipy, matplotlib.
"""
import argparse
import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--assoc_in", required=True)
ap.add_argument("--window", default="chr8:126650000")
ap.add_argument("--trait", default="NBR")
ap.add_argument("--trait2", default="DTA")
ap.add_argument("--label", default="nodes with brace roots")
ap.add_argument("--label2", default="days to anthesis")
ap.add_argument("--title", default="ZCN8")
ap.add_argument("--tables", default="data/introgression_fix1/zcn8", help="step-38 output folder (zcn8_effects.tsv)")
ap.add_argument("--out", required=True)
a = ap.parse_args()
COL = {"Zd": "#E69F00", "Zh": "#6A3D9A", "Zl": "#009E73", "Zv": "#CC79A7", "Zx_Chal": "#D55E00", "Zx_Dura": "#0072B2",
       "Zx_Mesa": "#56B4E9", "Zx_Nobo": "#8C564B"}
NAME = {"Zd": "diploperennis", "Zh": "huehuetenangensis", "Zl": "luxurians", "Zv": "parviglumis", "Zx_Chal": "mexicana Chalco",
        "Zx_Dura": "mexicana Durango", "Zx_Mesa": "mexicana Mesa Central", "Zx_Nobo": "mexicana Nobogame"}
ORDER = ["Zx_Nobo", "Zx_Mesa", "Zx_Dura", "Zx_Chal", "Zv", "Zl", "Zh", "Zd"]      # bottom to top

P = pd.read_csv(f"{a.assoc_in}/pheno.csv")
P = P.rename(columns={P.columns[0]: "line"})
idc = pd.read_csv(f"{a.assoc_in}/geno_dosage.csv", nrows=0).columns[0]
G = pd.read_csv(f"{a.assoc_in}/geno_dosage.csv", index_col=0, usecols=[idc, a.window])
P["dos"] = P.line.map(G[a.window])
P = P.dropna(subset=["dos"])
P["acc"] = P.line.str.split("_").str[0]
P["car"] = P.dos > 0
P["bc2"] = P.line.str.split(r"[._]").str[:5].str.join("_")


def groups(t):
    d0 = P.dropna(subset=[t])
    rows = []
    for g, d in d0.groupby("SUBTAXON"):
        rows.append({"g": g, "non": d.loc[~d.car, t].mean(), "wi": d.loc[d.car, t].mean(), "n_with": int(d.car.sum()),
                     "n_acc": d.loc[d.car, "acc"].nunique()})
    Gr = pd.DataFrame(rows).set_index("g").loc[ORDER].reset_index()
    fam = []
    for b, d in d0.groupby("bc2"):
        c, n = d[d.car], d[~d.car]
        if len(c) and len(n):
            fam.append({"g": d.SUBTAXON.iloc[0], "diff": c[t].mean() - n[t].mean()})
    Fm = pd.DataFrame(fam)
    s = Fm.groupby("g")["diff"].agg(["mean", "size"])
    Gr["sib"] = Gr.g.map(s["mean"]); Gr["n_fam"] = Gr.g.map(s["size"]).fillna(0).astype(int)
    return Gr, Fm


def dumbbell(ax, t, label, letter):
    Gr, Fm = groups(t)
    p_sib = wilcoxon(Fm["diff"]).pvalue
    vals = pd.concat([Gr.non, Gr.wi, Gr.non + Gr.sib.fillna(0)])
    span = vals.max() - vals.min()
    off = 0.03 * span
    for i, q in Gr.iterrows():
        c = COL.get(q.g, "#888")
        ax.plot([q.non, q.wi], [i, i], color=c, lw=2.2, alpha=0.6, zorder=1)
        ax.scatter(q.non, i, s=55, facecolor="white", edgecolor=c, lw=2, zorder=2)
        ax.scatter(q.wi, i, s=70, color=c, zorder=3)
        ax.text(max(q.non, q.wi) + off, i, f"{q.wi - q.non:+.2f}  ({q.n_with} plots)", va="center", fontsize=8, color="#3d3d3a")
        if q.n_fam:
            ys, x1 = i - 0.34, q.non + q.sib
            ax.plot([q.non, x1], [ys, ys], color="#c3c1b8", lw=2, zorder=1)
            ax.scatter(q.non, ys, s=26, facecolor="white", edgecolor="#a8a69d", lw=1.4, zorder=2)
            ax.scatter(x1, ys, s=32, color="#a8a69d", zorder=3)
            ax.text(max(q.non, x1) + off, ys, f"sibs {q.sib:+.2f}  ({q.n_fam} fam.)", va="center", fontsize=7.5,
                    color="#73726c", style="italic")
    ax.set_yticks(np.arange(len(Gr)), [f"{NAME.get(g, g)}  ·  {n} acc." for g, n in zip(Gr.g, Gr.n_acc)], fontsize=8.5)
    ax.set_xlabel(f"mean {label} (plot BLUE)", fontsize=9)
    ax.set_title(f"{letter}  {label[0].upper() + label[1:]}: plots without (○) and with (●) the segment\n"
                 f"    at {a.window}{(' (' + a.title + ')') if a.title else ''}\n"
                 f"    grey line below: the same within BC2 families (carrier minus non-carrier sibs);\n"
                 f"    {len(Fm)} families, mean {Fm['diff'].mean():+.2f}, Wilcoxon p = {p_sib:.2g}", fontsize=9, loc="left")
    ax.set_xlim(right=ax.get_xlim()[1] + 0.3 * span)
    ax.set_ylim(-0.75, len(Gr) - 0.5)
    return Gr, Fm


fig = plt.figure(figsize=(21, 11.6))
gs = fig.add_gridspec(2, 6, height_ratios=[0.8, 1.15], hspace=0.42, wspace=1.6)
top = [fig.add_subplot(gs[0, 0:3]), fig.add_subplot(gs[0, 3:6])]
ax = [fig.add_subplot(gs[1, 0:2]), fig.add_subplot(gs[1, 2:4]), fig.add_subplot(gs[1, 4:6])]

# top row: dose response
E = pd.read_csv(f"{a.tables}/zcn8_effects.tsv", sep="\t")
eff = E[E.model.str.startswith("1 ")].set_index("trait")


def adj_means(t):
    d = P.dropna(subset=[t]).copy()
    d["r"] = d[t] - d.groupby("TAXON")[t].transform("mean") + d[t].mean()
    return d.groupby("dos").r.agg(["mean", "sem", "size"])


def pfmt(p):
    m, e = f"{p:.1e}".split("e")
    return f"p = {p:.2g}" if p >= 1e-3 else f"p = {m}×10$^{{{int(e)}}}$"


s = adj_means(a.trait)
top[0].errorbar(s.index, s["mean"], yerr=s["sem"], color="#3B6E8F", marker="o", ms=8, lw=2, capsize=4)
for x, r in s.iterrows():
    top[0].text(x + 0.05, r["mean"] - r["sem"] - 0.004, f"{int(r['size'])} plots", fontsize=8, color="#73726c", va="top")
top[0].set_ylim(s["mean"].min() - s["sem"].max() - 0.03, None)
e = eff.loc[a.trait]
top[0].set_title(f"a  {a.label[0].upper() + a.label[1:]} by teosinte copies at {a.window}\n"
                 f"    {e.beta:+.2f} nodes per teosinte copy ({pfmt(e.p)}, mixed model)", fontsize=9, loc="left")
top[0].set_ylabel(f"{a.trait}, adjusted for taxon\n(mean ± SE over plots)", fontsize=9)
for t, mk, fill, dx in ((a.trait2, "o", True, -0.04), ("DTS", "s", False, 0.04)):
    s = adj_means(t)
    top[1].errorbar(s.index + dx, s["mean"], yerr=s["sem"], color="#B5502B", marker=mk, ms=7, lw=2 if fill else 1.3,
                    ls="-" if fill else "--", capsize=4, mfc="#B5502B" if fill else "white", label=t)
ea, es = eff.loc[a.trait2], eff.loc["DTS"]
top[1].set_title(f"b  Flowering by teosinte copies at {a.window}\n    days per teosinte copy: {a.trait2} {ea.beta:+.2f} "
                 f"({pfmt(ea.p)}),  DTS {es.beta:+.2f} ({pfmt(es.p)}), mixed model", fontsize=9, loc="left")
top[1].set_ylabel("days after planting, adjusted for taxon\n(mean ± SE over plots)", fontsize=9)
top[1].legend(frameon=False, fontsize=8.5, loc="upper left")
for z in top:
    z.set_xticks([0, 1, 2], ["0\nB73", "1\nsegregating", "2\ndonor"], fontsize=8.5)
    z.set_xlim(-0.3, 2.5)
    z.set_xlabel("teosinte copies at the window (pool call)", fontsize=9)
    z.grid(axis="y", color="#e4e4e0", lw=0.8); z.set_axisbelow(True)

G1, F1 = dumbbell(ax[0], a.trait, a.label, "c")
G2, F2 = dumbbell(ax[1], a.trait2, a.label2, "d")

acc = []
for k, d in P.groupby("acc"):
    c, n = d[d.car], d[~d.car]
    if len(c) >= 1 and len(n) >= 2:
        acc.append({"acc": k, "g": d.SUBTAXON.iloc[0], "n": len(c), "e1": c[a.trait].mean() - n[a.trait].mean(),
                    "e2": c[a.trait2].mean() - n[a.trait2].mean()})
A = pd.DataFrame(acc).dropna()
r, p = spearmanr(A.e2, A.e1)
cen = A.groupby("g")[["e1", "e2"]].transform(lambda s: s - s.mean())
rw, pw = spearmanr(cen.e2, cen.e1)
for g in ORDER[::-1]:
    d = A[A.g == g]
    if len(d):
        ax[2].scatter(d.e2, d.e1, s=18 + 14 * d.n, color=COL.get(g, "#888"), edgecolor="white", lw=0.8,
                      label=f"{NAME.get(g, g)} ({len(d)})", zorder=2)
ax[2].axhline(0, color="#d6d4cc", lw=1, zorder=0)
ax[2].axvline(0, color="#d6d4cc", lw=1, zorder=0)
ax[2].set_xlabel(f"effect of the segment on {a.label2}\n(carrier minus non-carrier plots, same accession)", fontsize=9)
ax[2].set_ylabel(f"effect of the segment on {a.label}\n(carrier minus non-carrier plots, same accession)", fontsize=9)
ax[2].set_title(f"e  One dot per accession (size = carrier plots)\n    {len(A)} accessions: {(A.e1 > 0).sum()} with more brace-root "
                f"nodes, {(A.e2 > 0).sum()} with later {a.trait2}\n    Spearman ρ = {r:.2f} (p = {p:.2g}); "
                f"within taxa and races ρ = {rw:.2f} (p = {pw:.2g})", fontsize=9, loc="left")
ax[2].legend(frameon=False, fontsize=8, title="taxon / race (accessions)", title_fontsize=8, loc="upper left", bbox_to_anchor=(1.0, 1.0))
for z in ax + top:
    z.spines[["top", "right"]].set_visible(False)
    z.tick_params(labelsize=8)
fig.subplots_adjust(left=0.095, right=0.875, top=0.94, bottom=0.07)
os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
fig.savefig(a.out, dpi=200)
print(G1.round(2).to_string(index=False)); print(G2.round(2).to_string(index=False))
for t, F in ((a.trait, F1), (a.trait2, F2)):
    print(f"{t} sibs: {len(F)} BC2 families, mean {F['diff'].mean():+.3f}, Wilcoxon p {wilcoxon(F['diff']).pvalue:.2g}")
print(f"{len(A)} accessions: Spearman rho {r:.3f} (p {p:.2g}); within taxa and races rho {rw:.3f} (p {pw:.2g})")
print("wrote", a.out)
