#!/usr/bin/env python
"""
38b — Figure for NBR Q6 / ZCN8 (step 38): brace-root nodes and flowering time by teosinte copies at the window,
and the effect per copy in each taxon and mexicana race.
  a  NBR by copies at the window (0 = B73, 1 = segregating pool, 2 = donor), taxon-adjusted mean +- SE over plots
  b  DTA and DTS, the same
  c  NBR effect per teosinte copy within each taxon / race (+-95% CI; step 38, y ~ TAXON + TAXON:dos + kinship)
  d  DTA and DTS, the same
Effects and p in the panel titles are from the mixed model of step 38 (y ~ TAXON + dos + LOCO kinship).

Usage (from the repo root; step-25 inputs as in step 38):
  python final_script/38b_zcn8_figure.py --assoc_in $HOME/bzea_assoc_fix1 --window chr8:126650000 \
     --tables data/introgression_fix1/zcn8 --out figures/fix1/NBR_Q6_ZCN8.png
Needs: numpy, pandas, matplotlib.
"""
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--assoc_in", required=True)
ap.add_argument("--window", default="chr8:126650000")
ap.add_argument("--tables", required=True, help="step-38 output folder (zcn8_effects.tsv, zcn8_by_taxon.tsv)")
ap.add_argument("--out", required=True)
a = ap.parse_args()

COL = {"Zd": "#E69F00", "Zh": "#6A3D9A", "Zl": "#009E73", "Zv": "#CC79A7", "Zx_Chal": "#D55E00", "Zx_Dura": "#0072B2",
       "Zx_Mesa": "#56B4E9", "Zx_Nobo": "#8C564B"}
NAME = {"Zd": "diploperennis", "Zh": "huehuetenangensis", "Zl": "luxurians", "Zv": "parviglumis", "Zx_Chal": "mexicana Chalco",
        "Zx_Dura": "mexicana Durango", "Zx_Mesa": "mexicana Mesa Central", "Zx_Nobo": "mexicana Nobogame"}
INK, MUTED, GRID = "#2b2b28", "#6b6b66", "#e4e4e0"
TCOL = {"NBR": "#3B6E8F", "DTA": "#B5502B", "DTS": "#B5502B"}

P = pd.read_csv(f"{a.assoc_in}/pheno.csv")
P = P.rename(columns={P.columns[0]: "line"})
id_col = pd.read_csv(f"{a.assoc_in}/geno_dosage.csv", nrows=0).columns[0]
G = pd.read_csv(f"{a.assoc_in}/geno_dosage.csv", usecols=[id_col, a.window], index_col=0)
P["dos"] = P.line.map(G[a.window])
P = P.dropna(subset=["dos"])
E = pd.read_csv(f"{a.tables}/zcn8_effects.tsv", sep="\t")
B = pd.read_csv(f"{a.tables}/zcn8_by_taxon.tsv", sep="\t")
eff = E[E.model.str.startswith("1 ")].set_index("trait")


def adj_means(t):
    d = P.dropna(subset=[t]).copy()
    d["r"] = d[t] - d.groupby("TAXON")[t].transform("mean") + d[t].mean()
    s = d.groupby("dos").r.agg(["mean", "sem", "size"])
    return s


def pfmt(p):
    m, e = f"{p:.1e}".split("e")
    return f"p = {p:.2g}" if p >= 1e-3 else f"p = {m}×10$^{{{int(e)}}}$"


fig, ax = plt.subplots(2, 2, figsize=(10.5, 8.2), gridspec_kw={"height_ratios": [1, 1.15]})
XL = ["0\nB73", "1\nsegregating", "2\ndonor"]

# a: NBR by copies
s = adj_means("NBR")
ax[0, 0].errorbar(s.index, s["mean"], yerr=s["sem"], color=TCOL["NBR"], marker="o", ms=8, lw=2, capsize=4)
for x, r in s.iterrows():
    ax[0, 0].text(x + 0.05, r["mean"] - r["sem"] - 0.004, f"{int(r['size'])} plots", fontsize=8, color=MUTED, va="top")
ax[0, 0].set_ylim(s["mean"].min() - s["sem"].max() - 0.03, None)
e = eff.loc["NBR"]
ax[0, 0].set_title(f"a  Brace-root nodes\n{e.beta:+.2f} nodes per teosinte copy ({pfmt(e.p)})", fontsize=9.5,
                   loc="left", color=INK)
ax[0, 0].set_ylabel("NBR (nodes), adjusted for taxon\nmean ± SE over plots", fontsize=9)

# b: flowering by copies
for t, mk, fill, dx in (("DTA", "o", True, -0.04), ("DTS", "s", False, 0.04)):
    s = adj_means(t)
    ax[0, 1].errorbar(s.index + dx, s["mean"], yerr=s["sem"], color=TCOL[t], marker=mk, ms=7, lw=2 if fill else 1.3,
                      ls="-" if fill else "--", capsize=4, mfc=TCOL[t] if fill else "white", label=t)
ea, es = eff.loc["DTA"], eff.loc["DTS"]
ax[0, 1].set_title(f"b  Flowering (days per teosinte copy)\nDTA {ea.beta:+.2f} ({pfmt(ea.p)}),  DTS {es.beta:+.2f} ({pfmt(es.p)})",
                   fontsize=9.5, loc="left", color=INK)
ax[0, 1].set_ylabel("days after planting, adjusted for taxon\nmean ± SE over plots", fontsize=9)
ax[0, 1].legend(frameon=False, fontsize=8.5, loc="upper left")
for x in ax[0]:
    x.set_xticks([0, 1, 2], XL, fontsize=8.5)
    x.set_xlim(-0.3, 2.55)
    x.set_xlabel("teosinte copies at chr8 126.6 Mb (ZCN8 window)", fontsize=9)

# c, d: per taxon / race
order = ["Zd", "Zh", "Zl", "Zv", "Zx_Chal", "Zx_Dura", "Zx_Mesa", "Zx_Nobo"]


def group_rows(t):
    b = B[(B.trait == t) & (((B.by == "TAXON") & B.group.isin(order[:4])) | ((B.by == "SUBTAXON") & B.group.isin(order[4:])))]
    return b.set_index("group").loc[order]


def forest(axx, t, dy=0.0, fill=True, label=None):
    b = group_rows(t)
    y = np.arange(len(order))[::-1] + dy
    for yi, (g, r) in zip(y, b.iterrows()):
        axx.plot([r.beta - 1.96 * r.se, r.beta + 1.96 * r.se], [yi, yi], color=COL[g], lw=2 if fill else 1.2,
                 ls="-" if fill else "--")
        axx.plot(r.beta, yi, "o" if fill else "s", ms=7 if fill else 6, color=COL[g], mfc=COL[g] if fill else "white",
                 mec=COL[g], zorder=3)
    if label:
        axx.plot([], [], "o" if fill else "s", color=MUTED, mfc=MUTED if fill else "white", label=label)
    e = eff.loc[t]
    yt = len(order) + 0.2 + dy
    axx.plot([e.beta - 1.96 * e.se, e.beta + 1.96 * e.se], [yt, yt], color=INK, lw=2 if fill else 1.2, ls="-" if fill else "--")
    axx.plot(e.beta, yt, "D" if fill else "d", ms=7, color=INK, mfc=INK if fill else "white", zorder=3)
    return b


b = forest(ax[1, 0], "NBR")
forest(ax[1, 1], "DTA", dy=0.13, label="DTA")
forest(ax[1, 1], "DTS", dy=-0.13, fill=False, label="DTS")
ylab = [f"{NAME[g]}  ({int(r.carrier_accessions)} acc., {int(r.carrier_plots)} plots)" for g, r in b.iterrows()][::-1]
for k, (axx, t) in enumerate(((ax[1, 0], "NBR"), (ax[1, 1], "DTA"))):
    axx.axvline(0, color=MUTED, lw=1)
    axx.set_yticks(list(range(len(order))) + [len(order) + 0.2], ylab + ["all carriers"] if k == 0 else
                   [""] * len(order) + [""], fontsize=8.5)
    axx.set_ylim(-0.7, len(order) + 0.9)
ax[1, 0].set_xlabel("brace-root nodes per teosinte copy (± 95% CI)", fontsize=9)
ax[1, 1].set_xlabel("days per teosinte copy (± 95% CI)", fontsize=9)
ax[1, 0].set_title("c  Brace-root nodes, effect within each taxon / race", fontsize=9.5, loc="left", color=INK)
ax[1, 1].set_title("d  Flowering, effect within each taxon / race", fontsize=9.5, loc="left", color=INK)
ax[1, 1].legend(frameon=False, fontsize=8.5, loc="upper right")
for x in ax.flat:
    x.grid(axis="x" if x in ax[1] else "y", color=GRID, lw=0.8)
    x.set_axisbelow(True)
    for sp in ("top", "right"):
        x.spines[sp].set_visible(False)
    x.tick_params(colors=INK, labelsize=8.5)
fig.suptitle("NBR Q6 (chr8 126.2–126.7 Mb, holds ZCN8/PEBP8 Zm00001eb353250): the teosinte tract adds brace-root nodes "
             "and delays flowering", fontsize=10.5, x=0.01, ha="left", color=INK)
fig.tight_layout(rect=(0, 0, 1, 0.96))
fig.savefig(a.out, dpi=200)
print("wrote", a.out)
