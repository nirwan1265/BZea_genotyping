#!/usr/bin/env python3
"""
31b: one trait's k-mer ancestry GWAS on the combined BLUEs and on each year's BLUEs (README section 3.6.1).

Draws the GridLMM scans one above the other (same x axis), each with its own permutation threshold, and marks
the joint-model QTL of each run (filled = holds after Bonferroni over that run's QTL, open = does not).
Usage:
  python final_script/31b_year_compare.py --trait DTS --base data/introgression_spread --out figures/DTS_by_year.png
  (reads assoc_gridlmm_spats{,_y2023,_y2025}/scan_<T>.tsv, summary_<T>.tsv and joint_qtl/noperm{,_y2023,_y2025}/qtl2_<T>_qtl.tsv)
Needs: numpy, pandas, matplotlib.
"""
import argparse
import math
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--trait", default="DTS")
ap.add_argument("--base", default="data/introgression_spread")
ap.add_argument("--out", required=True)
ap.add_argument("--runs", default=":Both years (combined BLUEs),_y2023:2023 only,_y2025:2025 only")
a = ap.parse_args()
runs = [r.split(":", 1) for r in a.runs.split(",")]

CHR = [f"chr{i}" for i in range(1, 11)]
fig, axes = plt.subplots(len(runs), 1, figsize=(13, 2.6 * len(runs)), sharex=True)
offs = None
for ax, (suf, label) in zip(axes, runs):
    s = pd.read_csv(f"{a.base}/assoc_gridlmm_spats{suf}/scan_{a.trait}.tsv", sep="\t")
    thr = float(pd.read_csv(f"{a.base}/assoc_gridlmm_spats{suf}/summary_{a.trait}.tsv", sep="\t").threshold.iloc[0])
    q = pd.read_csv(f"{a.base}/joint_qtl/noperm{suf}/qtl2_{a.trait}_qtl.tsv", sep="\t")
    if offs is None:
        ln = s.groupby("chr").pos.max().reindex(CHR)
        offs = dict(zip(CHR, np.r_[0, np.cumsum(ln.values[:-1])]))
        mids = [offs[c] + ln[c] / 2 for c in CHR]
    s["x"] = s.pos + s.chr.map(offs)
    for i, c in enumerate(CHR):
        d = s[s.chr == c]
        ax.scatter(d.x / 1e6, d.mlog10p_add, s=2, color="#3b6ea5" if i % 2 == 0 else "#9bb7d4", rasterized=True)
    ax.axhline(thr, color="#c0392b", lw=1, ls="--")
    cut = math.log10(len(q) / 0.05)
    for r in q.itertuples():
        x = (r.pos_Mb * 1e6 + offs[r.chr]) / 1e6
        hold = r.mlog10p_joint > cut
        ax.scatter([x], [r.gridlmm_mlog10p + 0.8], marker="v", s=40, color="#222222" if hold else "white",
                   edgecolor="#222222", zorder=3)
    held = int((q.mlog10p_joint > cut).sum())
    ax.set_title(f"{label}: {len(s.dropna()):,} windows; threshold {thr:.2f}; {len(q)} QTL, {held} hold in the joint model",
                 loc="left", fontsize=10)
    ax.set_ylabel("−log10 p")
    ax.set_ylim(0, max(14.5, s.mlog10p_add.max() + 1.5))
    ax.spines[["top", "right"]].set_visible(False)
axes[-1].set_xticks([m / 1e6 for m in mids]); axes[-1].set_xticklabels(CHR)
fig.text(0.99, 0.005, "dashed line: 5% genome-wide threshold (1,000 permutations); triangles: joint-model QTL "
         "(filled = holds, open = does not)", ha="right", fontsize=8, color="#555555")
fig.tight_layout(rect=(0, 0.02, 1, 1))
fig.savefig(a.out, dpi=200)
print("wrote", a.out)
