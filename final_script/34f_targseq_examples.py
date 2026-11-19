#!/usr/bin/env python3
"""
34f — Examples: our k-mer calls against the targeted sequencing, one chromosome per panel (README section 2.9.1).

Each panel has three rows, top to bottom:
  1. "k-mer call (pool)":       our HMM call of the 6-plant pool from 0.8x skim-seq (B73 / SEG / DONOR)
  2. "targeted seq. call":      the targeted-sequencing call per targeted region (SNPs within 200 kb of each other):
                                AB if its donor-carrying SNPs have a median fraction below 0.8, BB if >= 0.8;
                                a short light tick = "no donor allele seen" (fewer than 2 SNPs carry donor alleles:
                                AA, or the donor carries the B73 allele at these SNPs, so the region says nothing)
  3. "targeted seq. SNPs":      the donor allele fraction at each real SNP (homozygous alt in another sample, as in 34d)
Same colours in rows 1 and 2: blue = SEG / AB, orange = DONOR / BB; grey = k-mer B73.
The outlined tract is the one described in the title. The three bad samples of 2.9.1 are not used.

Usage: python final_script/34f_targseq_examples.py --sites data/introgression_spread/targseq/sites_used.tsv.gz \
         --tracts data/introgression_spread/tracts_all.tsv.gz --out figures/targseq_examples.png \
         [--example SAMPLE:CHR[@FROM-TO]:START_MB:END_MB:"title" ...]   (@FROM-TO: zoom the panel to FROM-TO Mb)
Needs: numpy, pandas, matplotlib.
"""
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

EX = [
    "Zd0010_P3_P3_P3411:chr10:81.28:152.45:we call SEG, the targeted sequencing says AB (agree)",
    "Zx0030_P2_P2_P2111:chr3:70.58:81.28:the chromosome agrees, except the outlined SEG tract between two DONOR tracts: it is BB (our error)",
    "Zx0510_P1_P3_P2111:chr4:149.16:219.23:we call SEG, the targeted sequencing says AB with a BB stretch at 192–200 Mb",
    "Zv0030_P1_P3_P1111:chr1@280-309:297.01:299.21:parviglumis, zoomed to 280–309 Mb: SEG between two DONOR tracts, BB in the targeted sequencing (our error)",
]
ap = argparse.ArgumentParser()
ap.add_argument("--sites", required=True)
ap.add_argument("--tracts", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--example", action="append", default=None)
a = ap.parse_args()
ex = [e.split(":", 4) for e in (a.example or EX)]
GREY, BLUE, ORANGE, INK, MUTED, PALE = "#cfcdc4", "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#a8a69e"
COL = {"B73": GREY, "SEG": BLUE, "DONOR": ORANGE}

S = pd.read_csv(a.sites, sep="\t")
S["hom"] = (S.af >= 0.9) & (S.dp >= 10)
S = S[(S.groupby(["chr", "pos"]).hom.transform("sum") - S.hom.astype(int)) >= 1].copy()
S["alt"] = (np.rint(S.af * S.dp) >= 3) & (S.af >= 0.2)
T = pd.read_csv(a.tracts, sep="\t")


def regions(x):
    """targeted regions (SNPs within 200 kb) and their call"""
    x = x.sort_values("pos")
    rid = (x.pos.diff() > 2e5).cumsum()
    out = []
    for _, g in x.groupby(rid):
        al = g[g.alt]
        if len(al) >= 2:
            m = al.af.median()
            call = "BB" if m >= 0.8 else ("AB" if 0.25 <= m <= 0.75 else "AB")
        else:
            call = "none"
        out.append((g.pos.min(), g.pos.max(), call, len(g)))
    return out


fig, axes = plt.subplots(len(ex), 1, figsize=(12.5, 2.75 * len(ex)))
for k, (ax, (smp, ch, s0, e0, title)) in enumerate(zip(np.atleast_1d(axes), ex)):
    s0, e0 = float(s0), float(e0)
    view = None
    if "@" in ch:                                   # chr1@280-309: zoom to 280-309 Mb
        ch, v = ch.split("@"); view = tuple(float(q) for q in v.split("-"))
    x = S[(S["sample"] == smp) & (S.chr == ch)]
    pool = x.pool.iloc[0]
    t = T[(T.SAMPLE == pool) & (T.CHROM == ch)]
    L = t.END.max() / 1e6
    Lw = (view[1] - view[0]) if view else L          # scale for minimum widths
    Y1, Y2, H = 1.62, 1.30, 0.20                  # row 1 (k-mer) and row 2 (targseq call) bands
    for r in t.itertuples():
        ax.add_patch(Rectangle((r.START / 1e6, Y1), (r.END - r.START) / 1e6, H, color=COL[r.STATE], lw=0))
    ax.add_patch(Rectangle((0, Y2 + H / 2 - 0.01), L, 0.02, color="#eeece6", lw=0))
    for p0, p1, call, n in sorted(regions(x), key=lambda r: r[2] != "none"):   # informative regions drawn last
        w = max((p1 - p0) / 1e6, 0.008 * Lw)
        x0 = (p0 + p1) / 2e6 - w / 2
        if call == "none":       # uninformative or AA: a short light tick
            ax.add_patch(Rectangle((x0, Y2 + 0.05), w, H - 0.10, color="#dddbd3", lw=0))
        else:
            ax.add_patch(Rectangle((x0, Y2), w, H, color={"AB": BLUE, "BB": ORANGE}[call], lw=0))
    ax.add_patch(Rectangle((s0, Y2 - 0.04), e0 - s0, Y1 + H - Y2 + 0.08, fill=False, ec=INK, lw=1.5))
    ax.axvspan(s0, e0, ymin=0, ymax=(1.08 + 0.06) / 1.95, color=BLUE, alpha=0.07, lw=0)
    ax.axhline(1.18, color="#e4e2da", lw=0.8)
    # row 3: SNPs
    ref = x[~x.alt]; het = x[x.alt & (x.af < 0.8)]; hom = x[x.alt & (x.af >= 0.8)]
    ax.scatter(ref.pos / 1e6, ref.af, s=6, color=PALE, alpha=0.6, lw=0)
    ax.scatter(het.pos / 1e6, het.af, s=15, color=BLUE, lw=0)
    ax.scatter(hom.pos / 1e6, hom.af, s=15, color=ORANGE, lw=0)
    for yy in (0.5, 1.0):
        ax.axhline(yy, color="#eeece6", lw=0.8, zorder=0)
    for yy, lab in [(Y1 + H / 2, "k-mer call (pool)"), (Y2 + H / 2, "targeted seq. call"), (0.5, "targeted seq. SNPs:\ndonor allele fraction")]:
        ax.text(-0.035, yy, lab, ha="right", va="center", fontsize=8, color=INK, transform=ax.get_yaxis_transform())
    ax.set_title(f"{k + 1}. {smp}, {ch}: {title}", loc="left", fontsize=9, color=INK)
    ax.set_xlim(*(view or (0, L))); ax.set_ylim(-0.06, 1.89)
    ax.set_yticks([0, 0.5, 1]); ax.set_yticklabels(["0", "0.5", "1"])
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(labelsize=8, colors=MUTED)
axes[-1].set_xlabel("position (Mb, B73 v5)", fontsize=9, color=MUTED)
h = [Rectangle((0, 0), 1, 1, color=c) for c in (GREY, BLUE, ORANGE, "#dddbd3")]
fig.legend(h, ["k-mer B73", "k-mer SEG  /  targeted seq. AB (one teosinte copy)",
               "k-mer DONOR  /  targeted seq. BB (two teosinte copies)", "targeted seq.: no donor allele seen (AA, or the donor shares the B73 allele)"],
           ncol=2, frameon=False, fontsize=8.5, loc="upper center", bbox_to_anchor=(0.5, 1.0))
fig.tight_layout(rect=(0.09, 0, 1, 0.965))
fig.savefig(a.out, dpi=200)
print("wrote", a.out)
