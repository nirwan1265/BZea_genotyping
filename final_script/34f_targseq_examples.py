#!/usr/bin/env python3
"""
34f — Examples: our k-mer calls drawn on top of the targeted-sequencing reads, one chromosome per panel
(README section 2.9.1).

Each panel shows one sample and one chromosome:
  top band:  our k-mer call of the pool along the chromosome (B73 / SEG / DONOR); the tract discussed is outlined
  points:    donor allele fraction of the targeted sequencing at each real SNP (homozygous alt in another sample,
             as in step 34d), at its position. Grey points at 0: the plant carries B73 there (or the donor shares
             the B73 allele); blue: intermediate (0.25-0.75, one teosinte copy); orange: >= 0.8 (two copies).
The defaults are SEG tracts chosen from step 34d's table to show each outcome (no donor alleles, intermediate,
homozygous, and mixtures inside one tract), with the three bad samples of 2.9.1 left out.
The panel SNPs come in clusters (one targeted region each, often 5-30 SNPs), so points stack at one position; the
fractions of SNPs in one region can differ (probably amplification bias), so a region is read by its majority.

Usage: python final_script/34f_targseq_examples.py --sites data/introgression_spread/targseq/sites_used.tsv.gz \
         --tracts data/introgression_spread/tracts_all.tsv.gz --out figures/targseq_examples.png \
         [--example SAMPLE:CHR:START_MB:END_MB:"title" ...]
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
    "Zd0010_P3_P3_P3211:chr10:120.16:129.36:we call SEG; B73 allele only (no teosinte in this plant)",
    "Zd0010_P3_P3_P3411:chr10:81.28:152.45:we call SEG; intermediate throughout (one teosinte copy, AB)",
    "Zx0030_P2_P2_P2111:chr3:70.58:81.28:we call SEG; ≥ 0.8 (two teosinte copies, BB), like the DONOR tracts around it",
    "Zx0510_P1_P3_P2111:chr4:149.16:219.23:we call SEG; mostly intermediate, with regions ≥ 0.8 at 187–200 Mb",
    "Zx0560_P4_P3_P5111:chr2:185.74:216.04:we call SEG; B73 allele only at 186–191 Mb, intermediate at 200–201 Mb",
]
ap = argparse.ArgumentParser()
ap.add_argument("--sites", required=True)
ap.add_argument("--tracts", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--example", action="append", default=None)
a = ap.parse_args()
ex = [e.split(":", 4) for e in (a.example or EX)]
COL = {"B73": "#d3d1c7", "SEG": "#2a78d6", "DONOR": "#eb6834"}
GREY, BLUE, ORANGE, INK, MUTED = "#a8a69e", "#2a78d6", "#eb6834", "#0b0b0b", "#52514e"

S = pd.read_csv(a.sites, sep="\t")
S["hom"] = (S.af >= 0.9) & (S.dp >= 10)
S = S[(S.groupby(["chr", "pos"]).hom.transform("sum") - S.hom.astype(int)) >= 1]
S["alt"] = (np.rint(S.af * S.dp) >= 3) & (S.af >= 0.2)
T = pd.read_csv(a.tracts, sep="\t")

fig, axes = plt.subplots(len(ex), 1, figsize=(12, 2.15 * len(ex)))
for ax, (smp, ch, s0, e0, title) in zip(np.atleast_1d(axes), ex):
    s0, e0 = float(s0), float(e0)
    x = S[(S["sample"] == smp) & (S.chr == ch)]
    pool = x.pool.iloc[0]
    t = T[(T.SAMPLE == pool) & (T.CHROM == ch)]
    L = t.END.max() / 1e6
    # our calls, as a band above the points
    for r in t.itertuples():
        ax.add_patch(Rectangle((r.START / 1e6, 1.12), (r.END - r.START) / 1e6, 0.16, color=COL[r.STATE], lw=0))
    ax.add_patch(Rectangle((s0, 1.10), e0 - s0, 0.20, fill=False, ec=INK, lw=1.4))
    ax.axvspan(s0, e0, ymin=0, ymax=1.06 / 1.36, color="#2a78d6", alpha=0.06, lw=0)
    # targseq allele fractions
    ref = x[~x.alt]; het = x[x.alt & (x.af < 0.8)]; hom = x[x.alt & (x.af >= 0.8)]
    ax.scatter(ref.pos / 1e6, ref.af, s=7, color=GREY, alpha=0.6, lw=0, label="B73 allele only")
    ax.scatter(het.pos / 1e6, het.af, s=16, color=BLUE, lw=0, label="intermediate (one teosinte copy)")
    ax.scatter(hom.pos / 1e6, hom.af, s=16, color=ORANGE, lw=0, label="≥ 0.8 (two teosinte copies)")
    for y in (0.5, 1.0):
        ax.axhline(y, color="#e4e2da", lw=0.8, zorder=0)
    ins = x[(x.pos >= s0 * 1e6) & (x.pos < e0 * 1e6)]
    ax.set_title(f"{smp}, {ch}: {title}   (outlined tract {s0:.1f}–{e0:.1f} Mb: {len(ins)} SNPs, "
                 f"{int((ins.alt & (ins.af < 0.8)).sum())} intermediate, {int((ins.alt & (ins.af >= 0.8)).sum())} ≥ 0.8)",
                 loc="left", fontsize=9, color=INK)
    ax.set_xlim(0, L); ax.set_ylim(-0.06, 1.36)
    ax.set_yticks([0, 0.5, 1]); ax.set_ylabel("donor allele\nfraction", fontsize=8, color=MUTED)
    ax.text(-0.006 * L, 1.20, "our call", ha="right", va="center", fontsize=8, color=MUTED)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=8, colors=MUTED)
axes[-1].set_xlabel("position (Mb, B73 v5)", fontsize=9, color=MUTED)
h = [Rectangle((0, 0), 1, 1, color=COL[k]) for k in COL]
l1 = axes[0].legend(h, [f"our {k}" for k in COL], ncol=3, frameon=False, fontsize=8, loc="lower left", bbox_to_anchor=(0, 1.22))
axes[0].add_artist(l1)
axes[0].legend(ncol=3, frameon=False, fontsize=8, loc="lower right", bbox_to_anchor=(1, 1.22), markerscale=1.5)
fig.tight_layout()
fig.savefig(a.out, dpi=200)
print("wrote", a.out)
