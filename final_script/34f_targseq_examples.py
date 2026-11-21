#!/usr/bin/env python3
"""
34f — Examples: our k-mer calls against the targeted sequencing, one chromosome per panel (README section 2.9.1).

Each panel has three rows, top to bottom. "Ours" = our k-mer call; "theirs" = the targeted sequencing (one plant of the plot?).
  1. "OUR call":     our HMM call of the 6-plant pool from 0.8x skim-seq: B73 / SEG / DONOR
  2. "THEIR call":   their call per targeted region (SNPs within 200 kb of each other): AB if its donor-carrying SNPs have
                     a median fraction below 0.8, BB if >= 0.8; a short light tick = no donor allele seen (fewer than 2 SNPs
                     carry donor alleles: AA, or the donor carries the B73 allele at these SNPs, so the region says nothing)
  3. "THEIR SNPs":   the donor allele fraction at each real SNP (homozygous alt in another sample, as in 34d). Dots come
                     in clumps: one clump is usually one amplicon (10-20 SNPs in ~150 bp, from the same reads)
Above row 1 (with --targseq_tracts), for each of our SEG/DONOR tracts that 34d calls clearly: "ours SEG -> theirs AB",
i.e. our call and the plant's genotype over the whole tract (34d, from all its SNPs). The verdict is per tract, not per dot.
Same colours: grey = B73 / AA, blue = SEG / AB, orange = DONOR / BB.
The outlined tract is the one described in the title.

Usage (the README figure):
  python final_script/34f_targseq_examples.py --sites data/introgression_fix1/targseq/sites_used.tsv.gz \
     --tracts data/introgression_fix1/tracts_all.tsv.gz --targseq_tracts data/introgression_fix1/targseq/tracts_targseq.tsv \
     --out figures/fix1/targseq_examples.png   [--example SAMPLE:CHR[@FROM-TO]:START_MB:END_MB:"title" ...]
     (@FROM-TO: zoom the panel to FROM-TO Mb; without --example the four README examples below)
  All samples, one figure each (every chromosome with an introgression): the same with --gallery_dir figures/fix1/targseq_all
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
    "Zv0370_P3_P4_P1111:chr5:28.99:116.57:ours SEG, theirs AB (agree)",
    "Zv0420_P3_P1_P1111:chr5:40.79:97.37:ours DONOR, theirs BB (agree)",
    "Zv0210_P1_P5_P1111:chr7:54.18:66.18:the outlined stretch was SEG in the old model; ours is now DONOR, theirs BB (fixed)",
    "Zv0030_P1_P3_P1111:chr1@280-309:297.01:299.21:zoomed to 280–309 Mb: ours SEG between two DONOR tracts, theirs BB (still our error)",
]
ap = argparse.ArgumentParser()
ap.add_argument("--sites", required=True)
ap.add_argument("--tracts", required=True)
ap.add_argument("--out", default="")
ap.add_argument("--gallery_dir", default="", help="one figure per sample (all its introgressed chromosomes) instead of the examples")
ap.add_argument("--targseq_tracts", default="", help="34d tracts_targseq.tsv: label each tract with its targeted-seq genotype (gallery)")
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


def draw(ax, smp, ch, s0, e0, title, lab=None):
    """one panel: k-mer call, targeted-seq call, targeted-seq SNPs. s0/e0: outlined tract (None = none);
    lab: (start, end, our call, plant genotype) of the tracts the targeted seq. calls clearly; printed above the k-mer
    row as "ours SEG -> theirs AB" so that our call and theirs are never mixed up"""
    view = None
    if "@" in ch:                                   # chr1@280-309: zoom to 280-309 Mb
        ch, v = ch.split("@"); view = tuple(float(q) for q in v.split("-"))
    x = S[(S["sample"] == smp) & (S.chr == ch)]
    pool = x.pool.iloc[0]
    t = T[(T.SAMPLE == pool) & (T.CHROM == ch)]
    L = t.END.max() / 1e6 if len(t) else x.pos.max() / 1e6
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
    if s0 is not None:
        ax.add_patch(Rectangle((s0, Y2 - 0.04), e0 - s0, Y1 + H - Y2 + 0.08, fill=False, ec=INK, lw=1.5))
        ax.axvspan(s0, e0, ymin=0, ymax=(1.08 + 0.06) / 1.95, color=BLUE, alpha=0.07, lw=0)
    for a0, a1, ours, g in (lab or []):
        w = (a1 - a0) / 1e6
        txt = f"ours {ours} → theirs {g}" if w >= 0.13 * Lw else (f"{ours}→{g}" if w >= 0.045 * Lw else "")
        if txt:
            ax.text((a0 + a1) / 2e6, Y1 + H + 0.07, txt, ha="center", va="bottom", fontsize=7, color=INK)
    ax.axhline(1.18, color="#e4e2da", lw=0.8)
    # row 3: SNPs
    ref = x[~x.alt]; het = x[x.alt & (x.af < 0.8)]; hom = x[x.alt & (x.af >= 0.8)]
    ax.scatter(ref.pos / 1e6, ref.af, s=6, color=PALE, alpha=0.6, lw=0)
    ax.scatter(het.pos / 1e6, het.af, s=15, color=BLUE, lw=0)
    ax.scatter(hom.pos / 1e6, hom.af, s=15, color=ORANGE, lw=0)
    for yy in (0.5, 1.0):
        ax.axhline(yy, color="#eeece6", lw=0.8, zorder=0)
    for yy, lb in [(Y1 + H / 2, "OUR call\n(k-mer, pool of 6 plants)"), (Y2 + H / 2, "THEIR call\n(targeted seq., per region)"),
                    (0.5, "THEIR SNPs:\ndonor allele fraction")]:
        ax.text(-0.035, yy, lb, ha="right", va="center", fontsize=8, color=INK, transform=ax.get_yaxis_transform())
    ax.set_title(title, loc="left", fontsize=9, color=INK)
    ax.set_xlim(*(view or (0, L))); ax.set_ylim(-0.06, 2.04 if lab else 1.89)
    ax.set_yticks([0, 0.5, 1]); ax.set_yticklabels(["0", "0.5", "1"])
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(labelsize=8, colors=MUTED)


def legend(fig):
    h = [Rectangle((0, 0), 1, 1, color=c) for c in (GREY, BLUE, ORANGE, "#dddbd3")]
    fig.legend(h, ["no teosinte: ours B73 (all plants AA)", "one teosinte copy: ours SEG (plot segregates)  /  theirs AB",
                   "two teosinte copies: ours DONOR (all plants BB)  /  theirs BB",
                   "theirs: no donor allele seen (AA, or the donor shares the B73 allele)"],
               ncol=2, frameon=False, fontsize=8.5, loc="upper center", bbox_to_anchor=(0.5, 1.0))


TT = pd.read_csv(a.targseq_tracts, sep="\t") if a.targseq_tracts else None


def labs(smp, ch):
    """our call and the plant's genotype (34d) for each of our SEG/DONOR tracts that the targeted seq. calls clearly"""
    if TT is None:
        return None
    q = TT[(TT["sample"] == smp) & (TT.chr == ch.split("@")[0]) & (TT.our_call != "B73") & TT.plant_genotype.isin(["AA", "AB", "BB"])]
    return [(r.start, r.end, r.our_call, r.plant_genotype) for r in q.itertuples()]


if a.gallery_dir:
    # one figure per sample: every chromosome with a k-mer SEG/DONOR call or a targeted-seq AB/BB region;
    # above the k-mer row, the genotype of each clearly called tract (34d: AA / AB / BB)
    import os
    os.makedirs(a.gallery_dir, exist_ok=True)
    # the three bad samples of README 2.9.1 point 6: flagged in the title so the figure is not read as a model error
    BAD = {"Zv0480_P3_P2_P1211": "BAD SAMPLE: contaminated with another sample's DNA (README 2.9.1, point 6)",
           "Zx0090_P2_P5_P1111": "BAD SAMPLE: swap or wrong plant (README 2.9.1, point 6)",
           "Zx0160_P2_P4_P1111": "BAD SAMPLE: swap or wrong plant (README 2.9.1, point 6)"}
    made = []
    for smp, xs in S.groupby("sample"):
        pool = xs.pool.iloc[0]
        chs = []
        for ch in sorted(xs.chr.unique(), key=lambda c: int(c[3:])):
            t = T[(T.SAMPLE == pool) & (T.CHROM == ch)]
            if (t.STATE != "B73").any() or any(r[2] != "none" for r in regions(xs[xs.chr == ch])):
                chs.append(ch)
        if not chs or not (T.SAMPLE == pool).any():
            continue
        fig, axes = plt.subplots(len(chs), 1, figsize=(12.5, 2.6 * len(chs) + 0.5), squeeze=False)
        for ax, ch in zip(axes[:, 0], chs):
            lab = labs(smp, ch)
            draw(ax, smp, ch, None, None, f"{smp} (pool {pool}), {ch}" + (f"  —  {BAD[smp]}" if smp in BAD else ""), lab)
        axes[-1, 0].set_xlabel("position (Mb, B73 v5)", fontsize=9, color=MUTED)
        legend(fig)
        fig.tight_layout(rect=(0.12, 0, 1, 1 - 0.45 / (2.6 * len(chs) + 0.5)))
        out = f"{a.gallery_dir}/{smp}.png"; fig.savefig(out, dpi=110); plt.close(fig)
        made.append((smp, pool, len(chs)))
    with open(f"{a.gallery_dir}/index.tsv", "w") as f:
        f.write("sample\tpool\tchromosomes\n"); f.writelines(f"{s}\t{p}\t{n}\n" for s, p, n in made)
    print(f"wrote {len(made)} sample figures to {a.gallery_dir}")
else:
    fig, axes = plt.subplots(len(ex), 1, figsize=(12.5, 2.75 * len(ex)))
    for k, (ax, (smp, ch, s0, e0, title)) in enumerate(zip(np.atleast_1d(axes), ex)):
        draw(ax, smp, ch, float(s0), float(e0), f"{k + 1}. {smp}, {ch.split('@')[0]}: {title}", labs(smp, ch))
    np.atleast_1d(axes)[-1].set_xlabel("position (Mb, B73 v5)", fontsize=9, color=MUTED)
    legend(fig)
    fig.tight_layout(rect=(0.12, 0, 1, 0.965))
    fig.savefig(a.out, dpi=200)
    print("wrote", a.out)
