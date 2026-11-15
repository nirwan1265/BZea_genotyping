#!/usr/bin/env python3
"""
24 — Result figures: observed vs expected BC2S3 proportions (PNG, for the README and slides).

Reads steps 21-22 and writes to <out_dir>:
  fig1_overall.png        all non-noisy donor pools: observed (bar, 95% family-bootstrap CI) vs expected (line)
  fig2_by_taxon.png       the same per taxon, one panel per measure
  fig3_by_family.png      per family: donor alleles and AB/(AB+BB), expected as a line
  fig4_by_chromosome.png  AB/(AB+BB) per chromosome with CI, expected as a line
  fig5_pool_share.png     distribution of donor share across pools, expected as a line

Pool states -> plant genotypes: AA = B73 + SEG/4, AB = SEG/2, BB = DONOR + SEG/4, donor alleles = BB + AB/2.
Expected BC2S3: SEG 6.25%, DONOR 9.375%, AB 3.125%, BB 10.9375%, donor alleles 12.5%, AB/(AB+BB) 22.2%.
Only non-noisy donor pools (FLAG ok). Family = first two '_' fields of the pedigree (as in steps 20b/22).
"""
import argparse
import collections
import csv
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

EXP = {"SEG": 6.25, "DONOR": 9.375, "AB": 3.125, "BB": 10.9375, "DONOR_ALLELES": 12.5, "AB_OF_DONOR_PLANTS": 100 * 0.03125 / 0.140625}
LAB = {"SEG": "SEG\n(pool)", "DONOR": "DONOR\n(pool)", "AB": "AB\n(plant)", "BB": "BB\n(plant)", "DONOR_ALLELES": "Donor\nalleles",
       "AB_OF_DONOR_PLANTS": "AB/(AB+BB)"}
TAXA = ["Zd", "Zh", "Zl", "Zv", "Zx"]
NAME = {"Zd": "Zd (diploperennis)", "Zh": "Zh (huehuetenangensis)", "Zl": "Zl (luxurians)", "Zv": "Zv (parviglumis)", "Zx": "Zx (mexicana)"}
COL = dict(zip(TAXA, ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]))   # validated categorical slots 1-5
BAR, INK, INK2, GRID, EXPC = "#2a78d6", "#0b0b0b", "#52514e", "#e6e5e1", "#0b0b0b"


def style(ax):
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"): ax.spines[sp].set_color("#b9b8b3"); ax.spines[sp].set_linewidth(0.8)
    ax.tick_params(colors=INK2, labelsize=9, length=3, width=0.8)
    ax.yaxis.grid(True, color=GRID, linewidth=0.8); ax.set_axisbelow(True)


def title(ax, t, sub):
    ax.set_title(t, loc="left", pad=24)
    ax.text(0, 1.015, sub, transform=ax.transAxes, fontsize=8.5, color=INK2, va="bottom")


def exp_line(ax, x0, x1, y, label=None, side="right"):
    ax.plot([x0, x1], [y, y], color=EXPC, lw=1.6, solid_capstyle="butt", zorder=4)
    if label: ax.annotate(label, (x1 if side == "right" else x0, y), xytext=(4 if side == "right" else -4, 0), textcoords="offset points",
                          ha="left" if side == "right" else "right", va="center", fontsize=8.5, color=INK)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stats_dir", required=True, help="step 21/22 output dir (pool_summary.tsv, stats_overall.tsv, stats_by_chrom.tsv)")
    ap.add_argument("--meta", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--label", default="", help="added to the figure titles, e.g. 'current method'")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "figure.dpi": 100, "savefig.dpi": 200, "axes.titlesize": 11,
                         "axes.titleweight": "bold", "axes.titlecolor": INK, "axes.labelcolor": INK2, "axes.labelsize": 9.5})
    sub = f" ({a.label})" if a.label else ""

    st = collections.defaultdict(dict)
    with open(f"{a.stats_dir}/stats_overall.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            st[r["SCOPE"]][r["MEASURE"]] = (100 * float(r["OBSERVED"]), 100 * float(r["CI_LOW"]), 100 * float(r["CI_HIGH"]),
                                            float(r["P"]), int(r["FAMILIES"]), int(r["POOLS"]))

    # ---- fig 1: overall
    ms = ["SEG", "AB", "DONOR", "BB", "DONOR_ALLELES"]
    fig, ax = plt.subplots(figsize=(7.5, 4.2)); style(ax)
    for i, m in enumerate(ms):
        o, lo, hi, p, nf, npool = st["ALL"][m]
        ax.bar(i, o, width=0.5, color=BAR, zorder=3)
        ax.errorbar(i, o, yerr=[[o - lo], [hi - o]], fmt="none", ecolor=INK2, elinewidth=1.2, capsize=4, zorder=5)
        exp_line(ax, i - 0.36, i + 0.36, EXP[m])
        ax.annotate(f"{o:.2f}%", (i, max(hi, EXP[m] if abs(EXP[m] - hi) < 0.6 else hi)), xytext=(0, 5), textcoords="offset points",
                    ha="center", fontsize=8.5, color=INK)
    ax.set_xticks(range(len(ms))); ax.set_xticklabels([f"{LAB[m]}\nexp. {EXP[m] + 1e-9:.2f}%" for m in ms], fontsize=9, color=INK)
    ax.set_ylabel("% of genome"); ax.set_ylim(0, 14.5)
    nf, npool = st["ALL"]["SEG"][4:6]
    title(ax, f"Observed vs expected BC2S3, all taxa{sub}",
          f"{npool:,} non-noisy pools, {nf} families. Bars: observed (95% family-bootstrap CI). Black line: expected.")
    fig.tight_layout(); fig.savefig(f"{a.out_dir}/fig1_overall.png"); plt.close(fig)

    # ---- fig 2: by taxon, one panel per measure
    ms2 = ["DONOR_ALLELES", "BB", "AB", "AB_OF_DONOR_PLANTS"]
    taxa = [t for t in TAXA if t in st]
    fig, axs = plt.subplots(1, 4, figsize=(12, 3.8))
    for ax, m in zip(axs, ms2):
        style(ax)
        for i, t in enumerate(taxa + ["ALL"]):
            o, lo, hi, *_ = st[t][m]
            ax.bar(i, o, width=0.62, color=COL.get(t, "#52514e"), zorder=3)
            ax.errorbar(i, o, yerr=[[o - lo], [hi - o]], fmt="none", ecolor=INK2, elinewidth=1, capsize=3, zorder=5)
        exp_line(ax, -0.5, len(taxa) + 0.5, EXP[m])
        ax.set_xticks(range(len(taxa) + 1)); ax.set_xticklabels(taxa + ["All"], fontsize=9, color=INK)
        ax.set_title(f"{LAB[m].replace(chr(10), ' ')}  ·  expected {EXP[m]:.1f}%", loc="left", fontsize=10)
        top = max(EXP[m], max(st[t][m][2] for t in taxa + ["ALL"])); ax.set_ylim(0, top * 1.18)
    for ax in axs[:3]: ax.set_ylabel("% of genome")
    axs[3].set_ylabel("% of donor-carrying plants")
    fig.suptitle(f"Observed vs expected by taxon{sub}  ·  bars: observed with 95% CI, black line: expected BC2S3",
                 x=0.01, ha="left", fontsize=11, fontweight="bold", color=INK)
    fig.tight_layout(); fig.savefig(f"{a.out_dir}/fig2_by_taxon.png"); plt.close(fig)

    # ---- per family (from pool_summary)
    fam = {}
    with open(a.meta, encoding="utf-8", errors="ignore") as f:
        next(f)
        for line in f:
            p = line.rstrip("\r\n").split("\t")
            if len(p) > 1: fam[p[0]] = "_".join(p[1].split("_")[:2])
    by = collections.defaultdict(list); share = []
    with open(f"{a.stats_dir}/pool_summary.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["ROLE"] != "donor" or r["FLAG"] != "ok": continue
            s, d = float(r["FRAC_SEG"]), float(r["FRAC_DONOR"])
            by[(r["GROUP"], fam.get(r["SAMPLE"], "unknown"))].append((s, d)); share.append((r["GROUP"], 100 * (d + s / 2)))
    rows = []
    for (t, fm), v in by.items():
        s, d = np.mean([x[0] for x in v]), np.mean([x[1] for x in v]); ab, bb = s / 2, d + s / 4
        rows.append((t, fm, len(v), 100 * (bb + ab / 2), 100 * ab / (ab + bb) if ab + bb > 0 else np.nan))
    rows.sort(key=lambda x: (TAXA.index(x[0]) if x[0] in TAXA else 9, -x[3]))
    with open(f"{a.out_dir}/family_values.tsv", "w") as f:
        f.write("TAXON\tFAMILY\tPOOLS\tDONOR_ALLELES_PCT\tAB_OF_DONOR_PLANTS_PCT\n")
        for r in rows: f.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]:.3f}\t{r[4]:.2f}\n")

    # ---- fig 3: by family
    x = np.arange(len(rows)); cols = [COL.get(r[0], "#52514e") for r in rows]
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(13, 6.6), sharex=True)
    for ax in (a1, a2): style(ax)
    a1.bar(x, [r[3] for r in rows], width=0.78, color=cols, zorder=3)
    exp_line(a1, -0.8, len(rows) - 0.2, EXP["DONOR_ALLELES"], "expected 12.5%")
    a1.axhline(st["ALL"]["DONOR_ALLELES"][0], color=INK2, lw=1, zorder=2)
    a1.annotate(f"all pools {st['ALL']['DONOR_ALLELES'][0]:.1f}%", (len(rows) - 0.2, st["ALL"]["DONOR_ALLELES"][0]), xytext=(4, 0),
                textcoords="offset points", va="center", fontsize=8.5, color=INK2)
    a1.set_ylabel("Donor alleles (% of genome)"); a1.set_ylim(0, 18)
    a1.set_title(f"Per family{sub}: donor alleles (top) and share of donor-carrying plants that are het (bottom)", loc="left")
    a2.bar(x, [r[4] for r in rows], width=0.78, color=cols, zorder=3)
    exp_line(a2, -0.8, len(rows) - 0.2, EXP["AB_OF_DONOR_PLANTS"], "expected 22.2%")
    a2.axhline(st["ALL"]["AB_OF_DONOR_PLANTS"][0], color=INK2, lw=1, zorder=2)
    a2.annotate(f"all pools {st['ALL']['AB_OF_DONOR_PLANTS'][0]:.1f}%", (len(rows) - 0.2, st["ALL"]["AB_OF_DONOR_PLANTS"][0]), xytext=(4, 0),
                textcoords="offset points", va="center", fontsize=8.5, color=INK2)
    a2.set_ylabel("AB/(AB+BB) (%)"); a2.set_ylim(0, 75)
    a2.set_xticks([]); a2.set_xlim(-1, len(rows))
    for t in TAXA:                                  # taxon group labels under the bars
        idx = [i for i, r in enumerate(rows) if r[0] == t]
        if not idx: continue
        n = sum(rows[i][2] for i in idx)
        a2.annotate(t, ((idx[0] + idx[-1]) / 2, 0), xytext=(0, -6), textcoords="offset points",
                    ha="center", va="top", fontsize=8.5, color=INK, annotation_clip=False)
        a2.plot([idx[0] - 0.4, idx[-1] + 0.4], [-1.5, -1.5], color=COL[t], lw=3, clip_on=False, solid_capstyle="butt")
    present = [t for t in TAXA if any(r[0] == t for r in rows)]
    hs = [plt.Rectangle((0, 0), 1, 1, color=COL[t]) for t in present]
    labs = [f"{NAME[t]}: {sum(1 for r in rows if r[0] == t)} fam., {sum(r[2] for r in rows if r[0] == t)} pools" for t in present]
    a1.legend(hs, labs, ncol=3, fontsize=8.5, frameon=False,
              loc="upper left", bbox_to_anchor=(0, 1.0))
    fig.text(0.01, 0.005, "Each bar is one family (mean of its non-noisy pools), sorted by donor alleles within taxon. "
             "Grey line: mean over all pools. Black line: expected BC2S3.", fontsize=8.5, color=INK2)
    fig.tight_layout(rect=(0, 0.06, 1, 1)); fig.savefig(f"{a.out_dir}/fig3_by_family.png"); plt.close(fig)

    # ---- fig 4: by chromosome
    ch = collections.defaultdict(dict)
    with open(f"{a.stats_dir}/stats_by_chrom.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            ch[r["CHROM"]][r["MEASURE"]] = (100 * float(r["OBSERVED"]), 100 * float(r["CI_LOW"]), 100 * float(r["CI_HIGH"]))
    cs = sorted(ch, key=lambda c: int(c.replace("chr", "")) if c.replace("chr", "").isdigit() else 99)
    fig, ax = plt.subplots(figsize=(8.5, 3.8)); style(ax)
    for i, c in enumerate(cs):
        o, lo, hi = ch[c]["AB_OF_DONOR_PLANTS"]
        ax.bar(i, o, width=0.6, color=BAR, zorder=3)
        ax.errorbar(i, o, yerr=[[o - lo], [hi - o]], fmt="none", ecolor=INK2, elinewidth=1, capsize=3, zorder=5)
    exp_line(ax, -0.5, len(cs) - 0.5, EXP["AB_OF_DONOR_PLANTS"], "expected 22.2%")
    ax.set_xticks(range(len(cs))); ax.set_xticklabels(cs, fontsize=9, color=INK); ax.set_xlim(-0.6, len(cs) + 0.9)
    ax.set_ylabel("AB/(AB+BB) (%)"); ax.set_ylim(0, 45)
    title(ax, f"Share of donor-carrying plants that are het, per chromosome{sub}",
          "Bars: observed (95% family-bootstrap CI). Black line: expected BC2S3.")
    fig.tight_layout(); fig.savefig(f"{a.out_dir}/fig4_by_chromosome.png"); plt.close(fig)

    # ---- fig 5: pool donor share distribution
    v = np.array([s for _, s in share])
    fig, ax = plt.subplots(figsize=(8.5, 3.8)); style(ax)
    ax.hist(v, bins=np.arange(0, max(32, v.max() + 1), 1), color=BAR, rwidth=0.88, zorder=3)
    top = ax.get_ylim()[1]
    ax.plot([12.5, 12.5], [0, top], color=EXPC, lw=1.6, zorder=4)
    ax.annotate("expected 12.5%", (12.5, top), xytext=(4, -10), textcoords="offset points", fontsize=8.5, color=INK)
    med = float(np.median(v)); ax.plot([med, med], [0, top], color=INK2, lw=1, zorder=4)
    ax.annotate(f"median {med:.1f}%", (med, top), xytext=(-4, -10), textcoords="offset points", ha="right", fontsize=8.5, color=INK2)
    ax.set_xlabel("Donor share of the pool (% of genome)"); ax.set_ylabel("Pools")
    title(ax, f"Donor share per pool{sub}", f"{len(v):,} non-noisy donor pools; {np.mean(v >= 12.5) * 100:.0f}% reach the expected 12.5%.")
    fig.tight_layout(); fig.savefig(f"{a.out_dir}/fig5_pool_share.png"); plt.close(fig)
    print(f"WROTE: {a.out_dir}/fig1_overall.png fig2_by_taxon.png fig3_by_family.png fig4_by_chromosome.png fig5_pool_share.png family_values.tsv")


if __name__ == "__main__":
    main()
