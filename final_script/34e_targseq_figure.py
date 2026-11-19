#!/usr/bin/env python3
"""
34e — Figure for the targeted-sequencing validation (README section 2.8.1, Check 4).

Panels:
  a. Donor-allele sites per 100 SNPs inside our B73, SEG and DONOR calls: all sites in the targseq VCF, and real SNPs only
     (homozygous alt in at least one other sample; step 34d). Shows how the amplification errors are removed.
  b. Plant genotype (AA / AB / BB) of our tracts in the targseq data, share of the clearly called tracts (step 34d).
  c. Median donor allele fraction of each tract that carries donor alleles: our SEG and DONOR tracts, with the
     values expected for one AB plant (0.5) and from the k-mer estimate inside SEG pools (0.46, step 22b).
Usage: python final_script/34e_targseq_figure.py --sites data/introgression_spread/targseq/sites_used.tsv.gz \
         --tracts data/introgression_spread/targseq/tracts_targseq.tsv --out figures/targseq_validation.png
Needs: numpy, pandas, matplotlib.
"""
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--sites", required=True)
ap.add_argument("--tracts", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--min_alt", type=int, default=3)
ap.add_argument("--min_af", type=float, default=0.2)
a = ap.parse_args()
CALLS = ["B73", "SEG", "DONOR"]
GREY, BLUE, ORANGE, INK, MUTED = "#a8a69e", "#2a78d6", "#eb6834", "#0b0b0b", "#52514e"

S = pd.read_csv(a.sites, sep="\t")
S["is_alt"] = (np.rint(S.af * S.dp) >= a.min_alt) & (S.af >= a.min_af)
S["hom"] = (S.af >= 0.9) & (S.dp >= 10)
S["real"] = (S.groupby(["chr", "pos"]).hom.transform("sum") - S.hom.astype(int)) >= 1
rate_all = [100 * S.is_alt[S.state == i].mean() for i in range(3)]
rate_real = [100 * S.is_alt[(S.state == i) & S.real].mean() for i in range(3)]
R = pd.read_csv(a.tracts, sep="\t")

fig, ax = plt.subplots(1, 3, figsize=(14, 4.2), gridspec_kw={"width_ratios": [1, 1.1, 1.3]})
for x in ax:
    x.spines[["top", "right"]].set_visible(False)
    x.tick_params(colors=MUTED, labelsize=9)

# a. alt-site rate, all sites vs real SNPs
xx = np.arange(3); w = 0.38
b1 = ax[0].bar(xx - w / 2 - 0.01, rate_all, w, color=GREY, label=f"all {S[['chr', 'pos']].drop_duplicates().shape[0]:,} sites")
b2 = ax[0].bar(xx + w / 2 + 0.01, rate_real, w, color=BLUE, label=f"real SNPs ({S[S.real][['chr', 'pos']].drop_duplicates().shape[0]:,})")
for bars in (b1, b2):
    for r in bars:
        ax[0].text(r.get_x() + r.get_width() / 2, r.get_height() + 0.6, f"{r.get_height():.1f}", ha="center", fontsize=8, color=INK)
ax[0].set_xticks(xx); ax[0].set_xticklabels([f"our {c}" for c in CALLS], fontsize=9)
ax[0].set_ylabel("sites with donor alleles per 100 SNPs", fontsize=9, color=MUTED)
ax[0].legend(frameon=False, fontsize=8, loc="upper left")
ax[0].set_title("a. Donor alleles by our call", loc="left", fontsize=10, color=INK)

# b. plant genotype of our tracts
cols = {"AA": GREY, "AB": BLUE, "BB": ORANGE}
for i, c in enumerate(CALLS):
    x = R[(R.our_call == c) & R.plant_genotype.isin(["AA", "AB", "BB"])]
    sh = x.plant_genotype.value_counts(normalize=True).reindex(["AA", "AB", "BB"], fill_value=0)
    left = 0
    for g in ["AA", "AB", "BB"]:
        ax[1].barh(i, sh[g], left=left, color=cols[g], edgecolor="white", linewidth=2, height=0.62)
        if sh[g] >= 0.06:
            ax[1].text(left + sh[g] / 2, i, f"{g}\n{100 * sh[g]:.0f}%", ha="center", va="center", fontsize=8,
                       color="white" if g != "AA" else INK)
        left += sh[g]
    ax[1].text(1.02, i, f"{len(x):,} tracts", va="center", fontsize=8, color=MUTED)
ax[1].set_yticks(range(3)); ax[1].set_yticklabels([f"our {c}" for c in CALLS], fontsize=9)
ax[1].invert_yaxis(); ax[1].set_xlim(0, 1.18); ax[1].set_xticks([0, 0.25, 0.5, 0.75, 1])
ax[1].set_xticklabels(["0", "25%", "50%", "75%", "100%"])
ax[1].set_xlabel("share of clearly called tracts (≥ 1 Mb)", fontsize=9, color=MUTED)
ax[1].set_title("b. Genotype in the targeted sequencing", loc="left", fontsize=10, color=INK)

# c. median donor allele fraction of carrier tracts
bins = np.linspace(0.2, 1.0, 17)
for c, col in [("SEG", BLUE), ("DONOR", ORANGE)]:
    x = R[(R.our_call == c) & R.plant_genotype.isin(["AB", "BB"])].median_alt_af.dropna()
    ax[2].hist(x, bins=bins, color=col, alpha=0.75, label=f"our {c} ({len(x)} tracts, median {x.median():.2f})", edgecolor="white")
ax[2].axvline(0.46, color=MUTED, lw=1.2, ls=":", label="0.46: k-mer estimate inside SEG pools")
ax[2].axvline(0.5, color=INK, lw=1, ls="--", label="0.5: one AB plant")
ax[2].set_xlabel("median donor allele fraction of the tract", fontsize=9, color=MUTED)
ax[2].set_ylabel("tracts", fontsize=9, color=MUTED)
ax[2].legend(frameon=True, facecolor="white", edgecolor="none", framealpha=1, fontsize=8, loc="upper left")
ax[2].set_title("c. Donor allele fraction in carrier tracts", loc="left", fontsize=10, color=INK)

fig.tight_layout()
fig.savefig(a.out, dpi=200)
print("wrote", a.out, "| rates all", np.round(rate_all, 2), "real", np.round(rate_real, 2))
