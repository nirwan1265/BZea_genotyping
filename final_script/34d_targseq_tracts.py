#!/usr/bin/env python3
"""
34d — Tract-level check of the k-mer calls with the targeted sequencing (targseq), one plant per plot.

Reads 34c's sites_used.tsv.gz (sample, pool, chr, pos, dp, af, state) and the k-mer tracts. Compares whole tracts,
not 1-Mb bins, so that each comparison rests on many SNPs.

Donor-allele evidence at a site ("alt site"): at least --min_alt alt reads and AF >= --min_af. Single noisy reads
(about 2% of sites in our B73 regions have AF 0.1-0.35) are not counted.

If the targseq DNA is one plant of the plot, the plant's genotype in a tract of ours should be:
  our B73   -> AA (no alt sites)
  our DONOR -> BB (alt sites at AF ~1)
  our SEG   -> AA, AB or BB: 1/4 : 1/2 : 1/4 under Mendel (the plant is a BC2S3 child of an AB mother)
Plant genotype of a tract (tracts >= --min_mb, with >= --min_sites covered sites):
  BB if >= 2 alt sites and their median AF >= 0.8; AB if >= 2 alt sites and median AF 0.25-0.75;
  AA if no more alt sites than noise allows (the sample's own noise rate in our B73 tracts, Poisson 99th
  percentile) and at least --min_sites covered sites; otherwise unclear.
  A DONOR region whose teosinte allele equals B73 at almost every panel SNP can look AA; --min_sites guards
  against that (see the AA share of our DONOR tracts: it is the false-AA rate).
Hidden introgressions in our B73 tracts: 2-Mb windows (1-Mb steps) with more alt sites than the sample's noise
rate allows (Poisson, p < 0.001); reported per sample against the number expected by chance.

Outputs (--out): tracts_targseq.tsv (one row per sample x tract), summary.tsv (plant genotype by our call,
by number and by Mb), hidden_b73_windows.tsv. Printed: the same summaries.
Usage: python final_script/34d_targseq_tracts.py --sites $D34/compare/sites_used.tsv.gz --tracts $D21/tracts_all.tsv.gz --out $D34/tracts
Needs: numpy, pandas, scipy.
"""
import argparse
import os
import numpy as np
import pandas as pd
from scipy.stats import poisson

ap = argparse.ArgumentParser()
ap.add_argument("--sites", required=True)
ap.add_argument("--tracts", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--min_alt", type=int, default=3)
ap.add_argument("--min_af", type=float, default=0.2)
ap.add_argument("--min_mb", type=float, default=1.0)
ap.add_argument("--min_sites", type=int, default=30)
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
ST = {0: "B73", 1: "SEG", 2: "DONOR"}

S = pd.read_csv(a.sites, sep="\t")
S["alt"] = np.rint(S.af * S.dp).astype(int)
S["is_alt"] = (S.alt >= a.min_alt) & (S.af >= a.min_af)
T = pd.read_csv(a.tracts, sep="\t")
T = T[T.SAMPLE.isin(S.pool.unique())].copy()
T["Mb"] = (T.END - T.START) / 1e6
print(f"{S['sample'].nunique()} samples, {len(S):,} sample x site; alt sites: >= {a.min_alt} alt reads and AF >= {a.min_af}")

# noise: alt-site rate in our B73 tracts, per sample (fallback: overall)
b = S[S.state == 0]
noise = b.groupby("sample").is_alt.mean()
noise_all = b.is_alt.mean()
print(f"noise (alt sites in our B73 regions): median {noise.median():.4f} per site (overall {noise_all:.4f})")

# ---- tract level
rows = []
for (smp, pool), g in S.groupby(["sample", "pool"]):
    nr = max(noise.get(smp, noise_all), 1e-4)
    tt = T[(T.SAMPLE == pool) & (T.Mb >= a.min_mb)]
    for ch, gc in g.groupby("chr"):
        pos = gc.pos.values; alt = gc.is_alt.values; af = gc.af.values
        for r in tt[tt.CHROM == ch].itertuples():
            m = (pos >= r.START) & (pos < r.END)
            n, k = int(m.sum()), int(alt[m].sum())
            med = float(np.median(af[m][alt[m]])) if k else np.nan
            if k >= 2 and med >= 0.8: geno = "BB"
            elif k >= 2 and 0.25 <= med <= 0.75: geno = "AB"
            elif n >= a.min_sites and k <= poisson.ppf(0.99, nr * n): geno = "AA"
            else: geno = "unclear"
            rows.append({"sample": smp, "pool": pool, "chr": ch, "start": r.START, "end": r.END, "Mb": round(r.Mb, 2),
                         "our_call": r.STATE, "sites": n, "alt_sites": k, "median_alt_af": round(med, 3) if k else np.nan,
                         "plant_genotype": geno})
R = pd.DataFrame(rows)
R.to_csv(os.path.join(a.out, "tracts_targseq.tsv"), sep="\t", index=False)
G = ["AA", "AB", "BB", "unclear"]
sm = []
for call in ["B73", "SEG", "DONOR"]:
    x = R[R.our_call == call]
    n = x.plant_genotype.value_counts().reindex(G, fill_value=0)
    mb = x.groupby("plant_genotype").Mb.sum().reindex(G, fill_value=0)
    clear = n[["AA", "AB", "BB"]].sum(); clear_mb = mb[["AA", "AB", "BB"]].sum()
    sm.append({"our_call": call, "tracts": len(x), **{f"n_{g}": int(n[g]) for g in G},
               **{f"share_{g}": round(n[g] / clear, 3) if clear else np.nan for g in ["AA", "AB", "BB"]},
               **{f"Mb_share_{g}": round(mb[g] / clear_mb, 3) if clear_mb else np.nan for g in ["AA", "AB", "BB"]}})
SM = pd.DataFrame(sm)
SM.to_csv(os.path.join(a.out, "summary.tsv"), sep="\t", index=False)
print(f"\nPlant genotype in our tracts (>= {a.min_mb:g} Mb; shares among clear tracts). Mendel for SEG: AA 0.25, AB 0.50, BB 0.25")
print(SM.to_string(index=False))

# ---- hidden introgressions in our B73 tracts: 2-Mb windows with too many alt sites
hid = []
for (smp, pool), g in b.groupby(["sample", "pool"]):
    nr = max(noise.get(smp, noise_all), 1e-4)
    for ch, gc in g.groupby("chr"):
        pos = gc.pos.values; alt = gc.is_alt.values
        if len(pos) == 0: continue
        for w0 in np.arange(0, pos.max() + 1, 1e6):
            m = (pos >= w0) & (pos < w0 + 2e6)
            n, k = int(m.sum()), int(alt[m].sum())
            if n < 10: continue
            p = poisson.sf(k - 1, nr * n)
            hid.append({"sample": smp, "chr": ch, "start_Mb": w0 / 1e6, "sites": n, "alt_sites": k, "p": p,
                        "median_alt_af": float(np.median(gc.af.values[m][alt[m]])) if k else np.nan})
H = pd.DataFrame(hid)
H["flag"] = H.p < 1e-3
H.to_csv(os.path.join(a.out, "hidden_b73_windows.tsv"), sep="\t", index=False)
print(f"\nOur B73 regions: {len(H):,} 2-Mb windows tested; {int(H.flag.sum())} with more alt sites than noise allows "
      f"(p < 0.001; about {len(H) * 1e-3:.0f} expected by chance), in {H[H.flag]['sample'].nunique()} samples")
if H.flag.any():
    f = H[H.flag]
    print("  median AF of their alt sites: " + ", ".join(f"{q:.2f}" for q in f.median_alt_af.quantile([0.25, 0.5, 0.75])) + " (25/50/75%)")
print("wrote", a.out)
