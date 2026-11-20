#!/usr/bin/env python3
"""
pilot fix 1 / step 3 — does counting only informative k-mers raise the DONOR vs SEG contrast?

For the pilot pools on chr3 + chr7, the score r of each window is computed twice with the SAME procedure as step 20
(pool / median, non-negative least-squares fit of the B73 controls, refit without low windows, r / median):
  old: all B73 single-copy k-mers (step 17 windows, step 19 denominators)
  new: only the informative k-mers of the pool's taxon (p1 steps 1-2: B73 k-mers absent from all the taxon's assemblies)
Then, per taxon and window, with the CURRENT calls as labels (B73 tracts are right 97%, DONOR tracts 86-90%):
  noise  = robust SD of r - 1 over the B73 controls (each scored against the other controls) and the pools called
           B73 there (the pilot pools carry teosinte at the hotspots, so the controls carry the noise estimate there)
  level  = median r over pools called DONOR there (tracts >= 2 Mb; not the SEG islands under test)
  gap    = (1 - level) x g(1/2)   (expected SEG score above DONOR, as in step 20)
  contrast = gap / noise
Reported for the two low-contrast hotspots (chr7 54-65, chr3 75-78 Mb) and for the rest of chr3/chr7, old vs new.
Direct test on the errors: in the SEG islands between two DONOR tracts at the hotspots, is each window's r closer to
the DONOR level or to the SEG level, old vs new? (If the islands are DONOR, the new counts should move them to DONOR.)
Figure: r along chr3/chr7 for a few island pools, rescaled so that DONOR level = 0 and B73 = 1, old vs new.

Usage (from new_model/pilot_fix1/, env with numpy, scipy, pandas, matplotlib), see README.md.
"""
import argparse
import glob
import gzip
import math
import os
import warnings
import numpy as np
import pandas as pd
from scipy.optimize import nnls

ap = argparse.ArgumentParser()
ap.add_argument("--pilot_dir", required=True, help="${KMER_DIR}/new_model/pilot_fix1 (sets/, windows/)")
ap.add_argument("--pools", default="pools.tsv")
ap.add_argument("--d16", required=True, help="step 17 windows (old k-mers): <SAMPLE>.tsv.gz")
ap.add_argument("--d19", required=True, help="step 19: B73.tsv.gz (old denominators)")
ap.add_argument("--hmm_dir", required=True, help="step 20 .npz of the current run (lam per pool)")
ap.add_argument("--tracts", required=True, help="current tracts_all.tsv.gz")
ap.add_argument("--out", required=True)
ap.add_argument("--min_old", type=float, default=2000)
ap.add_argument("--min_new", type=float, default=200)
ap.add_argument("--rule", default="all", choices=["all", "any"], help="which informative set (p1 RULE)")
ap.add_argument("--taxa", default="Zv,Zx")
ap.add_argument("--targseq", default="../../data/introgression_spread/targseq/tracts_targseq.tsv", help="34d tracts (optional)")
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
TAG = "_any" if a.rule == "any" else ""
WIN = f"{a.pilot_dir}/windows{TAG}"
warnings.simplefilter("ignore", RuntimeWarning)
CHRS = ("chr3", "chr7")
HOT = {"chr7 54-65 Mb": ("chr7", 54e6, 66e6), "chr3 75-78 Mb": ("chr3", 75e6, 79e6)}


def read_win(path):
    d = pd.read_csv(path, sep="\t", usecols=["seqname", "start", "observed_kmers"])
    d = d[d.seqname.isin(CHRS)]
    return d.set_index(["seqname", "start"]).observed_kmers


P = pd.read_csv(a.pools, sep="\t")
ctrl = sorted({os.path.basename(f)[:-7] for f in glob.glob(f"{WIN}/{a.taxa.split(',')[0]}/*.tsv.gz")} - set(P.SAMPLE))
print(f"pilot pools {len(P)}, controls with new windows {len(ctrl)}")
Bold = read_win(f"{a.d19}/B73.tsv.gz")
IDX = Bold.index
T = pd.read_csv(a.tracts, sep="\t"); T = T[T.CHROM.isin(CHRS)].sort_values(["SAMPLE", "CHROM", "START"])
T["prev"] = T.groupby(["SAMPLE", "CHROM"]).STATE.shift(1); T["next"] = T.groupby(["SAMPLE", "CHROM"]).STATE.shift(-1)
T["island"] = (T.STATE == "SEG") & (T.prev == "DONOR") & (T.next == "DONOR")
pos = np.array([s for _, s in IDX]) + 50000; chrom = np.array([c for c, _ in IDX])


def labels(sample):
    """per window: state code (0 B73, 1 SEG, 2 DONOR, -1 none), tract length (Mb), island flag"""
    st = np.full(len(IDX), -1); ln = np.zeros(len(IDX)); isl = np.zeros(len(IDX), bool)
    for r in T[T.SAMPLE == sample].itertuples():
        m = (chrom == r.CHROM) & (pos >= r.START) & (pos < r.END)
        st[m] = {"B73": 0, "SEG": 1, "DONOR": 2}[r.STATE]; ln[m] = (r.END - r.START) / 1e6; isl[m] = r.island
    return st, ln, isl


def scores(obs, ctrl_obs, B, min_k, cut=2.0):
    """step 20's normalisation on these windows"""
    N = obs / np.median(obs[obs > 0])
    X = np.vstack([c / np.median(c[c > 0]) for c in ctrl_obs]).T
    Rm = np.median(X, axis=1)
    ok = (B >= min_k) & (Rm > 0.2)
    keep = ok.copy()
    for _ in range(8):
        b, _ = nnls(X[keep], N[keep]); fit = X @ b
        r = N / np.where(fit > 0, fit, np.nan)
        sig = 1.4826 * np.nanmedian(np.abs(r[keep] - 1))
        new = ok & (r > 1 - cut * sig)
        if (new == keep).all(): break
        keep = new
    r = r / np.nanmedian(r[keep])
    return np.where(ok, r, np.nan)


lam = {}
for s in P.SAMPLE:
    f = f"{a.hmm_dir}/{s}.npz"
    if os.path.exists(f):
        with np.load(f, allow_pickle=True) as z: lam[s] = float(z["lam"])

res = {}
for tx in a.taxa.split(","):
    pools = [s for s in P[P.TAXON == tx].SAMPLE if os.path.exists(f"{WIN}/{tx}/{s}.tsv.gz") and os.path.exists(f"{a.d16}/{s}.tsv.gz")]
    if not pools: print(f"{tx}: no pools with both old and new windows"); continue
    Bnew = read_win(f"{a.pilot_dir}/sets/B73_inf_{tx}{TAG}.tsv.gz").reindex(IDX).fillna(0).values
    Co = [read_win(f"{a.d16}/{c}.tsv.gz").reindex(IDX).fillna(0).values for c in ctrl]
    Cn = [read_win(f"{WIN}/{tx}/{c}.tsv.gz").reindex(IDX).fillna(0).values for c in ctrl]
    Ro, Rn, ST, LN, ISL = [], [], [], [], []
    for s in pools:
        Ro.append(scores(read_win(f"{a.d16}/{s}.tsv.gz").reindex(IDX).fillna(0).values, Co, Bold.values, a.min_old))
        Rn.append(scores(read_win(f"{WIN}/{tx}/{s}.tsv.gz").reindex(IDX).fillna(0).values, Cn, Bnew, a.min_new))
        st, ln, isl = labels(s); ST.append(st); LN.append(ln); ISL.append(isl)
    Ro, Rn, ST, LN, ISL = map(np.vstack, (Ro, Rn, ST, LN, ISL))
    # controls scored against the other controls: B73 everywhere, for the noise
    CRo = np.vstack([scores(Co[i], [c for j, c in enumerate(Co) if j != i], Bold.values, a.min_old) for i in range(len(Co))])
    CRn = np.vstack([scores(Cn[i], [c for j, c in enumerate(Cn) if j != i], Bnew, a.min_new) for i in range(len(Cn))])
    g = (1 - math.exp(-0.5 * np.median([lam.get(s, 0.8) for s in pools]))) / (1 - math.exp(-np.median([lam.get(s, 0.8) for s in pools])))
    out = {}
    for lab, R, CR in (("old", Ro, CRo), ("new", Rn, CRn)):
        B73r = np.vstack([np.where(ST == 0, R, np.nan), CR])
        nb = np.sum(np.isfinite(B73r), axis=0); nd = np.sum((ST == 2) & (LN >= 2) & np.isfinite(R), axis=0)
        noise = 1.4826 * np.nanmedian(np.abs(B73r - 1), axis=0)
        level = np.nanmedian(np.where((ST == 2) & (LN >= 2), R, np.nan), axis=0)
        noise = np.where(nb >= 5, noise, np.nan); level = np.where(nd >= 2, level, np.nan)
        gap = (1 - level) * g; contrast = gap / noise
        # islands: closer to DONOR or to SEG?
        segl = level + gap
        dD = np.abs(R - level[None, :]); dS = np.abs(R - segl[None, :])
        out[lab] = dict(R=R, noise=noise, level=level, contrast=contrast, closer_donor=np.where(ISL & np.isfinite(dD) & np.isfinite(dS), dD < dS, np.nan),
                        nb=nb, nd=nd)
    res[tx] = dict(pools=pools, out=out, ST=ST, ISL=ISL, g=g)

rows = []
for tx, d in res.items():
    for reg, (c, s0, e0) in list(HOT.items()) + [("rest of chr3 + chr7", (None, 0, 0))]:
        m = (chrom == c) & (pos >= s0) & (pos < e0) if c else ~np.any([(chrom == cc) & (pos >= ss) & (pos < ee) for cc, ss, ee in HOT.values()], axis=0)
        r = {"taxon": tx, "region": reg, "pools": len(d["pools"])}
        for lab in ("old", "new"):
            o = d["out"][lab]
            r[f"contrast_{lab}"] = round(float(np.nanmedian(o["contrast"][m])), 2)
            r[f"noise_{lab}"] = round(float(np.nanmedian(o["noise"][m])), 3)
            r[f"donor_level_{lab}"] = round(float(np.nanmedian(o["level"][m])), 3)
            cd = o["closer_donor"][:, m]
            r[f"island_windows_closer_to_DONOR_{lab}"] = (f"{100 * np.nanmean(cd):.0f}% of {int(np.isfinite(cd).sum())}" if np.isfinite(cd).any() else "-")
        r["windows_with_estimate"] = int(np.isfinite(d["out"]["old"]["contrast"][m] + d["out"]["new"]["contrast"][m]).sum())
        rows.append(r)
S = pd.DataFrame(rows)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
print("\nContrast = (1 - DONOR level) x g(1/2) / noise, medians over windows; old = all B73 k-mers, new = informative only")
print(S.to_string(index=False))
S.to_csv(f"{a.out}/contrast_old_new.tsv", sep="\t", index=False)

# per pool at each hotspot: median score on a 0 = DONOR level, 1 = B73 scale (SEG expected at g(1/2) ~ 0.6),
# over the windows the pool calls SEG-island there (or all its windows there), with the targeted-seq genotype if any
TS = None
if a.targseq and os.path.exists(a.targseq):
    TS = pd.read_csv(a.targseq, sep="\t")
prow = []
for tx, d in res.items():
    for reg, (c, s0, e0) in HOT.items():
        m = (chrom == c) & (pos >= s0) & (pos < e0)
        for i, s in enumerate(d["pools"]):
            isl = d["ISL"][i] & m
            w = isl if isl.any() else m
            st = d["ST"][i][w]; st = st[st >= 0]
            call = "SEG island" if isl.any() else ({0: "B73", 1: "SEG", 2: "DONOR"}[int(np.bincount(st).argmax())] if len(st) else "-")
            r = {"taxon": tx, "region": reg, "pool": s, "call_there": call}
            for lab in ("old", "new"):
                o = d["out"][lab]; y = (o["R"][i] - o["level"]) / (1 - o["level"])
                r[f"score_{lab}"] = round(float(np.nanmedian(y[w])), 2) if np.isfinite(y[w]).any() else np.nan
            if TS is not None:
                t = TS[(TS.pool == s) & (TS.chr == c) & (TS.start < e0) & (TS.end > s0)]
                r["targseq"] = ",".join(f"{g}({n})" for g, n in zip(t.plant_genotype, t.sites)) if len(t) else ""
            prow.append(r)
PS = pd.DataFrame(prow)
PS.to_csv(f"{a.out}/pool_scores.tsv", sep="\t", index=False)
print(f"\nPer pool (score: 0 = DONOR level, ~{np.median([d['g'] for d in res.values()]):.2f} = SEG, 1 = B73; "
      "targseq = plant genotype(sites) of the targeted-seq tracts overlapping the hotspot)")
for (tx, reg), g in PS.groupby(["taxon", "region"], sort=False):
    print(f"\n{tx} {reg}")
    print(g[g.call_there != "B73"].drop(columns=["taxon", "region"]).to_string(index=False))

# figure: island pools, r rescaled to DONOR level = 0, B73 = 1
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
ex = []
for tx, d in res.items():
    for i, s in enumerate(d["pools"]):
        if d["ISL"][i].any(): ex.append((tx, i, s))
ex = ex[:4]
if ex:
    fig, axes = plt.subplots(len(ex), 1, figsize=(12, 2.4 * len(ex)), squeeze=False)
    for ax, (tx, i, s) in zip(axes[:, 0], ex):
        d = res[tx]; c = chrom == ("chr7" if d["ISL"][i][chrom == "chr7"].any() else "chr3")
        for lab, col in (("old", "#a8a69e"), ("new", "#2a78d6")):
            o = d["out"][lab]; y = (o["R"][i] - o["level"]) / (1 - o["level"])
            ax.plot(pos[c] / 1e6, y[c], lw=0.8, color=col, label=f"{lab} k-mers")
        isl = d["ISL"][i] & c
        if isl.any(): ax.axvspan(pos[isl].min() / 1e6, pos[isl].max() / 1e6, color="#eb6834", alpha=0.12, lw=0, label="SEG island (current call)")
        ax.axhline(0, color="#eb6834", lw=0.8, ls="--"); ax.axhline(d["g"], color="#2a78d6", lw=0.8, ls=":"); ax.axhline(1, color="#52514e", lw=0.8, ls="--")
        ax.set_ylim(-0.6, 1.6); ax.set_ylabel("0 = DONOR\n~0.6 = SEG\n1 = B73", fontsize=7)
        ax.set_title(f"{s} ({tx}), {'chr7' if c[chrom == 'chr7'].any() else 'chr3'}", loc="left", fontsize=9)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0, 0].legend(fontsize=8, frameon=False, ncol=3, loc="upper right")
    axes[-1, 0].set_xlabel("position (Mb)")
    fig.tight_layout(); fig.savefig(f"{a.out}/island_examples.png", dpi=150)
print("wrote", a.out)
