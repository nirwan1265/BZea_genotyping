#!/usr/bin/env python3
"""
new_model / contrast_check — can the HMM tell DONOR from SEG in each window?

In a window, the expected score of a DONOR pool is the level d_w and that of a SEG pool is about
d_w + (1 - d_w) g(1/2), g = the k-mer detection curve (step 20). The gap between them is (1 - d_w) g(1/2).
Contrast = gap / noise of the window (noise = median pool sigma x the window's noise factor, measured as in
20c from the spread of all pools' scores near 1). Where the teosinte keeps most of the B73 k-mers (d_w high),
the gap is small, DONOR and SEG overlap, and a DONOR stretch can be read as SEG whatever level is used.

Reports, per taxon: contrast quantiles over the genome and in the four places where SEG tracts between two DONOR
tracts recur (chr6 21-23, chr1 296-298, chr7 54-65, chr3 75-78 Mb), and the contrast of the windows where such
tracts occur in >= 5 pools, against all windows (with 300-570 pools nearly every window has at least one).
Output: contrast_<TAXON>.tsv.gz (window, chr, start, level, noise, contrast, pools with a sandwiched SEG there).
Usage (from new_model/, any env with numpy + pandas):
  D19=...; D20=...; D21=...   (config.sh of the current run)
  python contrast_check.py --hmm_dir $D20 --rho_file $D21/rho_loglik.tsv --ref_dir $D19 --assemblies "$ASMS" \
     --tracts $D21/tracts_all.tsv.gz --out_dir $DNEW/contrast
Needs: numpy, pandas.
"""
import argparse
import collections
import glob
import gzip
import math
import os
import warnings
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--hmm_dir", required=True); ap.add_argument("--rho_file", required=True)
ap.add_argument("--ref_dir", required=True); ap.add_argument("--assemblies", required=True)
ap.add_argument("--tracts", required=True); ap.add_argument("--out_dir", required=True)
ap.add_argument("--noisy_sigma", type=float, default=0.15)
a = ap.parse_args()
os.makedirs(a.out_dir, exist_ok=True)
HOT = [("chr6", 21, 23), ("chr1", 296, 298), ("chr7", 54, 65), ("chr3", 75, 78)]


def read(path, col):
    with gzip.open(path, "rt") as f:
        h = f.readline().rstrip("\n").split("\t"); i = h.index(col)
        return np.array([l.rstrip("\n").split("\t")[i] for l in f])


CH = read(f"{a.ref_dir}/B73.tsv.gz", "seqname"); ST = read(f"{a.ref_dir}/B73.tsv.gz", "start").astype(int)
B = read(f"{a.ref_dir}/B73.tsv.gz", "observed_kmers").astype(float)
lev = collections.defaultdict(list)
for item in a.assemblies.split(","):
    n, t = item.split(":")[:2]
    lev[t].append(read(f"{a.ref_dir}/{n}.tsv.gz", "observed_kmers").astype(float) / np.maximum(B, 1))
lev = {t: np.mean(v, axis=0) for t, v in lev.items()}

with open(a.rho_file) as f:
    cols = f.readline().rstrip("\n").split("\t"); rows = [l.rstrip("\n").split("\t") for l in f]
rho = float(max(rows, key=lambda r: float(r[cols.index("ALL") if "ALL" in cols else 1]))[cols.index("RHO")])
P = collections.defaultdict(list)
for f in sorted(glob.glob(f"{a.hmm_dir}/*.npz")):
    with np.load(f, allow_pickle=True) as z:
        if str(z["role"]) != "donor" or float(z["sigma"]) > a.noisy_sigma: continue
        t = str(z["taxon"])
        if t not in lev: continue
        j = int(np.argmin(np.abs(z["rho_grid"] - rho)))
        P[t].append((z["r"].astype(np.float32), z["ok"].astype(bool), float(z["sigma"]), float(z["lam"]), float(z["par_best"][j])))

T = pd.read_csv(a.tracts, sep="\t").sort_values(["SAMPLE", "CHROM", "START"])
T["prev"] = T.groupby(["SAMPLE", "CHROM"]).STATE.shift(1); T["next"] = T.groupby(["SAMPLE", "CHROM"]).STATE.shift(-1)
SW = T[(T.STATE == "SEG") & (T.prev == "DONOR") & (T.next == "DONOR")]

rows = []
for t, pl in sorted(P.items()):
    R = np.vstack([p[0] for p in pl]).astype(float); OK = np.vstack([p[1] for p in pl]) & np.isfinite(R)
    SIG = np.array([p[2] for p in pl])[:, None]
    lam = float(np.median([p[3] for p in pl])); scale = float(np.median([p[4] for p in pl]))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        Z = np.where(OK & (np.abs(R - 1) < 4 * SIG), (R - 1) / SIG, np.nan)
        wf = np.clip(1.4826 * np.nanmedian(np.abs(Z), axis=0), 0.5, 3)
    wf = np.where(np.isfinite(wf), wf, 1.0)
    noise = float(np.median(SIG)) * wf
    g_half = (1 - math.exp(-0.5 * lam)) / (1 - math.exp(-lam))
    d = np.clip(scale * lev[t], 0.01, 0.95)
    contrast = (1 - d) * g_half / noise
    # pools of this taxon with a sandwiched SEG tract covering each window
    sw = SW[SW.GROUP == t]; nsw = np.zeros(len(B), int)
    for r in sw.itertuples():
        nsw[(CH == r.CHROM) & (ST + 50000 >= r.START) & (ST + 50000 < r.END)] += 1
    ok = B >= 2000
    df = pd.DataFrame({"chr": CH, "start": ST, "level": d.round(4), "noise": noise.round(4), "contrast": contrast.round(2), "sandwiched_pools": nsw})
    df.to_csv(f"{a.out_dir}/contrast_{t}.tsv.gz", sep="\t", index=False, compression="gzip")
    q = lambda x: "/".join(f"{v:.1f}" for v in np.quantile(x, [0.1, 0.5, 0.9])) if len(x) else "-"
    row = {"taxon": t, "pools": len(pl), "genome (10/50/90%)": q(contrast[ok]),
           "windows with sandwiched SEG in >= 5 pools": q(contrast[ok & (nsw >= 5)]),
           "low contrast (<3) genome": f"{100 * np.mean(contrast[ok] < 3):.1f}%",
           "low contrast (<3) where sandwiched in >= 5 pools": f"{100 * np.mean(contrast[ok & (nsw >= 5)] < 3):.1f}%" if (ok & (nsw >= 5)).any() else "-"}
    for c, s, e in HOT:
        m = ok & (CH == c) & (ST >= s * 1e6) & (ST < (e + 1) * 1e6)
        row[f"{c} {s}-{e}"] = f"{np.median(contrast[m]):.1f} (level {np.median(d[m]):.2f})" if m.any() else "-"
    rows.append(row)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 20)
O = pd.DataFrame(rows)
print(f"Contrast = gap between SEG and DONOR expected scores / noise (rho {rho:g}); medians unless stated")
print(O.to_string(index=False))
O.to_csv(f"{a.out_dir}/contrast_summary.tsv", sep="\t", index=False)
print("wrote", a.out_dir)
