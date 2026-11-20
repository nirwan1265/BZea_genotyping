#!/usr/bin/env python3
"""
new_model/fix1 / 3 — choose, per taxon and window, between the old counts (all B73 single-copy k-mers) and the
new counts (informative k-mers only), and set the DONOR level and its spread for the windows that switch.

For every pool of the taxon and every B73 control, the score r of each window is computed twice with step 20's
normalisation (pool / median, non-negative least-squares fit of the controls, refit without low windows,
r / median): old (step 17 windows, step 19 denominators) and new (fix1/2 windows, fix1/1 denominators).
Then per window, with the CURRENT calls as labels (as in the pilot):
  noise    = robust SD of r - 1 over the controls (each scored against the other controls) + pools called B73 there
  level    = median r over pools called DONOR there (tracts >= 2 Mb); needs >= --min_pools such pools
  contrast = (1 - level) x g(1/2) / noise      (g = step 20's detection curve at the taxon's median coverage)
A window switches to the new counts when the taxon's informative k-mers in B73 >= --min_new, both contrasts can
be measured, and contrast_new >= --gain x contrast_old.
For the switched windows the assemblies cannot give the DONOR level (by construction they have none of these
k-mers), so it comes from our own pools:
  level_hmm = level / median per-pool scale of the taxon (the HMM multiplies the level by each pool's fitted scale)
  tau_hmm   = between-pool spread of the DONOR level = sqrt(max(0, SD(r of DONOR pools)^2 - noise^2)), smoothed over
              +-5 windows, / median scale (the HMM's --d_spread term for these windows)
The labels come from the current calls, so the level of a switched window is learned from pools called DONOR there
(86-90% right in the targeted sequencing); the SEG islands under test are not used.

Output (--out_dir): <TAXON>.npz (use_new, level_hmm, tau_hmm, contrast_old, contrast_new; read by fx1_kmer_hmm.py),
choice_<TAXON>.tsv.gz (per window), summary_<TAXON>.tsv (per chromosome and hotspot).
Needs: numpy, scipy, pandas.
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
ap.add_argument("--taxon", required=True)
ap.add_argument("--fx_dir", required=True, help="${FX}: sets/, windows/<T>/")
ap.add_argument("--d16", required=True, help="step 17 windows: <SAMPLE>.tsv.gz")
ap.add_argument("--d19", required=True, help="step 19: B73.tsv.gz")
ap.add_argument("--pools_tsv", required=True, help="step 11 pools.tsv (SAMPLE ROLE GROUP)")
ap.add_argument("--hmm_dir", required=True, help="current step 20 .npz (lam, scale, sigma per pool)")
ap.add_argument("--rho_file", required=True, help="current step 21 rho_loglik.tsv")
ap.add_argument("--tracts", required=True, help="current tracts_all.tsv.gz")
ap.add_argument("--out_dir", required=True)
ap.add_argument("--min_old", type=float, default=2000)
ap.add_argument("--min_new", type=float, default=200)
ap.add_argument("--min_pools", type=int, default=5)
ap.add_argument("--gain", type=float, default=1.2)
ap.add_argument("--noisy_sigma", type=float, default=0.15)
ap.add_argument("--smooth", type=int, default=5)
a = ap.parse_args()
warnings.simplefilter("ignore", RuntimeWarning)
os.makedirs(a.out_dir, exist_ok=True)
T = a.taxon
HOT = [("chr6", 21, 23), ("chr1", 296, 298), ("chr7", 54, 65), ("chr3", 75, 78)]


def read_win(path):
    d = pd.read_csv(path, sep="\t", usecols=["seqname", "start", "observed_kmers"])
    return d.seqname.values, d.start.values.astype(np.int64), d.observed_kmers.values.astype(float)


CH, ST, Bold = read_win(f"{a.d19}/B73.tsv.gz")
W = len(Bold)
idx_of = {(c, s): i for i, (c, s) in enumerate(zip(CH, ST))}


def aligned(path):
    c, s, o = read_win(path)
    if len(o) == W and (s == ST).all() and (c == CH).all():
        return o
    out = np.zeros(W)
    for ci, si, oi in zip(c, s, o):
        j = idx_of.get((ci, si))
        if j is not None: out[j] = oi
    return out


Bnew = aligned(f"{a.fx_dir}/sets/B73_inf_{T}.tsv.gz")
NEWD = f"{a.fx_dir}/windows/{T}"
pt = pd.read_csv(a.pools_tsv, sep="\t", usecols=[0, 1, 2]); pt.columns = ["SAMPLE", "ROLE", "GROUP"]
ctrl = [s for s in pt[pt.ROLE == "control"].SAMPLE if os.path.exists(f"{a.d16}/{s}.tsv.gz") and os.path.exists(f"{NEWD}/{s}.tsv.gz")]
if len(ctrl) < 5:
    raise SystemExit(f"{T}: only {len(ctrl)} controls with new windows in {NEWD}; need >= 5")
donors = [s for s in pt[(pt.ROLE == "donor") & (pt.GROUP == T)].SAMPLE]
missing = [s for s in donors if not os.path.exists(f"{NEWD}/{s}.tsv.gz")]
donors = [s for s in donors if s not in missing and os.path.exists(f"{a.d16}/{s}.tsv.gz")]
print(f"{T}: {len(donors)} pools with new windows, {len(missing)} without (they keep the old counts in the HMM), {len(ctrl)} controls")
if missing: print("  without new windows: " + ",".join(missing[:30]) + (" ..." if len(missing) > 30 else ""))

with open(a.rho_file) as f:
    cols = f.readline().rstrip("\n").split("\t"); rows = [l.rstrip("\n").split("\t") for l in f]
rho = float(max(rows, key=lambda r: float(r[cols.index("ALL") if "ALL" in cols else 1]))[cols.index("RHO")])
info = {}
for s in donors:
    p = f"{a.hmm_dir}/{s}.npz"
    if os.path.exists(p):
        with np.load(p, allow_pickle=True) as z:
            j = int(np.argmin(np.abs(z["rho_grid"] - rho)))
            info[s] = (float(z["lam"]), float(z["par_best"][j]), float(z["sigma"]))
use_pools = [s for s in donors if s in info and info[s][2] <= a.noisy_sigma]
lam = float(np.median([info[s][0] for s in use_pools])) if use_pools else 0.8
scale = float(np.median([info[s][1] for s in use_pools])) if use_pools else 1.0
g = (1 - math.exp(-0.5 * lam)) / (1 - math.exp(-lam))
print(f"{T}: {len(use_pools)} pools (sigma <= {a.noisy_sigma}) used to measure levels and noise; best rho {rho:g}, "
      f"median coverage lam {lam:.3f} (g(1/2) = {g:.3f}), median scale {scale:.2f}")


def scores(obs, ctrl_obs, B, min_k, by_share, cut=2.0):
    """step 20's normalisation. Windows dropped where the controls see little: step 20 (old counts) compares the
    controls' count with the genome median (Rm > 0.2); for the informative k-mers, whose number per window varies
    much more (few where the teosinte keeps most B73 k-mers), the SHARE of the window's k-mers the controls see is
    compared with the typical share instead (by_share)."""
    N = obs / np.median(obs[obs > 0])
    X = np.vstack([c / np.median(c[c > 0]) for c in ctrl_obs]).T
    if by_share:
        Q = np.vstack([c / np.maximum(B, 1) for c in ctrl_obs]).T
        Rm = np.median(Q / np.median(Q[B >= min_k], axis=0)[None, :], axis=1)
    else:
        Rm = np.median(X, axis=1)
    ok = (B >= min_k) & (Rm > 0.2)
    keep = ok.copy(); r = np.full(len(obs), np.nan)
    for _ in range(8):
        b, _ = nnls(X[keep], N[keep]); fit = X @ b
        r = N / np.where(fit > 0, fit, np.nan)
        sig = 1.4826 * np.nanmedian(np.abs(r[keep] - 1))
        new = ok & (r > 1 - cut * sig)
        if (new == keep).all(): break
        keep = new
    r = r / np.nanmedian(r[keep])
    return np.where(ok, r, np.nan).astype(np.float32)


# labels from the current calls: 0 B73, 1 SEG, 2 DONOR (tract >= 2 Mb), -1 other
TR = pd.read_csv(a.tracts, sep="\t")
TR = TR[TR.SAMPLE.isin(use_pools)]
mid = ST + 50000
chr_idx = {c: np.where(CH == c)[0] for c in np.unique(CH)}
lab = {s: np.full(W, -1, np.int8) for s in use_pools}
for r in TR.itertuples():
    ii = chr_idx.get(r.CHROM)
    if ii is None: continue
    m = mid[ii]; lo, hi = np.searchsorted(m, r.START), np.searchsorted(m, r.END)
    code = {"B73": 0, "SEG": 1, "DONOR": 2}.get(r.STATE, -1)
    if code == 2 and (r.END - r.START) < 2e6: code = -1
    lab[r.SAMPLE][ii[lo:hi]] = code

out = {}
for kind, B, dirn, mink in (("old", Bold, a.d16, a.min_old), ("new", Bnew, NEWD, a.min_new)):
    C = [aligned(f"{dirn}/{c}.tsv.gz") for c in ctrl]
    bs = kind == "new"
    CR = np.vstack([scores(C[i], [c for j, c in enumerate(C) if j != i], B, mink, bs) for i in range(len(C))])
    R = np.vstack([scores(aligned(f"{dirn}/{s}.tsv.gz"), C, B, mink, bs) for s in use_pools])
    L = np.vstack([lab[s] for s in use_pools])
    B73r = np.vstack([np.where(L == 0, R, np.nan), CR])
    nb = np.isfinite(B73r).sum(0)
    noise = np.where(nb >= 5, 1.4826 * np.nanmedian(np.abs(B73r - 1), axis=0), np.nan)
    D = np.where(L == 2, R, np.nan); nd = np.isfinite(D).sum(0)
    level = np.where(nd >= a.min_pools, np.nanmedian(D, axis=0), np.nan)
    sdD = np.where(nd >= a.min_pools, 1.4826 * np.nanmedian(np.abs(D - level[None, :]), axis=0), np.nan)
    contrast = (1 - level) * g / noise
    out[kind] = dict(noise=noise, level=level, sdD=sdD, contrast=contrast, nd=nd)
    print(f"{T} {kind}: windows with a measured contrast {int(np.isfinite(contrast).sum())} of {W}; "
          f"median contrast {np.nanmedian(contrast):.2f}")

co, cn = out["old"]["contrast"], out["new"]["contrast"]
use = (Bnew >= a.min_new) & np.isfinite(co) & np.isfinite(cn) & (out["new"]["level"] < 0.95) & (cn >= a.gain * co)
lev_hmm = np.where(use, np.clip(out["new"]["level"], 0.01, 0.95) / scale, np.nan)
tau2 = np.where(use, np.maximum(0.0, out["new"]["sdD"] ** 2 - out["new"]["noise"] ** 2), np.nan)
tau = np.full(W, np.nan); k = a.smooth
for i in np.where(use)[0]:
    w = tau2[max(0, i - k): i + k + 1]; tau[i] = math.sqrt(np.nanmean(w)) / scale
np.savez_compressed(f"{a.out_dir}/{T}.npz", use_new=use, level_hmm=lev_hmm.astype(np.float32), tau_hmm=tau.astype(np.float32),
                    contrast_old=co.astype(np.float32), contrast_new=cn.astype(np.float32), chrom=CH, start=ST)
pd.DataFrame({"chr": CH, "start": ST, "B73_kmers_old": Bold.astype(int), "B73_kmers_new": Bnew.astype(int),
              "donor_pools": out["new"]["nd"], "level_old": np.round(out["old"]["level"], 3), "level_new": np.round(out["new"]["level"], 3),
              "noise_old": np.round(out["old"]["noise"], 3), "noise_new": np.round(out["new"]["noise"], 3),
              "contrast_old": np.round(co, 2), "contrast_new": np.round(cn, 2), "use_new": use.astype(int),
              "level_hmm": np.round(lev_hmm, 3), "tau_hmm": np.round(tau, 3)}
             ).to_csv(f"{a.out_dir}/choice_{T}.tsv.gz", sep="\t", index=False, compression="gzip")

rows = []
q = lambda x: f"{np.nanmedian(x):.2f}" if np.isfinite(x).any() else "-"
regions = [(c, None, None) for c in sorted(np.unique(CH), key=lambda c: int(c[3:]) if c[3:].isdigit() else 99)] + HOT + [("genome", None, None)]
for c, s0, e0 in regions:
    m = np.ones(W, bool) if c == "genome" else (CH == c)
    if s0 is not None: m &= (ST >= s0 * 1e6) & (ST < (e0 + 1) * 1e6)
    meas = m & np.isfinite(co)
    best = np.where(use, cn, co)
    rows.append({"taxon": T, "region": c if s0 is None else f"{c} {s0}-{e0} Mb", "windows": int(m.sum()),
                 "switched": int((m & use).sum()), "switched_%": round(100 * (m & use).sum() / max(1, meas.sum()), 1),
                 "contrast_old": q(co[meas]), "contrast_after": q(best[meas]),
                 "low_contrast_<3_old_%": round(100 * np.mean(co[meas] < 3), 1) if meas.any() else np.nan,
                 "low_contrast_<3_after_%": round(100 * np.mean(best[meas] < 3), 1) if meas.any() else np.nan})
S = pd.DataFrame(rows)
S.to_csv(f"{a.out_dir}/summary_{T}.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 20)
print(S.to_string(index=False))
print("wrote", a.out_dir)
