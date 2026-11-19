#!/usr/bin/env python3
"""
new_model / 20c — DONOR level per window and taxon, learned from ALL pools of the taxon (no calls used).

Why: the targeted sequencing (README 2.9.1, point 7) shows that in parviglumis and mexicana, SEG tracts
between two DONOR tracts are often really DONOR, and that these errors sit at the same windows in many
pools (chr6 21-23, chr1 297, chr7 54-65, chr3 75-78 Mb). The DONOR level d_w comes from other plants'
assemblies; where our donors keep more B73 k-mers than those assemblies, DONOR reads too high and looks
like SEG. Step 20b re-estimates levels only from pools CALLED DONOR, which never sees these windows.

What this does (in the spirit of RTIGER, which fits its emission parameters across all samples at once):
in each window, the scores r of all pools of a taxon are a mixture of three groups with the same window
level d (a = the pool's fitted scale, lam = its k-mer coverage, sigma = its noise; from step 20):
    B73   : r ~ N(1, s)
    SEG   : r ~ N(a d + (1 - a d) g(f), s),  f = k/12, k ~ Binomial(12, 1/2), k = 1..11 (as in step 20)
    DONOR : r ~ N(a d, s)
    + 2% outliers (uniform 0-2)
s = pool sigma x the window's noise factor (robust spread of r among the window's pools near 1).
For each window, d is chosen on a grid by maximum likelihood, with the three mixing weights of the window
fitted by EM at each d (profile likelihood). No HMM calls are used, so the fit is not circular.
The estimate is shrunk to the assembly level by the evidence: n = expected number of introgressed
pools in the window (sum of SEG + DONOR responsibilities),  level = (n d_hat + K d_asm) / (n + K).
Windows with almost no introgressed pools keep the assembly level.

Output (--out_dir): <TAXON>.npy (float32, one level per window; step 20n reads it), levels_<TAXON>.tsv.gz
(window, assembly level, fitted level, n, log-likelihood gain), summary.tsv.
Usage: python new_model/20c_pool_levels.py --hmm_dir $D20 --rho_file $D21/rho_loglik.tsv --ref_dir $D19 \
          --assemblies "$ASMS" --out_dir $DNEW/levels [--taxa Zv,Zx] [--k 5]
Needs: numpy.
"""
import argparse
import collections
import glob
import gzip
import math
import os
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--hmm_dir", required=True, help="step 20 output (<SAMPLE>.npz) of the current run")
ap.add_argument("--rho_file", required=True, help="step 21 rho_loglik.tsv of the same run")
ap.add_argument("--ref_dir", required=True, help="step 19: B73.tsv.gz + <NAME>.tsv.gz")
ap.add_argument("--assemblies", required=True, help="NAME:TAXON,...")
ap.add_argument("--out_dir", required=True)
ap.add_argument("--taxa", default="", help="comma list (default: every taxon with an assembly)")
ap.add_argument("--plants", type=int, default=6)
ap.add_argument("--noisy_sigma", type=float, default=0.15)
ap.add_argument("--k", type=float, default=5.0, help="pseudo-pools for shrinking toward the assembly level")
ap.add_argument("--d_grid", default="0.02:0.95:0.01")
ap.add_argument("--eps", type=float, default=0.02)
ap.add_argument("--em_iter", type=int, default=15)
ap.add_argument("--chunk", type=int, default=1500)
a = ap.parse_args()
os.makedirs(a.out_dir, exist_ok=True)


def read_obs(path):
    with gzip.open(path, "rt") as f:
        head = f.readline().rstrip("\n").split("\t"); i = head.index("observed_kmers")
        return np.array([float(l.rstrip("\n").split("\t")[i]) for l in f])


with open(a.rho_file) as f:
    next(f); rows = [l.rstrip("\n").split("\t") for l in f]
cols = None
with open(a.rho_file) as f:
    cols = f.readline().rstrip("\n").split("\t")
iR, iL = cols.index("RHO"), cols.index("ALL") if "ALL" in cols else 1
rho = float(max(rows, key=lambda r: float(r[iL]))[iR])
print(f"current run: best rho = {rho:g}")

B = read_obs(f"{a.ref_dir}/B73.tsv.gz")
d_asm = collections.defaultdict(list)
for item in a.assemblies.split(","):
    name, tax = item.split(":")[:2]
    d_asm[tax].append(read_obs(f"{a.ref_dir}/{name}.tsv.gz") / np.maximum(B, 1))
d_asm = {t: np.mean(v, axis=0) for t, v in d_asm.items()}
taxa = a.taxa.split(",") if a.taxa else sorted(d_asm)

N2 = 2 * a.plants
K = np.arange(1, N2)
wk = np.array([math.comb(N2, k) for k in K], float); wk /= wk.sum()
lo, hi, st = (float(x) for x in a.d_grid.split(":"))
DG = np.round(np.arange(lo, hi + 1e-9, st), 4)
LOG_OUT = math.log(a.eps / 2.0)

pools = collections.defaultdict(list)          # taxon -> [(sample, r, ok, a, lam, sigma)]
for f in sorted(glob.glob(f"{a.hmm_dir}/*.npz")):
    with np.load(f, allow_pickle=True) as z:
        if str(z["role"]) != "donor" or float(z["sigma"]) > a.noisy_sigma: continue
        t = str(z["taxon"])
        if t not in taxa: continue
        j = int(np.argmin(np.abs(z["rho_grid"] - rho)))
        pools[t].append((os.path.basename(f)[:-4], z["r"].astype(np.float32), z["ok"].astype(bool),
                         float(z["par_best"][j]), float(z["lam"]), float(z["sigma"])))


def npdf(x, m, s):
    return np.exp(-0.5 * ((x - m) / s) ** 2) / (s * 2.5066282746)


summ = []
for t in taxa:
    P = pools.get(t, [])
    if not P: print(f"{t}: no pools"); continue
    R = np.vstack([p[1] for p in P]).astype(float)                 # pools x windows
    OK = np.vstack([p[2] for p in P]) & np.isfinite(R)
    A = np.array([p[3] for p in P])[:, None]; LAM = np.array([p[4] for p in P]); SIG = np.array([p[5] for p in P])[:, None]
    G = (1 - np.exp(-(1 - K[None, :] / N2) * LAM[:, None])) / (1 - np.exp(-LAM[:, None]))   # pools x k
    W = R.shape[1]
    # window noise factor: robust spread of r among pools near 1 (B73-like), relative to pool sigma
    with np.errstate(invalid="ignore"):
        Z = np.where(OK & (np.abs(R - 1) < 4 * SIG), (R - 1) / SIG, np.nan)
        wf = np.clip(1.4826 * np.nanmedian(np.abs(Z), axis=0), 0.5, 3)
    wf = np.where(np.isfinite(wf), wf, 1.0)
    best_d = d_asm[t].copy(); n_intro = np.zeros(W); gain = np.zeros(W)
    for c0 in range(0, W, a.chunk):
        sl = slice(c0, min(W, c0 + a.chunk))
        r = R[:, sl]; ok = OK[:, sl]; s = SIG * wf[None, sl]
        r0 = np.where(ok, r, 1.0)
        fB = (1 - a.eps) * npdf(r0, 1.0, s) + a.eps / 2
        bestll = np.full(r.shape[1], -np.inf); bd = np.full(r.shape[1], np.nan); bn = np.zeros(r.shape[1])
        for d in DG:
            ad = np.clip(A * d, 0.01, 0.95)                                        # pools x 1
            fD = (1 - a.eps) * npdf(r0, ad, s) + a.eps / 2
            mu = ad[:, :, None] + (1 - ad[:, :, None]) * G[:, None, :]             # pools x 1 x k
            fS = (1 - a.eps) * np.sum(wk[None, None, :] * npdf(r0[:, :, None], mu, s[:, :, None]), axis=2) + a.eps / 2
            pi = np.tile(np.array([0.85, 0.06, 0.09])[:, None], (1, r.shape[1]))   # per-window weights, EM
            for _ in range(a.em_iter):
                num = np.stack([pi[0] * fB, pi[1] * fS, pi[2] * fD])               # 3 x pools x win
                tot = num.sum(0); resp = np.where(ok, num / tot, 0.0)
                pi = np.clip(resp.sum(1) / np.maximum(ok.sum(0), 1), 1e-4, None); pi /= pi.sum(0)
            ll = np.where(ok, np.log((pi[0] * fB + pi[1] * fS + pi[2] * fD)), 0.0).sum(0)
            better = ll > bestll
            bestll[better] = ll[better]; bd[better] = d
            bn[better] = (resp[1] + resp[2]).sum(0)[better]
        # likelihood at the assembly level, for the gain
        ad = np.clip(A * d_asm[t][None, sl], 0.01, 0.95)
        fD = (1 - a.eps) * npdf(r0, ad, s) + a.eps / 2
        mu = ad[:, :, None] + (1 - ad[:, :, None]) * G[:, None, :]
        fS = (1 - a.eps) * np.sum(wk[None, None, :] * npdf(r0[:, :, None], mu, s[:, :, None]), axis=2) + a.eps / 2
        pi = np.tile(np.array([0.85, 0.06, 0.09])[:, None], (1, r.shape[1]))
        for _ in range(a.em_iter):
            num = np.stack([pi[0] * fB, pi[1] * fS, pi[2] * fD]); resp = np.where(ok, num / num.sum(0), 0.0)
            pi = np.clip(resp.sum(1) / np.maximum(ok.sum(0), 1), 1e-4, None); pi /= pi.sum(0)
        ll0 = np.where(ok, np.log(pi[0] * fB + pi[1] * fS + pi[2] * fD), 0.0).sum(0)
        fit = np.where(np.isfinite(bd), bd, d_asm[t][sl])
        best_d[sl] = (bn * fit + a.k * d_asm[t][sl]) / (bn + a.k)
        n_intro[sl] = bn; gain[sl] = bestll - ll0
    best_d = np.clip(best_d, 0.01, 0.95)
    np.save(f"{a.out_dir}/{t}.npy", best_d.astype(np.float32))
    with gzip.open(f"{a.out_dir}/levels_{t}.tsv.gz", "wt") as f:
        f.write("window\tassembly_level\tpool_level\tintrogressed_pools\tloglik_gain\n")
        for i in range(W):
            f.write(f"{i}\t{d_asm[t][i]:.4f}\t{best_d[i]:.4f}\t{n_intro[i]:.1f}\t{gain[i]:.2f}\n")
    diff = best_d - d_asm[t]; used = n_intro >= 3
    summ.append({"taxon": t, "pools": len(P), "windows": W, "windows_with_>=3_introgressed_pools": int(used.sum()),
                 "median_change": round(float(np.median(diff[used])), 4) if used.any() else 0,
                 "windows_raised_>0.05": int((diff > 0.05).sum()), "windows_lowered_>0.05": int((diff < -0.05).sum()),
                 "median_loglik_gain": round(float(np.median(gain[used])), 2) if used.any() else 0})
    print(summ[-1])
with open(f"{a.out_dir}/summary.tsv", "w") as f:
    keys = list(summ[0].keys()) if summ else []
    f.write("\t".join(keys) + "\n")
    for s_ in summ: f.write("\t".join(str(s_[k]) for k in keys) + "\n")
print("wrote", a.out_dir)
