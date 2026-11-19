#!/usr/bin/env python3
"""
new_model / 20d — DONOR level per window and ACCESSION, from the accession's own pools (no calls used).

Why: with one level per taxon (20c), the SEG-between-DONOR errors stayed (README of this folder, results).
In the same window, one family's teosinte can keep more B73 k-mers than another's: the level depends on the
donor plant (haplotype), not only on the taxon. Each accession's lines come from its own teosinte plants.

Same fit as 20c, on the pools of one accession: in each window their scores r are a mixture of
    B73 ~ N(1, s),  SEG ~ N(a d + (1 - a d) g(f), s) over f = k/12 (k ~ Binomial(12, 1/2)),  DONOR ~ N(a d, s)
(+ 2% outliers), with the window level d chosen on a grid by maximum likelihood. An accession has few pools
(median 11; in a window often 1-3 carry teosinte), so two things keep the fit stable:
  - the three group sizes are fitted by EM with PI_PSEUDO pseudo-pools at the taxon-wide mix (0.85, 0.06, 0.09);
  - the level is shrunk to the taxon level from 20c by the evidence:
        level = (n d_hat + K d_taxon) / (n + K),  n = expected number of the accession's pools carrying teosinte
    in the window, K (default 2). With n = 0 the accession keeps the taxon level.
  Simulation (12 pools per accession, a 0.18 shift in one accession): K = 1 / 2 / 3 leave -0.05 / -0.07 / -0.09
  of the shift, with 90% of other windows within +-0.058 / 0.047 / 0.039 of the truth (other accession: unchanged).
Accessions with fewer than --min_pools pools are not written (20n then uses the taxon level).
Accession = first part of the pool name (e.g. Zx.0040 in Zx.0040_P1_P4_P3.1.1.1), from SAMPLE_META.

Output (--out_dir): <ACCESSION>.npy (float32, one level per window), summary.tsv (per accession).
Usage: python new_model/20d_accession_levels.py --hmm_dir $D20 --rho_file $D21/rho_loglik.tsv \
          --taxon_levels $DNEW/levels --meta $SAMPLE_META --out_dir $DNEW/levels_acc [--taxa Zx] [--k 2]
Needs: numpy.
"""
import argparse
import collections
import glob
import math
import os
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--hmm_dir", required=True, help="step 20 output (<SAMPLE>.npz) of the current run")
ap.add_argument("--rho_file", required=True, help="step 21 rho_loglik.tsv of the same run")
ap.add_argument("--taxon_levels", required=True, help="20c output: <TAXON>.npy (the prior for each accession)")
ap.add_argument("--meta", required=True, help="SAMPLE_META (Sequencing_ID, new_genotype)")
ap.add_argument("--out_dir", required=True)
ap.add_argument("--taxa", default="")
ap.add_argument("--plants", type=int, default=6)
ap.add_argument("--noisy_sigma", type=float, default=0.15)
ap.add_argument("--min_pools", type=int, default=4)
ap.add_argument("--k", type=float, default=2.0)
ap.add_argument("--pi_pseudo", type=float, default=5.0)
ap.add_argument("--d_grid", default="0.02:0.95:0.01")
ap.add_argument("--eps", type=float, default=0.02)
ap.add_argument("--em_iter", type=int, default=15)
ap.add_argument("--chunk", type=int, default=3000)
a = ap.parse_args()
os.makedirs(a.out_dir, exist_ok=True)

with open(a.rho_file) as f:
    cols = f.readline().rstrip("\n").split("\t"); rows = [l.rstrip("\n").split("\t") for l in f]
iR, iL = cols.index("RHO"), (cols.index("ALL") if "ALL" in cols else 1)
rho = float(max(rows, key=lambda r: float(r[iL]))[iR])

acc_of = {}
with open(a.meta, encoding="utf-8", errors="ignore") as f:
    next(f)
    for line in f:
        p = line.rstrip("\r\n").split("\t")
        if len(p) > 1: acc_of[p[0]] = p[1].split("_")[0]

N2 = 2 * a.plants; K = np.arange(1, N2)
wk = np.array([math.comb(N2, k) for k in K], float); wk /= wk.sum()
lo, hi, st = (float(x) for x in a.d_grid.split(":"))
DG = np.round(np.arange(lo, hi + 1e-9, st), 4)
PI0 = np.array([0.85, 0.06, 0.09])
taxa = set(a.taxa.split(",")) if a.taxa else None

groups = collections.defaultdict(list)      # accession -> pools
tax_of = {}
for f in sorted(glob.glob(f"{a.hmm_dir}/*.npz")):
    s = os.path.basename(f)[:-4]
    with np.load(f, allow_pickle=True) as z:
        if str(z["role"]) != "donor" or float(z["sigma"]) > a.noisy_sigma: continue
        t = str(z["taxon"])
        if taxa and t not in taxa: continue
        if s not in acc_of or not os.path.exists(f"{a.taxon_levels}/{t}.npy"): continue
        j = int(np.argmin(np.abs(z["rho_grid"] - rho)))
        groups[acc_of[s]].append((z["r"].astype(np.float32), z["ok"].astype(bool), float(z["par_best"][j]),
                                  float(z["lam"]), float(z["sigma"])))
        tax_of[acc_of[s]] = t
print(f"best rho {rho:g}; {len(groups)} accessions, {sum(map(len, groups.values()))} pools")


def npdf(x, m, s):
    return np.exp(-0.5 * ((x - m) / s) ** 2) / (s * 2.5066282746)


# window noise factor per taxon (robust spread of r near 1 over all the taxon's pools), as in 20c
wf_tax = {}
for t in set(tax_of.values()):
    R = np.vstack([p[0] for acc, P in groups.items() if tax_of[acc] == t for p in P]).astype(float)
    OK = np.vstack([p[1] for acc, P in groups.items() if tax_of[acc] == t for p in P]) & np.isfinite(R)
    SIG = np.array([p[4] for acc, P in groups.items() if tax_of[acc] == t for p in P])[:, None]
    with np.errstate(invalid="ignore"):
        Z = np.where(OK & (np.abs(R - 1) < 4 * SIG), (R - 1) / SIG, np.nan)
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            w = np.clip(1.4826 * np.nanmedian(np.abs(Z), axis=0), 0.5, 3)
    wf_tax[t] = np.where(np.isfinite(w), w, 1.0)

summ = []
for acc, P in sorted(groups.items()):
    t = tax_of[acc]
    if len(P) < a.min_pools:
        summ.append({"accession": acc, "taxon": t, "pools": len(P), "written": "no (too few pools: taxon level)"}); continue
    prior = np.load(f"{a.taxon_levels}/{t}.npy").astype(float)
    R = np.vstack([p[0] for p in P]).astype(float); OK = np.vstack([p[1] for p in P]) & np.isfinite(R)
    A = np.array([p[2] for p in P])[:, None]; LAM = np.array([p[3] for p in P]); SIG = np.array([p[4] for p in P])[:, None]
    G = (1 - np.exp(-(1 - K[None, :] / N2) * LAM[:, None])) / (1 - np.exp(-LAM[:, None]))
    W = R.shape[1]; lev = prior.copy(); nint = np.zeros(W)
    for c0 in range(0, W, a.chunk):
        sl = slice(c0, min(W, c0 + a.chunk))
        r = R[:, sl]; ok = OK[:, sl]; s = SIG * wf_tax[t][None, sl]; r0 = np.where(ok, r, 1.0)
        fB = (1 - a.eps) * npdf(r0, 1.0, s) + a.eps / 2
        bestll = np.full(r.shape[1], -np.inf); bd = np.full(r.shape[1], np.nan); bn = np.zeros(r.shape[1])
        nok = ok.sum(0)
        for d in DG:
            ad = np.clip(A * d, 0.01, 0.95)
            fD = (1 - a.eps) * npdf(r0, ad, s) + a.eps / 2
            mu = ad[:, :, None] + (1 - ad[:, :, None]) * G[:, None, :]
            fS = (1 - a.eps) * np.sum(wk[None, None, :] * npdf(r0[:, :, None], mu, s[:, :, None]), axis=2) + a.eps / 2
            pi = np.tile(PI0[:, None], (1, r.shape[1]))
            for _ in range(a.em_iter):
                num = np.stack([pi[0] * fB, pi[1] * fS, pi[2] * fD]); resp = np.where(ok, num / num.sum(0), 0.0)
                pi = (resp.sum(1) + a.pi_pseudo * PI0[:, None]) / (nok + a.pi_pseudo)
            ll = np.where(ok, np.log(pi[0] * fB + pi[1] * fS + pi[2] * fD), 0.0).sum(0)
            better = ll > bestll
            bestll[better] = ll[better]; bd[better] = d; bn[better] = (resp[1] + resp[2]).sum(0)[better]
        fit = np.where(np.isfinite(bd), bd, prior[sl])
        lev[sl] = (bn * fit + a.k * prior[sl]) / (bn + a.k); nint[sl] = bn
    lev = np.clip(lev, 0.01, 0.95)
    np.save(f"{a.out_dir}/{acc}.npy", lev.astype(np.float32))
    diff = lev - prior; used = nint >= 1
    summ.append({"accession": acc, "taxon": t, "pools": len(P), "written": "yes",
                 "windows_with_>=1_introgressed_pool": int(used.sum()),
                 "median_abs_change_vs_taxon": round(float(np.median(np.abs(diff[used]))), 4) if used.any() else 0,
                 "windows_changed_>0.05": int((np.abs(diff) > 0.05).sum())})
    print(summ[-1])
keys = ["accession", "taxon", "pools", "written", "windows_with_>=1_introgressed_pool", "median_abs_change_vs_taxon", "windows_changed_>0.05"]
with open(f"{a.out_dir}/summary.tsv", "w") as f:
    f.write("\t".join(keys) + "\n")
    for s_ in summ: f.write("\t".join(str(s_.get(k, "")) for k in keys) + "\n")
print("wrote", a.out_dir)
