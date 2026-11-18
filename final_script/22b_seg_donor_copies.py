#!/usr/bin/env python3
"""
22b — How many donor copies do SEG pools carry? (README section 2.8, "Checking the deficit").

A SEG pool comes from a heterozygous (AB) BC2S2 mother. Its BC2S3 plants were grown and pooled without
selection, so under Mendel the donor allele share in the pool is 1/2: donor copies k ~ Binomial(2n, 1/2).
If teosinte alleles are lost within a generation (segregation distortion, or AB/BB seeds and seedlings
doing worse than AA), SEG pools carry fewer donor copies than 1/2. Selection on which plant is advanced
(the mother) does not change this share, so the check separates the two explanations:
  share ~ 1/2 -> no loss within the last generation; the DONOR deficit is not from distortion/viability
  share < 1/2 -> teosinte alleles are lost within a generation, which would also have acted in earlier ones

Estimate (method of moments, per pool): in the windows called SEG, the expected r under a donor allele
share p is   E[r] = d_w + (1 - d_w) * sum_k w_k(p) g(k/2n),   w_k(p) = Binomial(2n, p) over k = 1..2n-1
(rescaled; k = 0 and 2n are not SEG), g = the k-mer detection curve of step 20, d_w = a * assembly level.
p is solved so that the mean expected r equals the mean observed r. Checks on the same pools:
  DONOR windows: mean r - a*d_w   (is the DONOR level right? an offset here would also bias p)
  B73 windows:   mean r - 1
--selftest simulates SEG windows at p = 0.5 with the same noise and shows that the estimator returns 0.5.

Usage (cluster, after steps 20-21):
  python final_script/22b_seg_donor_copies.py --hmm_dir $D20 --ref_dir $D19 --pool_summary $D21/pool_summary.tsv \
      --assemblies "$ASMS" --plants 6 --out $D21/seg_donor_copies.tsv [--selftest]
  ASMS = NAME:TAXON,... as in 20_run_kmer_hmm.sh
Outputs: per-pool table (--out) and a per-taxon summary printed to the log (also <out>.taxon.tsv).
Needs: numpy, pandas, scipy.
"""
import argparse
import gzip
import math
import numpy as np
import pandas as pd
from scipy.optimize import brentq

ap = argparse.ArgumentParser()
ap.add_argument("--hmm_dir", required=True, help="step 20 output: <SAMPLE>.npz")
ap.add_argument("--ref_dir", required=True, help="step 19: B73.tsv.gz + <NAME>.tsv.gz")
ap.add_argument("--pool_summary", required=True, help="step 21 pool_summary.tsv (FLAG, PARAM)")
ap.add_argument("--assemblies", required=True, help="NAME:TAXON,...")
ap.add_argument("--plants", type=int, default=6)
ap.add_argument("--rho_loglik", default="", help="step 21 rho_loglik.tsv (default: <pool_summary dir>/rho_loglik.tsv)")
ap.add_argument("--out", required=True)
ap.add_argument("--selftest", action="store_true")
ap.add_argument("--seed", type=int, default=1)
a = ap.parse_args()
rng = np.random.default_rng(a.seed)
N2 = 2 * a.plants
K = np.arange(1, N2)


def read_obs(path):
    with gzip.open(path, "rt") as f:
        head = f.readline().rstrip("\n").split("\t"); i = head.index("observed_kmers")
        return np.array([float(l.rstrip("\n").split("\t")[i]) for l in f])


B = read_obs(f"{a.ref_dir}/B73.tsv.gz")
asm = {}
for item in a.assemblies.split(","):
    name, tax = item.split(":")[:2]
    asm.setdefault(tax, []).append(read_obs(f"{a.ref_dir}/{name}.tsv.gz") / np.maximum(B, 1))
d_tax = {t: np.mean(v, axis=0) for t, v in asm.items()}


def wk(p):
    w = np.array([math.comb(N2, k) * p ** k * (1 - p) ** (N2 - k) for k in K]); return w / w.sum()


def g(f, lam):
    return (1 - np.exp(-(1 - f) * lam)) / (1 - math.exp(-lam))


def expected_r(p, dw, lam):
    gm = np.sum(wk(p) * g(K / N2, lam))           # mean detection of the B73 part over k
    return dw + (1 - dw) * gm


def solve_p(r, dw, lam):
    """p such that mean expected r over the windows = mean observed r (clipped to the reachable range)"""
    target = np.mean(r); f = lambda p: np.mean(expected_r(p, dw, lam)) - target
    lo, hi = 0.02, 0.98
    if f(lo) < 0: return lo                         # E[r] falls as p rises; target above reach -> lo
    if f(hi) > 0: return hi
    return brentq(f, lo, hi)


ps = pd.read_csv(a.pool_summary, sep="\t")
ps = ps[(ps.ROLE == "donor") & (ps.FLAG == "ok") & ps.GROUP.isin(d_tax)]
rl = a.rho_loglik or a.pool_summary.rsplit("/", 1)[0] + "/rho_loglik.tsv"
rho_best = float(pd.read_csv(rl, sep="\t").sort_values("ALL").RHO.iloc[-1])

rows = []
for s in ps.itertuples():
    with np.load(f"{a.hmm_dir}/{s.SAMPLE}.npz", allow_pickle=True) as z:
        j = int(np.flatnonzero(np.isclose(z["rho_grid"], rho_best))[0])
        st, r, ok, lam = z["states_R1"][j], z["r"].astype(float), z["ok"], float(z["lam"])
        par = float(z["par_best"][j]); sig = float(z["sigma"])
    dw = np.clip(par * d_tax[s.GROUP], 0.01, 0.95)
    seg, don, b73 = ok & (st == 1), ok & (st == 2), ok & (st == 0)
    if seg.sum() < 20: continue
    row = {"SAMPLE": s.SAMPLE, "GROUP": s.GROUP, "lam": round(lam, 3), "sigma": round(sig, 3), "scale_a": par,
           "n_seg": int(seg.sum()), "n_donor": int(don.sum()),
           "p_seg": solve_p(r[seg], dw[seg], lam),
           "r_seg_obs": r[seg].mean(), "r_seg_exp_half": np.mean(expected_r(0.5, dw[seg], lam)),
           "donor_offset": (r[don] - dw[don]).mean() if don.sum() >= 20 else np.nan,
           "b73_offset": (r[b73] - 1).mean()}
    if a.selftest:   # simulate the SEG windows at p = 0.5 with the pool's noise, re-estimate
        k = rng.binomial(N2, 0.5, seg.sum())
        while ((k == 0) | (k == N2)).any():
            bad = (k == 0) | (k == N2); k[bad] = rng.binomial(N2, 0.5, bad.sum())
        rs = dw[seg] + (1 - dw[seg]) * g(k / N2, lam) + rng.normal(0, sig, seg.sum())
        row["p_selftest"] = solve_p(rs, dw[seg], lam)
    rows.append(row)

D = pd.DataFrame(rows)
D.to_csv(a.out, sep="\t", index=False, float_format="%.4f")


def summ(x):
    boot = [np.median(rng.choice(x.p_seg.values, len(x))) for _ in range(2000)]
    out = {"pools": len(x), "SEG_windows": int(x.n_seg.sum()),
           "median_p_seg": round(x.p_seg.median(), 3), "CI95": f"{np.percentile(boot, 2.5):.3f}-{np.percentile(boot, 97.5):.3f}",
           "mean_p_seg_weighted": round(np.average(x.p_seg, weights=x.n_seg), 3),
           "median_donor_offset": round(x.donor_offset.median(), 3), "median_b73_offset": round(x.b73_offset.median(), 3)}
    if "p_selftest" in x: out["median_p_selftest"] = round(x.p_selftest.median(), 3)
    return out


T = pd.DataFrame([{"taxon": t, **summ(x)} for t, x in D.groupby("GROUP")] + [{"taxon": "All", **summ(D)}])
T.to_csv(a.out + ".taxon.tsv", sep="\t", index=False)
print(f"rho = {rho_best:g}; {len(D)} pools with >= 20 SEG windows; plants per pool = {a.plants} (expected p = 0.5)")
print(T.to_string(index=False))
