#!/usr/bin/env python3
"""
20b — Family DONOR level per window, measured from the family's own pools.

All pools of a BC1 family carry the same teosinte haplotype (one F1 gamete). Where pools of the
family are called DONOR in a window, their score r there is that family's real DONOR level. This
replaces the assembly level (other plants of the same taxon) in step 20.

Per family and window (non-noisy pools, calls of the previous round at its best rho):
  n = number of the family's pools called DONOR in the window, m = their median r
  level = (n * m + K * prior) / (n + K)
  prior = family scale x taxon assembly level (scale = median fitted scale of the family's pools)
K (default 3) shrinks windows with few DONOR pools toward the assembly level; with n = 0 the
level is the prior. Leave-one-out: each pool's level is computed from the OTHER pools of its
family, so a pool's own calls never set its own DONOR level.
Output: <out_dir>/<SAMPLE>.npy (float32, one value per window) + summary.tsv (per family).
"""
import argparse
import collections
import glob
import gzip
import os
import warnings

import numpy as np


def read_obs(path):
    with gzip.open(path, "rt") as f:
        head = f.readline().rstrip("\n").split("\t"); i = head.index("observed_kmers")
        return np.array([float(l.rstrip("\n").split("\t")[i]) for l in f])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hmm_dir", required=True, help="previous round's step-20 .npz")
    ap.add_argument("--rho_file", required=True, help="previous round's rho_loglik.tsv (step 21)")
    ap.add_argument("--meta", required=True)
    ap.add_argument("--ref_dir", required=True)
    ap.add_argument("--assemblies", required=True, help="NAME:TAXON,...")
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--noisy_sigma", type=float, default=0.15)
    ap.add_argument("--k", type=float, default=3.0)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)

    with open(a.rho_file) as f:
        next(f); rows = [l.split("\t") for l in f]
    rho = float(max(rows, key=lambda r: float(r[1]))[0])
    print(f"previous round best rho = {rho:g}")

    B = read_obs(f"{a.ref_dir}/B73.tsv.gz")
    d_tax = collections.defaultdict(list)
    for item in a.assemblies.split(","):
        name, tax = item.split(":")[:2]
        d_tax[tax].append(read_obs(f"{a.ref_dir}/{name}.tsv.gz") / np.maximum(B, 1))
    d_tax = {t: np.mean(v, axis=0) for t, v in d_tax.items()}

    fam = {}
    with open(a.meta, encoding="utf-8", errors="ignore") as f:
        next(f)
        for line in f:
            p = line.rstrip("\r\n").split("\t")
            if len(p) > 1: fam[p[0]] = "_".join(p[1].split("_")[:2])

    by = collections.defaultdict(list)          # family -> [(r, states, scale)]
    tax_of = {}
    for f in sorted(glob.glob(f"{a.hmm_dir}/*.npz")):
        with np.load(f, allow_pickle=True) as z:
            if str(z["role"]) != "donor" or float(z["sigma"]) > a.noisy_sigma: continue
            t = str(z["taxon"])
            if t not in d_tax: continue
            j = int(np.argmin(np.abs(z["rho_grid"] - rho)))
            s = os.path.basename(f)[:-4]; fm = fam.get(s)
            if fm is None: continue
            by[fm].append((s, z["r"].astype(float), z["states_R1"][j].copy(), float(z["par_best"][j]))); tax_of[fm] = t

    with open(f"{a.out_dir}/summary.tsv", "w") as out:
        out.write("FAMILY\tTAXON\tPOOLS\tSCALE\tWINDOWS_WITH_DONOR_POOLS\tMEDIAN_DONOR_POOLS\tMEAN_LEVEL\tMEAN_ASSEMBLY_LEVEL\n")
        for fm, lst in sorted(by.items()):
            t = tax_of[fm]; R = np.vstack([x[1] for x in lst]); S = np.vstack([x[2] for x in lst])
            scale = float(np.median([x[3] for x in lst]))
            prior = np.clip(scale * d_tax[t], 0.01, 0.95)
            D = np.where(S == 2, R, np.nan)
            levels = []
            for i, (sample, *_rest) in enumerate(lst):
                Di = np.delete(D, i, axis=0)                     # the other pools of the family
                n = np.sum(np.isfinite(Di), axis=0)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", RuntimeWarning)
                    m = np.where(n > 0, np.nanmedian(Di, axis=0) if len(Di) else 0.0, 0.0)
                lev = np.clip((n * m + a.k * prior) / (n + a.k), 0.01, 0.95).astype(np.float32)
                np.save(f"{a.out_dir}/{sample}.npy", lev); levels.append(lev)
            nall = np.sum(S == 2, axis=0); cov = nall > 0
            out.write(f"{fm}\t{t}\t{len(lst)}\t{scale:g}\t{int(cov.sum())}\t{np.median(nall[cov]) if cov.any() else 0:g}\t"
                      f"{np.mean(levels):.4f}\t{prior.mean():.4f}\n")
    print(f"families: {len(by)}  pools: {sum(map(len, by.values()))}  WROTE: {a.out_dir}/<SAMPLE>.npy, summary.tsv")


if __name__ == "__main__":
    main()
