#!/usr/bin/env python3
"""
new_model / 20n — copy of final_script/20_kmer_hmm.py with one addition: --level_dir.

  --level_dir DIR: DONOR level per window from <DIR>/<TAXON>.npy (new_model/20c: learned from all pools of
  the taxon) instead of the mean of the taxon's assemblies. Everything else is unchanged: the per-pool scale
  grid, d_spread (tau from the assemblies), SEG emission, transitions, rho grid and decoding.
  --acc_level_dir DIR (with --meta): per-accession levels from new_model/20d, <ACCESSION>.npy; a pool uses its
  accession's level if the file exists, otherwise the taxon level (--level_dir, or the assemblies).
  Without --level_dir and --acc_level_dir this script gives the same results as step 20.

--- original description (step 20) ---
20 — Introgression HMM on k-mer windows (Nirwan's BC2S3 HMM; new emission).

Unchanged from 14_hmm_pooled.py (imported, not copied): NAM genetic map -> Haldane theta,
transition matrix build_A (switch s = rho * theta, target state by prior pi), forward-backward,
posterior + hysteresis decoding with the boundary fix.

New: the observation in each 100-kb window is r = (B73 single-copy k-mers seen in the pool) /
(what a B73 pool of the same depth sees), from steps 17 (pools) and 19 (assemblies).

  Reference: non-negative least-squares fit of the B73 controls to the pool (a control is never in
    its own reference), refit without windows below 1 - CUT*sigma until stable.
  B73   : r ~ Normal(1, sigma_w)
  DONOR : r ~ Normal(a * d_w, sigma_w)
          d_w = share of the window's B73 single-copy k-mers present in the taxon's genome
          assemblies (step 19, mean of the taxon's assemblies); a = per-pool scale, fitted.
          Taxa without an assembly (Zm): d_w = constant d, fitted.
  SEG   : the pool is a mix of plants from a segregating family. With n plants, donor copies in
          the pool ~ Binomial(2n, 1/2) (sum of n 1:2:1 plants), so donor fraction f = k/2n.
          k-mers the donor lacks come only from the B73 part of the pool:
          r ~ Normal(d + (1 - d) * g(f), sigma_w),  g(f) = (1 - exp(-(1-f)lam)) / (1 - exp(-lam)),
          lam = the pool's k-mer coverage (Lander-Waterman: share seen = 1 - exp(-lam)).
          Averaged over k = 1..2n-1 (binomial weights) and n = PLANTS (equal weights).
  Every state: 2% chance the window is an outlier (uniform 0-2).
  --d_spread (option, off by default): our lines' donor plant is not an assembly plant, so its level
          in a window can differ from the assembly mean. tau_w = how much the taxon's assemblies
          disagree in that window (SD of the two assemblies' d, x sqrt(1.5) for a new plant vs the mean
          of two, averaged over +-D_SPREAD_HALF windows). DONOR sd = sqrt(sigma_w^2 + (a tau_w)^2);
          SEG sd = sqrt(sigma_w^2 + ((1 - g(f)) a tau_w)^2) (SEG moves with d by (1 - g)).
          Taxa with one assembly (Zh, Zl): tau = 0. B73 controls: tau of the widest taxon (worst case).
  sigma_w = pool sigma (robust SD of r in B73 windows) x how much the controls disagree in that
          window relative to a typical window (clipped 0.5-3).

  rho is not fixed: each pool is run over RHO_GRID and its log-likelihood saved; step 21 picks the
  rho with the highest total log-likelihood over all pools and writes the calls for that rho.
  Calls are saved with hysteresis R = 1 (none) and R = 3 windows, to test whether R is needed.

Output per pool: <out_dir>/<SAMPLE>.npz
"""
import argparse
import gzip
import importlib.util
import math
import os
import sys

import numpy as np
from scipy.optimize import nnls

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("hmm14", os.path.join(HERE, "..", "final_script", "14_hmm_pooled.py"))
H = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(H)


def read_windows(path):
    chrom, start, obs = [], [], []
    with gzip.open(path, "rt") as f:
        head = f.readline().rstrip("\n").split("\t")
        ic, istart, iobs = head.index("seqname"), head.index("start"), head.index("observed_kmers")
        for line in f:
            p = line.rstrip("\n").split("\t")
            chrom.append(p[ic]); start.append(int(p[istart])); obs.append(float(p[iobs]))
    return np.array(chrom), np.array(start), np.array(obs)


def seg_fractions(plants):
    """Donor fractions of a SEG pool and their weights: k/2n, k ~ Binomial(2n, 1/2), k = 1..2n-1."""
    w = {}
    for n in plants:
        tot = sum(math.comb(2 * n, k) for k in range(1, 2 * n))
        for k in range(1, 2 * n):
            f = round(k / (2 * n), 6)
            w[f] = w.get(f, 0.0) + math.comb(2 * n, k) / tot / len(plants)
    fs = np.array(sorted(w)); return fs, np.array([w[f] for f in fs])


def norm_logpdf(x, m, s):
    return -0.5 * ((x - m) / s) ** 2 - np.log(s) - 0.9189385332


def emission_ll(r, sig, dw, lam, fs, fw, eps, tau=None):
    """log P(r | state) per window, columns B73, SEG, DONOR. tau: per-window SD of the DONOR level (or None)."""
    out = math.log(eps / 2.0)
    mix = lambda lp: np.logaddexp(math.log(1 - eps) + lp, out)
    tau = np.zeros_like(r) if tau is None else tau
    lB = mix(norm_logpdf(r, 1.0, sig))
    lD = mix(norm_logpdf(r, dw, np.sqrt(sig ** 2 + tau ** 2)))
    g = (1 - np.exp(-(1 - fs) * lam)) / (1 - math.exp(-lam))                     # per fraction
    mu = dw[:, None] + (1 - dw[:, None]) * g[None, :]                            # windows x fractions
    sS = np.sqrt(sig[:, None] ** 2 + (tau[:, None] * (1 - g[None, :])) ** 2)
    lS = np.log(np.sum(fw[None, :] * np.exp(norm_logpdf(r[:, None], mu, sS)), axis=1) + 1e-300)
    lS = mix(lS)
    return np.vstack([lB, lS, lD]).T


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", required=True, help="comma-separated SAMPLE names")
    ap.add_argument("--pools_tsv", required=True, help="step 11 pools.tsv (SAMPLE ROLE GROUP ...)")
    ap.add_argument("--win_dir", required=True, help="step 17 windows: <SAMPLE>.tsv.gz")
    ap.add_argument("--ref_dir", required=True, help="step 19: B73.tsv.gz + <NAME>.tsv.gz per assembly")
    ap.add_argument("--assemblies", required=True, help="NAME:TAXON,... (DONOR_ASSEMBLIES)")
    ap.add_argument("--map_dir", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--rho_grid", default="0.25,0.5,1,1.5,2,3,4,8")
    ap.add_argument("--scale_grid", default="0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.6,1.8")
    ap.add_argument("--const_d_grid", default="0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8", help="taxa without an assembly")
    ap.add_argument("--plants", default="3,4,5,6")
    ap.add_argument("--prior", default="0.844,0.0625,0.094", help="B73,SEG,DONOR (BC2S3 expectation)")
    ap.add_argument("--cut_sigma", type=float, default=2.0)
    ap.add_argument("--min_b73_kmers", type=float, default=2000)
    ap.add_argument("--eps", type=float, default=0.02)
    ap.add_argument("--rigidity", type=int, default=3)
    ap.add_argument("--family_dir", default="", help="family DONOR levels (20b): <SAMPLE>.npy (leave-one-out); replaces the assembly level")
    ap.add_argument("--meta", default="", help="SAMPLE_META, for pool -> family (with --family_dir)")
    ap.add_argument("--family_scale_grid", default="0.8,0.9,1.0,1.1,1.2,1.3")
    ap.add_argument("--d_spread", action="store_true", help="widen DONOR/SEG by the between-assembly spread of d (see top)")
    ap.add_argument("--d_spread_half", type=int, default=5, help="windows each side to average the spread over")
    ap.add_argument("--level_dir", default="", help="new_model/20c levels: <TAXON>.npy replaces the assembly mean")
    ap.add_argument("--acc_level_dir", default="", help="new_model/20d levels: <ACCESSION>.npy (needs --meta)")
    a = ap.parse_args()

    rho_grid = [float(x) for x in a.rho_grid.split(",")]
    scale_grid = [float(x) for x in a.scale_grid.split(",")]
    const_grid = [float(x) for x in a.const_d_grid.split(",")]
    plants = [int(x) for x in a.plants.split(",")]
    pi = [float(x) for x in a.prior.split(",")]
    fs, fw = seg_fractions(plants)

    role, group = {}, {}
    with open(a.pools_tsv) as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            role[p[0]], group[p[0]] = p[1], p[2]
    controls = [s for s in role if role[s] == "control" and os.path.exists(f"{a.win_dir}/{s}.tsv.gz")]
    if len(controls) < 3:
        sys.exit(f"need >= 3 B73 controls with windows, found {len(controls)}")

    chrom, start, B = read_windows(f"{a.ref_dir}/B73.tsv.gz")
    taxon_asm = {}
    for item in a.assemblies.split(","):
        name, tax = item.split(":")[:2]
        _c, _s, obs = read_windows(f"{a.ref_dir}/{name}.tsv.gz")
        taxon_asm.setdefault(tax, []).append(obs / np.maximum(B, 1))
    d_tax = {t: np.mean(v, axis=0) for t, v in taxon_asm.items()}
    if a.level_dir:
        for t in list(d_tax):
            p = f"{a.level_dir}/{t}.npy"
            if os.path.exists(p):
                lv = np.load(p).astype(float)
                if len(lv) != len(B): sys.exit(f"{p}: {len(lv)} windows, expected {len(B)}")
                print(f"{t}: DONOR level from {p} (mean change vs assemblies {np.mean(lv - d_tax[t]):+.4f})")
                d_tax[t] = lv
    d_all = np.mean([v for vs in taxon_asm.values() for v in vs], axis=0)       # for controls
    tau_tax = {}
    if a.d_spread:
        okB = B >= a.min_b73_kmers
        for t, v in taxon_asm.items():
            if len(v) < 2: tau_tax[t] = np.zeros(len(B)); continue
            var = np.var(np.vstack(v), axis=0, ddof=1) * (1 + 1 / len(v))     # a new plant vs the mean of n
            var = np.where(okB, var, np.nan); k = a.d_spread_half; sm = np.empty(len(B))
            for i in range(len(B)):                                            # smooth along the genome
                w = var[max(0, i - k): i + k + 1]; sm[i] = np.nanmean(w) if np.isfinite(w).any() else 0.0
            tau_tax[t] = np.sqrt(sm)
        widest = max(tau_tax, key=lambda t: np.nanmean(tau_tax[t]))
        tau_tax["__control__"] = tau_tax[widest]
        print("d_spread: mean tau per taxon " + ", ".join(f"{t} {np.mean(v):.3f}" for t, v in tau_tax.items() if t != "__control__")
              + f" | controls use {widest}")

    def load(s):
        c, st, o = read_windows(f"{a.win_dir}/{s}.tsv.gz")
        if len(o) != len(B) or not (st == start).all():
            sys.exit(f"{s}: windows do not match B73.tsv.gz")
        return o
    C = {s: load(s) for s in controls}
    Cn = {s: v / np.median(v) for s, v in C.items()}

    acc_of = {}
    if a.acc_level_dir:
        if not a.meta: sys.exit("--acc_level_dir needs --meta")
        with open(a.meta, encoding="utf-8", errors="ignore") as f:
            next(f)
            for line in f:
                p = line.rstrip("\r\n").split("\t")
                if len(p) > 1: acc_of[p[0]] = p[1].split("_")[0]
    family = {}
    if a.family_dir:
        with open(a.meta, encoding="utf-8", errors="ignore") as f:
            next(f)
            for line in f:
                p = line.rstrip("\r\n").split("\t")
                if len(p) > 1: family[p[0]] = "_".join(p[1].split("_")[:2])
    fam_grid = [float(x) for x in a.family_scale_grid.split(",")]

    maps = {c: H.read_map(f"{a.map_dir}/{c}.map") for c in np.unique(chrom)}
    cm = np.array([H.interp_cm(*maps[c], int(p) + 50000) for c, p in zip(chrom, start)])
    os.makedirs(a.out_dir, exist_ok=True)

    for S in a.samples.split(","):
        if not os.path.exists(f"{a.win_dir}/{S}.tsv.gz"):
            print(f"{S}: no windows (step 17)"); continue
        O = load(S); N = O / np.median(O)
        refs = [s for s in controls if s != S]
        X = np.vstack([Cn[s] for s in refs]).T
        Rm = np.median(X, axis=1)
        spread = np.median(np.abs(X / np.where(Rm > 0, Rm, 1)[:, None] - 1), axis=1) * 1.4826
        ok = (B >= a.min_b73_kmers) & (Rm > 0.2) & np.isfinite(spread)
        # reference: fit controls to the pool, refit without windows below 1 - cut*sigma
        keep = ok.copy(); fit = None
        for _ in range(8):
            b, _ = nnls(X[keep], N[keep]); fit = X @ b
            r_all = N / np.where(fit > 0, fit, np.nan)
            sig = 1.4826 * np.nanmedian(np.abs(r_all[keep] - 1))
            new = ok & (r_all > 1 - a.cut_sigma * sig)
            if (new == keep).all(): break
            keep = new
        r = r_all / np.nanmedian(r_all[keep])
        sig_p = 1.4826 * np.nanmedian(np.abs(r[keep] - 1))
        wscale = np.clip(spread / np.nanmedian(spread[ok]), 0.5, 3)
        frac = np.median(O[ok] / B[ok]); lam = -math.log(max(1e-6, 1 - min(frac, 0.99)))
        tax = group.get(S, "NA")
        famfile = f"{a.family_dir}/{S}.npy" if a.family_dir else ""
        if famfile and os.path.exists(famfile):
            dw0, grid, kind = np.load(famfile).astype(float), fam_grid, "family_scale"
        elif tax in d_tax:
            dw0, grid, kind = d_tax[tax], scale_grid, "assembly_scale"
            accf = f"{a.acc_level_dir}/{acc_of.get(S, '')}.npy" if a.acc_level_dir else ""
            if accf and acc_of.get(S) and os.path.exists(accf):
                dw0 = np.load(accf).astype(float)
        elif role.get(S) == "control": dw0, grid, kind = d_all, scale_grid, "assembly_scale"
        else:              dw0, grid, kind = None, const_grid, "constant_d"

        idx_by_chr = {c: np.where((chrom == c) & ok)[0] for c in np.unique(chrom)}
        th_by_chr = {c: [0.0] + [H.haldane_theta(max(0.0, (cm[i[k]] - cm[i[k - 1]]) / 100.0)) for k in range(1, len(i))]
                     for c, i in idx_by_chr.items()}
        sig_w = sig_p * wscale
        tau0 = None
        if a.d_spread and kind == "assembly_scale":
            tau0 = tau_tax.get(tax) if tax in tau_tax else tau_tax["__control__"]
        em_by_par = {}
        for par in grid:
            dw = np.clip(par * dw0 if dw0 is not None else np.full(len(B), par), 0.01, 0.95)
            tau = None if tau0 is None else np.where(ok, par * tau0, 0.0)
            em_by_par[par] = emission_ll(np.where(ok, r, 1.0), np.where(ok, sig_w, 1.0), dw, lam, fs, fw, a.eps, tau)

        W = len(B)
        ll = np.zeros(len(rho_grid)); par_best = np.zeros(len(rho_grid))
        st1 = np.full((len(rho_grid), W), -1, np.int8); stR = np.full((len(rho_grid), W), -1, np.int8)
        for j, rho in enumerate(rho_grid):
            best = None
            for par in grid:
                e = em_by_par[par]; tot = 0.0; posts = {}
                for c, idx in idx_by_chr.items():
                    if len(idx) == 0: continue
                    post, l = H.forward_backward(e[idx].tolist(), th_by_chr[c], pi, rho)
                    tot += l; posts[c] = post
                if best is None or tot > best[0]: best = (tot, par, posts)
            ll[j], par_best[j] = best[0], best[1]
            for c, post in best[2].items():
                idx = idx_by_chr[c]; mx = [int(np.argmax(p)) for p in post]
                st1[j, idx] = mx
                stR[j, idx] = H.decode_hysteresis(mx, a.rigidity)
        np.savez_compressed(f"{a.out_dir}/{S}.npz", rho_grid=np.array(rho_grid), loglik=ll, par_best=par_best,
                            param_kind=kind, states_R1=st1, states_R=stR, rigidity=a.rigidity, r=r.astype(np.float32),
                            ok=ok, sigma=sig_p, lam=lam, taxon=tax, role=role.get(S, "NA"), d_spread=bool(tau0 is not None),
                            level_source=("accession" if (a.acc_level_dir and acc_of.get(S) and os.path.exists(f"{a.acc_level_dir}/{acc_of.get(S)}.npy"))
                                          else "pools" if (a.level_dir and os.path.exists(f"{a.level_dir}/{tax}.npy")) else "assemblies"))
        jb = int(np.argmax(ll)); v = st1[jb][st1[jb] >= 0]
        print(f"{S} ({tax}): sigma {sig_p:.3f} | best rho for this pool {rho_grid[jb]} | {kind} {par_best[jb]} | "
              f"SEG {100*np.mean(v==1):.2f}% DONOR {100*np.mean(v==2):.2f}%")


if __name__ == "__main__":
    main()
