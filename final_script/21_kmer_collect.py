#!/usr/bin/env python3
"""
21 — Collect the k-mer HMM (step 20): pick rho, write the calls, draw block figures.

rho: summed log-likelihood over donor pools that are not Zm (no assembly) and not NOISY; the rho
     with the highest total is used for every pool. rho_loglik.tsv shows the whole curve, also
     per taxon.
Calls: hysteresis R = 1 (none). pool_summary.tsv also gives the share with R = step-20 rigidity,
       to show whether R changes anything.

Outputs in <out_dir>:
  rho_loglik.tsv     rho, total log-likelihood (relative to best), and per taxon
  pool_summary.tsv   SAMPLE ROLE GROUP SIGMA FLAG PARAM_KIND PARAM FRAC_B73 FRAC_SEG FRAC_DONOR
                     DONOR_SHARE DONOR_SHARE_R N_TRACTS
  tracts_all.tsv.gz  SAMPLE GROUP CHROM START END STATE (non-B73 and B73 runs)
  blocks_<GROUP>.pdf one page per pool (blue B73, purple SEG, red DONOR), if matplotlib is present
"""
import argparse
import glob
import gzip
import os

import numpy as np

STATES = ("B73", "SEG", "DONOR")


def read_windows(path):
    chrom, start = [], []
    with gzip.open(path, "rt") as f:
        head = f.readline().rstrip("\n").split("\t")
        ic, istart = head.index("seqname"), head.index("start")
        for line in f:
            p = line.rstrip("\n").split("\t"); chrom.append(p[ic]); start.append(int(p[istart]))
    return np.array(chrom), np.array(start)


def runs(states, chrom, start, win):
    """Contiguous runs per chromosome; unscored windows (-1) take the state before them."""
    out = []
    for c in dict.fromkeys(chrom):
        idx = np.where(chrom == c)[0]; x = states[idx].copy()
        for k in range(len(x)):
            if x[k] < 0: x[k] = x[k - 1] if k > 0 else 0
        b0 = 0
        for k in range(1, len(x) + 1):
            if k == len(x) or x[k] != x[b0]:
                out.append((c, int(start[idx[b0]]), int(start[idx[k - 1]]) + win, int(x[b0]))); b0 = k
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hmm_dir", required=True)
    ap.add_argument("--ref_dir", required=True, help="step 19 dir (B73.tsv.gz gives the windows)")
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--window", type=int, default=100000)
    ap.add_argument("--noisy_sigma", type=float, default=0.15)
    ap.add_argument("--no_plots", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)

    chrom, start = read_windows(f"{a.ref_dir}/B73.tsv.gz")
    P = {os.path.basename(f)[:-4]: np.load(f, allow_pickle=True) for f in sorted(glob.glob(f"{a.hmm_dir}/*.npz"))}
    if not P: raise SystemExit(f"no .npz in {a.hmm_dir}")
    rho = P[next(iter(P))]["rho_grid"]
    info = {s: dict(role=str(z["role"]), group=str(z["taxon"]), sigma=float(z["sigma"])) for s, z in P.items()}
    for s in info: info[s]["flag"] = "NOISY" if info[s]["sigma"] > a.noisy_sigma else "ok"

    use = [s for s in P if info[s]["role"] == "donor" and info[s]["group"] != "Zm" and info[s]["flag"] == "ok"]
    L = np.sum([P[s]["loglik"] for s in use], axis=0); j = int(np.argmax(L))
    groups = sorted({info[s]["group"] for s in use})
    with open(f"{a.out_dir}/rho_loglik.tsv", "w") as f:
        f.write("RHO\tALL\t" + "\t".join(groups) + "\n")
        per = {g: np.sum([P[s]["loglik"] for s in use if info[s]["group"] == g], axis=0) for g in groups}
        for k, r in enumerate(rho):
            f.write(f"{r:g}\t{L[k]-L.max():.1f}\t" + "\t".join(f"{per[g][k]-per[g].max():.1f}" for g in groups) + "\n")
    print(f"rho fitted on {len(use)} donor pools (no Zm, no NOISY): best rho = {rho[j]:g}")
    for g in groups: print(f"  {g}: best rho alone = {rho[int(np.argmax(per[g]))]:g}")

    summ = open(f"{a.out_dir}/pool_summary.tsv", "w")
    summ.write("SAMPLE\tROLE\tGROUP\tSIGMA\tFLAG\tPARAM_KIND\tPARAM\tFRAC_B73\tFRAC_SEG\tFRAC_DONOR\tDONOR_SHARE\tDONOR_SHARE_R\tN_TRACTS\n")
    tr = gzip.open(f"{a.out_dir}/tracts_all.tsv.gz", "wt"); tr.write("SAMPLE\tGROUP\tCHROM\tSTART\tEND\tSTATE\n")
    blocks = {}
    for s, z in P.items():
        st, stR = z["states_R1"][j], z["states_R"][j]
        v, vR = st[st >= 0], stR[stR >= 0]
        fB, fS, fD = (np.mean(v == k) for k in range(3))
        sh, shR = fD + fS / 2, np.mean(vR == 2) + np.mean(vR == 1) / 2
        bl = runs(st, chrom, start, a.window); blocks[s] = bl
        nt = sum(1 for b in bl if b[3] > 0)
        i = info[s]
        summ.write(f"{s}\t{i['role']}\t{i['group']}\t{i['sigma']:.3f}\t{i['flag']}\t{z['param_kind']}\t{z['par_best'][j]:g}\t"
                   f"{fB:.4f}\t{fS:.4f}\t{fD:.4f}\t{sh:.4f}\t{shR:.4f}\t{nt}\n")
        for c, b0, b1, k in bl: tr.write(f"{s}\t{i['group']}\t{c}\t{b0}\t{b1}\t{STATES[k]}\n")
    summ.close(); tr.close()

    print("\nGROUP   N  NOISY  median share  mean share (ok pools)  median tracts")
    for g in sorted({i["group"] for i in info.values()}):
        ss = [s for s in P if info[s]["group"] == g]; ok = [s for s in ss if info[s]["flag"] == "ok"]
        shs = [np.mean(P[s]["states_R1"][j] == 2) / np.mean(P[s]["states_R1"][j] >= 0) +
               0.5 * np.mean(P[s]["states_R1"][j] == 1) / np.mean(P[s]["states_R1"][j] >= 0) for s in ok]
        nts = [sum(1 for b in blocks[s] if b[3] > 0) for s in ok]
        if ok: print(f"{g:5s} {len(ss):4d} {len(ss)-len(ok):6d} {100*np.median(shs):12.2f} {100*np.mean(shs):21.2f} {np.median(nts):14.0f}")
    print(f"\nWROTE: {a.out_dir}/rho_loglik.tsv, pool_summary.tsv, tracts_all.tsv.gz")

    if a.no_plots: return
    try:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        from matplotlib.patches import Patch
    except ImportError:
        print("matplotlib not available: no figures"); return
    col = {0: "#2850C0", 1: "#4F1F87", 2: "#C0202A"}
    chs = list(dict.fromkeys(chrom)); clen = {c: start[chrom == c].max() + a.window for c in chs}
    leg = [Patch(color=col[0], label="AA (B73/B73)"), Patch(color=col[1], label="AB (Het / SEG)"), Patch(color=col[2], label="BB (Teo/Teo)")]
    for g in sorted({i["group"] for i in info.values()}):
        with PdfPages(f"{a.out_dir}/blocks_{g}.pdf") as pdf:
            for s in sorted(x for x in P if info[x]["group"] == g):
                fig, ax = plt.subplots(figsize=(11, 6.5))
                for n, c in enumerate(chs):
                    y = len(chs) - 1 - n
                    ax.add_patch(plt.Rectangle((0, y - 0.35), clen[c] / 1e6, 0.7, color=col[0], lw=0))
                    for cc, b0, b1, k in blocks[s]:
                        if cc == c and k > 0: ax.add_patch(plt.Rectangle((b0 / 1e6, y - 0.35), (b1 - b0) / 1e6, 0.7, color=col[k], lw=0))
                st = P[s]["states_R1"][j]; v = st[st >= 0]
                flag = "  [NOISY]" if info[s]["flag"] == "NOISY" else ""
                ax.set_title(f"{s} ({g}) - k-mer HMM, donor share {100*(np.mean(v==2)+np.mean(v==1)/2):.1f}%{flag}", fontweight="bold")
                ax.set_yticks(range(len(chs))); ax.set_yticklabels(chs[::-1]); ax.set_xlim(0, max(clen.values()) / 1e6 + 5)
                ax.set_ylim(-0.8, len(chs) - 0.2); ax.set_xlabel("Position (Mb)"); ax.tick_params(axis="y", length=0)
                for sp in ("top", "right", "left"): ax.spines[sp].set_visible(False)
                ax.legend(handles=leg, loc="lower right", fontsize=9); fig.tight_layout(); pdf.savefig(fig); plt.close(fig)
    print(f"WROTE: {a.out_dir}/blocks_<GROUP>.pdf")


if __name__ == "__main__":
    main()
