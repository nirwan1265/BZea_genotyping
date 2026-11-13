#!/usr/bin/env python3
"""
14_hmm_pooled.py - Nirwan's BC2S3 HMM (old_scripts/HMM_introgression/hmm_viterbi_bc2s3.py)
adapted to pooled read counts. Same structure: per-marker HMM, NAM genetic-map
transitions (Haldane, scaled by rho), sticky/stationary transition matrix,
forward-backward posteriors, Viterbi or posterior + rigidity (hysteresis) decoding,
optional min-run cleanup.

What changed from the original:
  1. Emissions: each sample is a POOL of plants from one BC2S3 family, so the diploid
     genotype likelihoods (RR/RH/HH) are replaced by read counts. States are now
       B73   family fixed B73
       SEG   family segregating (some plants carry the donor)
       DONOR family fixed donor
     At a marker with n reads, k of them ALT, sequencing/background error e, and
     c = chance the donor carries the ALT allele there (teosinte panel AC/AN):
       B73   : Bin(k; n, e)
       DONOR : c * Bin(k; n, 1-e)              + (1-c) * Bin(k; n, e)
       SEG   : c * mean_f Bin(k; n, f(1-e)+(1-f)e) + (1-c) * Bin(k; n, e)
               f = share of donor copies in the pool, averaged over 1/4, 1/2, 3/4
     (the binomial coefficient is the same for every state and is left out)
     The GL fixes (eta_hh_from_rh, eta_rr_from_rh, rh_penalty) are no longer needed.
  2. Each chromosome is its own chain (the original ran chr1..chr10 as one chain).
  3. thetas[t] is the recombination fraction between markers t-1 and t (the original
     list started with an extra 0.0, shifting every transition by one marker).
  4. Tracts cover each chromosome end to end (boundaries at the midpoint between the
     last marker of one state and the first of the next), so lengths sum correctly.

Input : <in_dir>/<chr>/<sample>.tsv.gz   CHROM POS B73 ALT REF_READS ALT_READS PANEL_AC PANEL_AN
Output: <out_prefix>.statepath.tsv.gz, <out_prefix>.tracts.bed.gz, <out_prefix>.summary.tsv
"""
import argparse
import gzip
import math
import os
from bisect import bisect_right

STATES = ["B73", "SEG", "DONOR"]


# -----------------------------
# Genetic map utilities (unchanged)
# -----------------------------
def read_map(map_path: str):
    bp, cm = [], []
    with open(map_path, "r") as f:
        for line in f:
            line = line.strip()
            if (not line) or line.startswith("#") or line.startswith("locus"):
                continue
            p = line.split()
            if len(p) < 2:
                continue
            try:
                b = int(float(p[0]))
                c = float(p[-1])
            except ValueError:
                continue
            bp.append(b)
            cm.append(c)
    if len(bp) < 2:
        raise RuntimeError(f"Map too small: {map_path}")
    z = sorted(zip(bp, cm))
    bp = [x for x, _ in z]
    cm = [y for _, y in z]
    return bp, cm


def interp_cm(bp_map, cm_map, pos: int) -> float:
    i = bisect_right(bp_map, pos) - 1
    if i < 0:
        x0, x1 = bp_map[0], bp_map[1]
        y0, y1 = cm_map[0], cm_map[1]
        return y0 if x1 == x0 else y0 + (pos - x0) * (y1 - y0) / (x1 - x0)
    if i >= len(bp_map) - 1:
        x0, x1 = bp_map[-2], bp_map[-1]
        y0, y1 = cm_map[-2], cm_map[-1]
        return y1 if x1 == x0 else y0 + (pos - x0) * (y1 - y0) / (x1 - x0)
    x0, x1 = bp_map[i], bp_map[i + 1]
    y0, y1 = cm_map[i], cm_map[i + 1]
    return y0 if x1 == x0 else y0 + (pos - x0) * (y1 - y0) / (x1 - x0)


def haldane_theta(d_morgan: float) -> float:
    return 0.5 * (1.0 - math.exp(-2.0 * d_morgan))


# -----------------------------
# Log-space helpers (unchanged)
# -----------------------------
def logsumexp3(a, b, c):
    m = max(a, b, c)
    return m + math.log(math.exp(a - m) + math.exp(b - m) + math.exp(c - m))


def slog(x):
    return -1e300 if x <= 0 else math.log(x)


def log_mix(c, la, lb):
    # log(c*exp(la) + (1-c)*exp(lb))
    m = max(la, lb)
    return m + math.log(c * math.exp(la - m) + (1.0 - c) * math.exp(lb - m))


# -----------------------------
# NEW: pooled read-count emissions
# -----------------------------
class PoolEmission:
    def __init__(self, err, seg_freqs=(0.25, 0.5, 0.75), c_min=0.01):
        self.le, self.l1e = math.log(err), math.log(1.0 - err)
        self.seg = [(math.log(f * (1 - err) + (1 - f) * err), math.log(1 - (f * (1 - err) + (1 - f) * err)))
                    for f in seg_freqs]
        self.lnseg = math.log(len(seg_freqs))
        self.c_min = c_min

    def __call__(self, n, k, ac, an):
        c = ac / an if an > 0 else 0.5
        c = min(max(c, self.c_min), 1.0 - self.c_min)
        l_b73 = k * self.le + (n - k) * self.l1e                     # read ALT only by error
        l_alt = k * self.l1e + (n - k) * self.le                     # donor carries ALT, fixed
        segs = [k * la + (n - k) * lr for la, lr in self.seg]        # donor carries ALT, pool freq f
        m = max(segs)
        l_seg_alt = m + math.log(sum(math.exp(x - m) for x in segs)) - self.lnseg
        return [l_b73, log_mix(c, l_seg_alt, l_b73), log_mix(c, l_alt, l_b73)]


# -----------------------------
# HMM transitions (unchanged: sticky + stationary)
# -----------------------------
def build_A(theta, pi, rho):
    s = rho * theta
    s = max(1e-12, min(0.25, s))
    A = [[0.0] * 3 for _ in range(3)]
    for i in range(3):
        A[i][i] = 1.0 - s
        denom = 1.0 - pi[i]
        for j in range(3):
            if j == i:
                continue
            A[i][j] = s * (pi[j] / denom)
    return A


# -----------------------------
# Forward-backward (unchanged; thetas[t] = transition t-1 -> t)
# -----------------------------
def forward_backward(obs_ll, thetas, pi, rho):
    n = len(obs_ll)
    alpha = [[-1e300] * 3 for _ in range(n)]
    c = [0.0] * n

    for s in range(3):
        alpha[0][s] = slog(pi[s]) + obs_ll[0][s]
    c[0] = logsumexp3(alpha[0][0], alpha[0][1], alpha[0][2])
    for s in range(3):
        alpha[0][s] -= c[0]

    for t in range(1, n):
        A = build_A(thetas[t], pi, rho)
        for j in range(3):
            v0 = alpha[t - 1][0] + slog(A[0][j])
            v1 = alpha[t - 1][1] + slog(A[1][j])
            v2 = alpha[t - 1][2] + slog(A[2][j])
            alpha[t][j] = logsumexp3(v0, v1, v2) + obs_ll[t][j]
        c[t] = logsumexp3(alpha[t][0], alpha[t][1], alpha[t][2])
        for j in range(3):
            alpha[t][j] -= c[t]

    beta = [[0.0] * 3 for _ in range(n)]
    for t in range(n - 2, -1, -1):
        A = build_A(thetas[t + 1], pi, rho)
        for i in range(3):
            v0 = slog(A[i][0]) + obs_ll[t + 1][0] + beta[t + 1][0]
            v1 = slog(A[i][1]) + obs_ll[t + 1][1] + beta[t + 1][1]
            v2 = slog(A[i][2]) + obs_ll[t + 1][2] + beta[t + 1][2]
            beta[t][i] = logsumexp3(v0, v1, v2) - c[t + 1]

    post = [[0.0] * 3 for _ in range(n)]
    for t in range(n):
        v0 = alpha[t][0] + beta[t][0]
        v1 = alpha[t][1] + beta[t][1]
        v2 = alpha[t][2] + beta[t][2]
        z = logsumexp3(v0, v1, v2)
        post[t][0] = math.exp(v0 - z)
        post[t][1] = math.exp(v1 - z)
        post[t][2] = math.exp(v2 - z)
    return post


# -----------------------------
# Standard Viterbi (unchanged)
# -----------------------------
def viterbi_standard(obs_ll, thetas, pi, rho):
    n = len(obs_ll)
    dp = [[-1e300] * 3 for _ in range(n)]
    ptr = [[0] * 3 for _ in range(n)]

    for s in range(3):
        dp[0][s] = slog(pi[s]) + obs_ll[0][s]

    for t in range(1, n):
        A = build_A(thetas[t], pi, rho)
        for j in range(3):
            best = -1e300
            best_i = 0
            for i in range(3):
                val = dp[t - 1][i] + slog(A[i][j])
                if val > best:
                    best = val
                    best_i = i
            dp[t][j] = best + obs_ll[t][j]
            ptr[t][j] = best_i

    last = max(range(3), key=lambda s: dp[n - 1][s])
    path = [last]
    for t in range(n - 1, 0, -1):
        last = ptr[t][last]
        path.append(last)
    path.reverse()
    return path


# -----------------------------
# RTIGER-ish rigidity via hysteresis (unchanged; called per chromosome)
# -----------------------------
def decode_hysteresis(best_state, R):
    n = len(best_state)
    if n == 0:
        return []
    if R <= 1:
        return best_state[:]
    out = [best_state[0]]
    cur = best_state[0]
    run_state = None
    run_len = 0
    for t in range(1, n):
        s = best_state[t]
        if s == cur:
            run_state = None
            run_len = 0
            out.append(cur)
            continue
        if run_state is None or s != run_state:
            run_state = s
            run_len = 1
        else:
            run_len += 1
        if run_len >= R:
            cur = run_state
            run_state = None
            run_len = 0
        out.append(cur)
    return out


# -----------------------------
# Post-run min-run filter (unchanged logic; per chromosome)
# -----------------------------
def run_filter(path, min_run_donor, min_run_seg):
    out = path[:]
    n = len(out)
    i = 0
    while i < n:
        j = i
        while j + 1 < n and out[j + 1] == out[i]:
            j += 1
        run_len = j - i + 1
        st = out[i]
        if st == 2 and run_len < min_run_donor:  # DONOR -> SEG
            for k in range(i, j + 1):
                out[k] = 1
        if st == 1 and run_len < min_run_seg:    # SEG -> B73
            for k in range(i, j + 1):
                out[k] = 0
        i = j + 1
    return out


# -----------------------------
# Main
# -----------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--in_dir", required=True, help="step-13 folder: <in_dir>/<chr>/<sample>.tsv.gz")
    ap.add_argument("--map_dir", required=True, help="NAM genetic maps: <map_dir>/<chr>.map")
    ap.add_argument("--fai", required=True, help="reference .fai (chromosome lengths for tracts)")
    ap.add_argument("--out_prefix", required=True)
    ap.add_argument("--chroms", default=",".join(f"chr{i}" for i in range(1, 11)))

    ap.add_argument("--err", type=float, default=0.002, help="background ALT rate (B73 controls, step 12)")
    ap.add_argument("--prior_b73", type=float, default=0.844)
    ap.add_argument("--prior_seg", type=float, default=0.0625)
    ap.add_argument("--prior_donor", type=float, default=0.094)

    ap.add_argument("--min_morgan", type=float, default=1e-8)
    ap.add_argument("--rho", type=float, default=0.2)

    ap.add_argument("--decode", choices=["posterior_hysteresis", "posterior", "viterbi"], default="posterior_hysteresis")
    ap.add_argument("--rigidity", "-R", type=int, default=50, help="R consecutive markers before switching (posterior_hysteresis)")
    ap.add_argument("--min_run_donor", type=int, default=1)
    ap.add_argument("--min_run_seg", type=int, default=1)
    args = ap.parse_args()

    s = args.prior_b73 + args.prior_seg + args.prior_donor
    pi = [args.prior_b73 / s, args.prior_seg / s, args.prior_donor / s]
    emis = PoolEmission(args.err)
    chroms = args.chroms.split(",")
    chr_len = {}
    with open(args.fai) as f:
        for line in f:
            p = line.split("\t")
            chr_len[p[0]] = int(p[1])

    sp = gzip.open(args.out_prefix + ".statepath.tsv.gz", "wt")
    sp.write("CHROM\tPOS\tB73\tALT\tREF_READS\tALT_READS\tPANEL_AC\tPANEL_AN\tSTATE\tP_B73\tP_SEG\tP_DONOR\tTHETA\n")
    tr = gzip.open(args.out_prefix + ".tracts.bed.gz", "wt")
    tr.write("chrom\tstart\tend\tstate\n")
    tot = {"B73": 0, "SEG": 0, "DONOR": 0}
    n_markers = 0
    n_tracts = 0

    for chrom in chroms:
        fn = os.path.join(args.in_dir, chrom, args.sample + ".tsv.gz")
        rows = []
        if os.path.exists(fn):
            with gzip.open(fn, "rt") as f:
                for line in f:
                    p = line.rstrip("\n").split("\t")
                    if len(p) < 8:
                        continue
                    rows.append((int(p[1]), p[2], p[3], int(p[4]), int(p[5]), int(p[6]), int(p[7])))
        rows.sort()
        L = chr_len[chrom]
        if not rows:                                      # no reads: leave chromosome uncalled
            continue

        bp_map, cm_map = read_map(os.path.join(args.map_dir, chrom + ".map"))
        obs_ll, thetas = [], []
        prev_cm = None
        for pos, b73, alt, r, a, ac, an in rows:
            obs_ll.append(emis(r + a, a, ac, an))
            cm_here = interp_cm(bp_map, cm_map, pos)
            if prev_cm is None:
                thetas.append(0.0)                         # first marker: no transition into it
            else:
                d_m = max(abs(cm_here - prev_cm) / 100.0, args.min_morgan)
                thetas.append(max(haldane_theta(d_m), 1e-8))
            prev_cm = cm_here

        post = forward_backward(obs_ll, thetas, pi, args.rho)
        if args.decode == "viterbi":
            path = viterbi_standard(obs_ll, thetas, pi, args.rho)
        else:
            best = [max(range(3), key=lambda k: post[t][k]) for t in range(len(post))]
            path = decode_hysteresis(best, args.rigidity) if args.decode == "posterior_hysteresis" else best
        path = run_filter(path, args.min_run_donor, args.min_run_seg)

        for i, (pos, b73, alt, r, a, ac, an) in enumerate(rows):
            pB, pS, pD = post[i]
            sp.write(f"{chrom}\t{pos}\t{b73}\t{alt}\t{r}\t{a}\t{ac}\t{an}\t{STATES[path[i]]}\t{pB:.4g}\t{pS:.4g}\t{pD:.4g}\t{thetas[i]:.4g}\n")
        n_markers += len(rows)

        # tracts cover the whole chromosome: boundary = midpoint between neighbouring markers
        start = 0
        for i in range(1, len(rows) + 1):
            if i == len(rows) or path[i] != path[i - 1]:
                end = L if i == len(rows) else (rows[i - 1][0] + rows[i][0]) // 2
                st = STATES[path[i - 1]]
                tr.write(f"{chrom}\t{start}\t{end}\t{st}\n")
                tot[st] += end - start
                n_tracts += 1
                start = end

    sp.close()
    tr.close()
    g = sum(tot.values())
    with open(args.out_prefix + ".summary.tsv", "w") as f:
        f.write("SAMPLE\tMARKERS\tCALLED_BP\tFRAC_B73\tFRAC_SEG\tFRAC_DONOR\tDONOR_SHARE\tN_TRACTS\n")
        if g:
            f.write(f"{args.sample}\t{n_markers}\t{g}\t{tot['B73']/g:.4f}\t{tot['SEG']/g:.4f}\t{tot['DONOR']/g:.4f}\t"
                    f"{(tot['DONOR'] + 0.5 * tot['SEG'])/g:.4f}\t{n_tracts}\n")
        else:
            f.write(f"{args.sample}\t0\t0\tNA\tNA\tNA\tNA\t0\n")


if __name__ == "__main__":
    main()
