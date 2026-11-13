#!/usr/bin/env python3
"""
14_hmm_pooled.py - Nirwan's BC2S3 HMM (old_scripts/HMM_introgression/hmm_viterbi_bc2s3.py)
adapted to pooled read counts. Same structure: per-marker HMM, NAM genetic-map
transitions (Haldane, scaled by rho), sticky/stationary transition matrix,
forward-backward posteriors, Viterbi or posterior + rigidity (hysteresis) decoding,
optional min-run cleanup.

What changed from the original:
  1. Emissions, still from the ANGSD GLs. Each sample is a POOL of plants from one BC2S3
     family, so the GLs are read as likelihoods of the reads at ALT allele fraction
     0 (GL_0), 1/2 (GL_1) and 1 (GL_2) rather than as genotypes of one diploid plant.
     The likelihood at any other fraction p comes from the curve through those three
     points (exact for sites with 1-2 reads, i.e. almost all sites at ~0.8x):
       L(p) = L0 (2p-1)(p-1) + L1 4p(1-p) + L2 p(2p-1)
     States are
       B73   family fixed B73            : L(0)
       SEG   family segregating          : c * mean_f L(f) + (1-c) * L(0),  f = 1/4, 1/2, 3/4
       DONOR family fixed donor          : c * L(1)        + (1-c) * L(0)
     c = chance the donor carries the ALT allele at the site = kappa x (teosinte panel AC/AN).
     kappa is learned per pool by EM (until it changes < --em_tol, max --em_rounds): donors that are under-represented in
     the panel (e.g. diploperennis) carry panel ALT alleles less often than the panel
     frequency says; without kappa their fixed-donor regions look "segregating".
     Each round: run forward-backward, then
       kappa = sum_t w_t r_t / sum_t w_t a_t   (+ a weak prior at kappa = 1)
     over SEG/DONOR markers, w_t = P(SEG or DONOR), r_t = posterior that the donor carries
     ALT at that marker, a_t = panel AC/AN.
     Sequencing error is already inside the GLs (base quality).
     The GL fixes (eta_hh_from_rh, eta_rr_from_rh, rh_penalty) are no longer needed.
  2. Each chromosome is its own chain (the original ran chr1..chr10 as one chain).
  3. thetas[t] is the recombination fraction between markers t-1 and t (the original
     list started with an extra 0.0, shifting every transition by one marker).
  4. Rigidity (hysteresis): the switch is placed at the first marker of the agreeing run,
     not the R-th, so boundaries are no longer shifted by R-1 markers.
  5. Tracts cover each chromosome end to end (boundaries at the midpoint between the
     last marker of one state and the first of the next), so lengths sum correctly.

Input : <in_dir>/<chr>/<sample>.tsv.gz   CHROM POS B73 ALT GL_0 GL_1 GL_2 PANEL_AC PANEL_AN
        (step 13; GLs are log10, already in B73 orientation)
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
# NEW: pool-state emissions from diploid GLs
# -----------------------------
LN10 = math.log(10.0)


class GLPoolEmission:
    def __init__(self, seg_freqs=(0.25, 0.5, 0.75), c_min=0.005, floor=1e-6):
        # Lagrange weights of L0, L1, L2 at each SEG fraction f
        self.w = [((2 * f - 1) * (f - 1), 4 * f * (1 - f), f * (2 * f - 1)) for f in seg_freqs]
        self.c_min = c_min
        self.floor = floor        # keeps L(f) > 0 (interpolation can dip below 0 at sites with 3+ reads)

    def prepare(self, gl0, gl1, gl2):
        """GLs (log10) -> (L0, Lseg, L2): relative likelihood at ALT fraction 0, SEG mean, 1."""
        m = max(gl0, gl1, gl2)
        L0, L1, L2 = (math.exp((g - m) * LN10) for g in (gl0, gl1, gl2))
        lf = [max(a * L0 + b * L1 + d * L2, self.floor) for a, b, d in self.w]
        return (max(L0, 1e-300), sum(lf) / len(lf), L2)

    def c_of(self, a, kappa, fc=None):
        """fc = family donor-allele probability (PALT, step 13c) if known, else panel freq x kappa"""
        c = fc if fc is not None else kappa * a
        return min(max(c, self.c_min), 1.0 - self.c_min)

    def __call__(self, L, c):
        L0, Ls, L2 = L
        return [math.log(L0), math.log(c * Ls + (1 - c) * L0), math.log(c * L2 + (1 - c) * L0)]


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
# Forward-backward (unchanged except it also returns the log-likelihood;
#                   thetas[t] = transition t-1 -> t)
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
    return post, sum(c)


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
# RTIGER-ish rigidity via hysteresis (called per chromosome). Change from the original:
# when R markers in a row agree on a new state, the switch is placed at the first marker
# of that run (the original placed it at the R-th, so every boundary sat R-1 markers late)
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
            # boundary fix: the switch starts at the FIRST marker of the run, not the R-th
            for k in range(t - run_len + 1, t):
                out[k] = cur
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
    ap.add_argument("--donor_dir", default=None, help="step-13c folder: <donor_dir>/<family>/<chr>.tsv.gz (family-informed run)")
    ap.add_argument("--family", default=None, help="this pool's BC1 family, e.g. Zd.0010_P2")

    ap.add_argument("--prior_b73", type=float, default=0.844)
    ap.add_argument("--prior_seg", type=float, default=0.0625)
    ap.add_argument("--prior_donor", type=float, default=0.094)

    ap.add_argument("--min_morgan", type=float, default=1e-8)
    ap.add_argument("--rho", type=float, default=0.2)

    ap.add_argument("--em_rounds", type=int, default=30, help="max EM rounds to learn kappa (donor match to the panel)")
    ap.add_argument("--em_tol", type=float, default=0.005, help="stop EM when kappa changes less than this")
    ap.add_argument("--em_starts", type=float, nargs="+", default=[1.0, 0.5, 0.25], help="EM starting values for kappa")
    ap.add_argument("--decode", choices=["posterior_hysteresis", "posterior", "viterbi"], default="posterior_hysteresis")
    ap.add_argument("--rigidity", "-R", type=int, default=50, help="R consecutive markers before switching (posterior_hysteresis)")
    ap.add_argument("--min_run_donor", type=int, default=1)
    ap.add_argument("--min_run_seg", type=int, default=1)
    args = ap.parse_args()

    s = args.prior_b73 + args.prior_seg + args.prior_donor
    pi = [args.prior_b73 / s, args.prior_seg / s, args.prior_donor / s]
    emis = GLPoolEmission()
    chroms = args.chroms.split(",")
    chr_len = {}
    with open(args.fai) as f:
        for line in f:
            p = line.split("\t")
            chr_len[p[0]] = int(p[1])

    sp = gzip.open(args.out_prefix + ".statepath.tsv.gz", "wt")
    sp.write("CHROM\tPOS\tB73\tALT\tGL_0\tGL_1\tGL_2\tPANEL_AC\tPANEL_AN\tSTATE\tP_B73\tP_SEG\tP_DONOR\tTHETA\tC\tC_SOURCE\n")
    tr = gzip.open(args.out_prefix + ".tracts.bed.gz", "wt")
    tr.write("chrom\tstart\tend\tstate\n")
    tot = {"B73": 0, "SEG": 0, "DONOR": 0}
    n_markers = 0
    n_tracts = 0

    # ---- load every chromosome (markers, GL-derived likelihoods, panel freq, thetas) ----
    data = {}
    n_family = 0
    for chrom in chroms:
        fn = os.path.join(args.in_dir, chrom, args.sample + ".tsv.gz")
        rows = []
        if os.path.exists(fn):
            with gzip.open(fn, "rt") as f:
                for line in f:
                    p = line.rstrip("\n").split("\t")
                    if len(p) < 9:
                        continue
                    rows.append((int(p[1]), p[2], p[3], float(p[4]), float(p[5]), float(p[6]), int(p[7]), int(p[8])))
        rows.sort()
        # family donor alleles (family-informed run): PALT where REF/ALT match this marker
        fam = {}
        if args.donor_dir and args.family:
            ff = os.path.join(args.donor_dir, args.family, chrom + ".tsv.gz")
            if os.path.exists(ff):
                want = {r[0] for r in rows}
                with gzip.open(ff, "rt") as f:
                    for line in f:
                        p = line.split("\t")
                        pos = int(p[0])
                        if pos in want:
                            fam[pos] = (p[1], p[2], float(p[3]))
        kept, fcs = [], []
        for r in rows:
            pos, b73, alt, ac, an = r[0], r[1], r[2], r[6], r[7]
            fx = fam.get(pos)
            fc = fx[2] if (fx is not None and fx[0] == b73 and fx[1] == alt) else None
            if fc is None and an <= 0:
                continue                                  # non-panel site without family info: no use
            kept.append(r)
            fcs.append(fc)
        rows = kept
        if not rows:                                      # no data: leave chromosome uncalled
            continue
        n_family += sum(fc is not None for fc in fcs)

        bp_map, cm_map = read_map(os.path.join(args.map_dir, chrom + ".map"))
        Ls, af, thetas = [], [], []
        prev_cm = None
        for pos, b73, alt, g0, g1, g2, ac, an in rows:
            Ls.append(emis.prepare(g0, g1, g2))
            af.append(ac / an if an > 0 else 0.5)
            cm_here = interp_cm(bp_map, cm_map, pos)
            if prev_cm is None:
                thetas.append(0.0)                         # first marker: no transition into it
            else:
                d_m = max(abs(cm_here - prev_cm) / 100.0, args.min_morgan)
                thetas.append(max(haldane_theta(d_m), 1e-8))
            prev_cm = cm_here
        data[chrom] = (rows, Ls, af, thetas, fcs)

    # ---- EM for kappa (how often this pool's donor carries the panel's ALT alleles) ----
    # SEG at kappa k and DONOR at kappa k/2 look alike at 1 read per site, so EM can stop
    # at the wrong one. Run it from several starts and keep the best log-likelihood
    # (sites with 2+ reads tell the two apart).
    PRIOR_W = 200.0                                       # weak prior: 200 markers' worth at kappa = 1

    def run_all(kap):
        out, ll = {}, 0.0
        for chrom, (rows, Ls, af, thetas, fcs) in data.items():
            cs = [emis.c_of(a, kap, fc) for a, fc in zip(af, fcs)]
            post, l = forward_backward([emis(L, c) for L, c in zip(Ls, cs)], thetas, pi, args.rho)
            out[chrom] = (cs, post); ll += l
        return out, ll

    def em(k0):
        kap, used = k0, 0
        for rnd in range(args.em_rounds):
            res, ll = run_all(kap)
            num, den = 0.0, 0.0
            for chrom, (cs, post) in res.items():
                Ls, af, fcs = data[chrom][1], data[chrom][2], data[chrom][4]
                for (L0, Lseg, L2), c, a, fc, (pB, pS, pD) in zip(Ls, cs, af, fcs, post):
                    if fc is not None:
                        continue                          # family-known sites do not inform kappa
                    num += pS * c * Lseg / (c * Lseg + (1 - c) * L0) + pD * c * L2 / (c * L2 + (1 - c) * L0)
                    den += (pS + pD) * a
            new = min(max((num + PRIOR_W * 0.3) / (den + PRIOR_W * 0.3), 0.02), 1.5)
            used = rnd + 1
            if abs(new - kap) < args.em_tol:
                kap = new
                break
            kap = new
        return kap, run_all(kap)[1], used

    best = None
    for k0 in args.em_starts:
        k, ll, used = em(k0)
        if best is None or ll > best[1]:
            best = (k, ll, used, k0)
    kappa, _, em_used, kappa_start = best

    # ---- final run, decode, write ----
    for chrom, (rows, Ls, af, thetas, fcs) in data.items():
        L = chr_len[chrom]
        cs = [emis.c_of(a, kappa, fc) for a, fc in zip(af, fcs)]
        obs_ll = [emis(Lx, c) for Lx, c in zip(Ls, cs)]
        post, _ = forward_backward(obs_ll, thetas, pi, args.rho)
        if args.decode == "viterbi":
            path = viterbi_standard(obs_ll, thetas, pi, args.rho)
        else:
            best = [max(range(3), key=lambda k: post[t][k]) for t in range(len(post))]
            path = decode_hysteresis(best, args.rigidity) if args.decode == "posterior_hysteresis" else best
        path = run_filter(path, args.min_run_donor, args.min_run_seg)

        for i, (pos, b73, alt, g0, g1, g2, ac, an) in enumerate(rows):
            pB, pS, pD = post[i]
            sp.write(f"{chrom}\t{pos}\t{b73}\t{alt}\t{g0:g}\t{g1:g}\t{g2:g}\t{ac}\t{an}\t{STATES[path[i]]}\t{pB:.4g}\t{pS:.4g}\t{pD:.4g}\t{thetas[i]:.4g}\t{cs[i]:.4g}\t{'family' if fcs[i] is not None else 'panel'}\n")
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
        f.write("SAMPLE\tMARKERS\tCALLED_BP\tFRAC_B73\tFRAC_SEG\tFRAC_DONOR\tDONOR_SHARE\tN_TRACTS\tKAPPA\tEM_ROUNDS\tKAPPA_START\tFAMILY\tFAMILY_MARKERS\n")
        if g:
            f.write(f"{args.sample}\t{n_markers}\t{g}\t{tot['B73']/g:.4f}\t{tot['SEG']/g:.4f}\t{tot['DONOR']/g:.4f}\t"
                    f"{(tot['DONOR'] + 0.5 * tot['SEG'])/g:.4f}\t{n_tracts}\t{kappa:.3f}\t{em_used}\t{kappa_start:g}\t{args.family or "NA"}\t{n_family}\n")
        else:
            f.write(f"{args.sample}\t0\t0\tNA\tNA\tNA\tNA\t0\tNA\t0\tNA\t{args.family or "NA"}\t0\n")


if __name__ == "__main__":
    main()
