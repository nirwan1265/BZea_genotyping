#!/usr/bin/env python3
"""
HMM Viterbi with RTIGER-style Rigidity

This version implements rigidity (R) INSIDE the Viterbi algorithm:
- Requires R consecutive markers to support a state change
- Prevents noisy single-SNP state switches
- More biologically realistic for BC2S3 (few recombination events)

Key difference from standard Viterbi:
- Standard: Can switch state at any position
- Rigidity: Must have R consecutive markers supporting new state before switching
"""

import argparse
import gzip
import math
from bisect import bisect_right

STATES = ["RR", "RH", "HH"]


def read_map(map_path: str):
    bp, cm = [], []
    with open(map_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("locus"):
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


def logsumexp3(a, b, c):
    m = max(a, b, c)
    return m + math.log(math.exp(a - m) + math.exp(b - m) + math.exp(c - m))


def mix_log(a, b, eta):
    if eta <= 0:
        return a
    if eta >= 1:
        return b
    m = max(a, b)
    return m + math.log((1 - eta) * math.exp(a - m) + eta * math.exp(b - m))


def apply_emission_adjustments(lRR, lRH, lHH, eta_hh_from_rh, eta_rr_from_rh, rh_penalty):
    lHH2 = mix_log(lHH, lRH, eta_hh_from_rh)
    lRR2 = mix_log(lRR, lRH, eta_rr_from_rh)
    lRH2 = lRH - rh_penalty
    return lRR2, lRH2, lHH2


def build_A(theta, pi, rho):
    s = rho * theta
    s = max(1e-12, min(0.25, s))
    A = [[0.0] * 3 for _ in range(3)]
    for i in range(3):
        A[i][i] = 1.0 - s
        denom = 1.0 - pi[i]
        for j in range(3):
            if j != i:
                A[i][j] = s * (pi[j] / denom)
    return A


def slog(x):
    return -1e300 if x <= 0 else math.log(x)


# ============================================================================
# RIGIDITY VITERBI - RTIGER-style
# ============================================================================
def viterbi_rigidity(obs_ll, thetas, pi, rho, rigidity, chroms):
    """
    Viterbi algorithm with rigidity constraint.

    Rigidity (R): Minimum number of consecutive markers required to switch states.

    At each position t, we compare:
    1. STAY path: Continue in same state (standard 1-step transition)
    2. SWITCH path: Switch to new state, but only if supported by R markers

    The switch path requires that the cumulative emission likelihood
    over the last R markers favors the new state.
    """
    n = len(obs_ll)
    R = rigidity

    if n == 0:
        return []

    if R <= 1:
        # Fall back to standard Viterbi if R=1
        return viterbi_standard(obs_ll, thetas, pi, rho)

    # Precompute cumulative emission log-likelihoods
    # cumsum_ll[t][s] = sum of obs_ll[0:t+1][s]
    cumsum_ll = [[0.0] * 3 for _ in range(n)]
    for s in range(3):
        cumsum_ll[0][s] = obs_ll[0][s]
    for t in range(1, n):
        for s in range(3):
            # Reset at chromosome boundaries
            if chroms[t] != chroms[t-1]:
                cumsum_ll[t][s] = obs_ll[t][s]
            else:
                cumsum_ll[t][s] = cumsum_ll[t-1][s] + obs_ll[t][s]

    # Get cumulative emission for window [t-R+1, t] (R markers)
    def get_window_ll(t, s):
        if t < R - 1:
            return cumsum_ll[t][s]
        # Check chromosome boundary
        start_t = t - R + 1
        if chroms[start_t] != chroms[t]:
            # Find where this chromosome starts
            for k in range(t, start_t - 1, -1):
                if chroms[k] != chroms[t]:
                    return cumsum_ll[t][s] - cumsum_ll[k][s]
            return cumsum_ll[t][s]
        return cumsum_ll[t][s] - (cumsum_ll[t-R][s] if t >= R else 0)

    # DP with rigidity
    # dp[t][s] = best log-probability to be in state s at position t
    dp = [[-1e300] * 3 for _ in range(n)]
    ptr = [[(-1, -1)] * 3 for _ in range(n)]  # (prev_state, switch_position)

    # Initialize first R positions with standard transitions
    for s in range(3):
        dp[0][s] = slog(pi[s]) + obs_ll[0][s]
        ptr[0][s] = (-1, 0)

    for t in range(1, n):
        A = build_A(thetas[t], pi, rho)
        new_chr = (chroms[t] != chroms[t-1])

        for j in range(3):
            # Option 1: STAY in same state (or switch from adjacent state)
            # Standard 1-step transition
            best_stay = -1e300
            best_stay_from = 0
            for i in range(3):
                val = dp[t-1][i] + slog(A[i][j]) + obs_ll[t][j]
                if val > best_stay:
                    best_stay = val
                    best_stay_from = i

            # Option 2: SWITCH with rigidity support
            # Only consider if we have R markers of evidence
            best_switch = -1e300
            best_switch_from = 0

            if t >= R and not new_chr:
                # Check if last R markers support state j
                window_ll_j = get_window_ll(t, j)

                for i in range(3):
                    if i == j:
                        continue

                    # Compare: staying in i vs switching to j
                    window_ll_i = get_window_ll(t, i)

                    # Only switch if j is clearly better over R markers
                    if window_ll_j > window_ll_i:
                        # Transition penalty for switching
                        trans_penalty = slog(A[i][j]) * R
                        val = dp[t-R][i] + window_ll_j + trans_penalty
                        if val > best_switch:
                            best_switch = val
                            best_switch_from = i

            # Take the better of stay vs switch
            if best_stay >= best_switch:
                dp[t][j] = best_stay
                ptr[t][j] = (best_stay_from, t)
            else:
                dp[t][j] = best_switch
                ptr[t][j] = (best_switch_from, t - R)

    # Backtrace
    last = max(range(3), key=lambda s: dp[n-1][s])
    path = [0] * n
    path[n-1] = last

    t = n - 1
    while t > 0:
        prev_state, switch_pos = ptr[t][path[t]]
        if prev_state == -1:
            break
        # Fill in the path from switch_pos to t with current state
        for k in range(switch_pos, t):
            if k >= 0:
                path[k] = path[t]
        if switch_pos > 0:
            path[switch_pos - 1] = prev_state
        t = switch_pos - 1

    # Fill any remaining positions
    current_state = path[0] if path[0] != 0 or dp[0][0] > max(dp[0][1], dp[0][2]) else max(range(3), key=lambda s: dp[0][s])
    for t in range(n):
        if path[t] == 0 and dp[t][0] < max(dp[t][1], dp[t][2]):
            path[t] = max(range(3), key=lambda s: dp[t][s])

    return path


def viterbi_standard(obs_ll, thetas, pi, rho):
    """Standard Viterbi without rigidity."""
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
                val = dp[t-1][i] + slog(A[i][j])
                if val > best:
                    best = val
                    best_i = i
            dp[t][j] = best + obs_ll[t][j]
            ptr[t][j] = best_i

    last = max(range(3), key=lambda s: dp[n-1][s])
    path = [last]
    for t in range(n-1, 0, -1):
        last = ptr[t][last]
        path.append(last)
    path.reverse()
    return path


# ============================================================================
# SIMPLIFIED RIGIDITY: Majority voting over R markers
# ============================================================================
def viterbi_majority_rigidity(obs_ll, chroms, rigidity):
    """
    Simpler rigidity implementation using majority voting.

    For each position, look at R markers centered on it.
    Assign the state that has the highest total probability across the window.
    """
    n = len(obs_ll)
    R = rigidity
    half_R = R // 2

    path = [0] * n

    for t in range(n):
        # Define window [start, end]
        start = max(0, t - half_R)
        end = min(n - 1, t + half_R)

        # Adjust for chromosome boundaries
        while start < t and chroms[start] != chroms[t]:
            start += 1
        while end > t and chroms[end] != chroms[t]:
            end -= 1

        # Sum log-likelihoods across window for each state
        scores = [0.0, 0.0, 0.0]
        for k in range(start, end + 1):
            for s in range(3):
                scores[s] += obs_ll[k][s]

        # Pick best state
        path[t] = max(range(3), key=lambda s: scores[s])

    return path


def forward_backward(obs_ll, thetas, pi, rho):
    """Standard forward-backward for posterior probabilities."""
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
            v0 = alpha[t-1][0] + slog(A[0][j])
            v1 = alpha[t-1][1] + slog(A[1][j])
            v2 = alpha[t-1][2] + slog(A[2][j])
            alpha[t][j] = logsumexp3(v0, v1, v2) + obs_ll[t][j]
        c[t] = logsumexp3(alpha[t][0], alpha[t][1], alpha[t][2])
        for j in range(3):
            alpha[t][j] -= c[t]

    beta = [[0.0] * 3 for _ in range(n)]
    for t in range(n-2, -1, -1):
        A = build_A(thetas[t+1], pi, rho)
        for i in range(3):
            v0 = slog(A[i][0]) + obs_ll[t+1][0] + beta[t+1][0]
            v1 = slog(A[i][1]) + obs_ll[t+1][1] + beta[t+1][1]
            v2 = slog(A[i][2]) + obs_ll[t+1][2] + beta[t+1][2]
            beta[t][i] = logsumexp3(v0, v1, v2) - c[t+1]

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


def main():
    ap = argparse.ArgumentParser(description="HMM Viterbi with RTIGER-style Rigidity")
    ap.add_argument("--in_gl_tsv_gz", required=True, help="Input GL file")
    ap.add_argument("--map_dir", required=True, help="Directory with genetic maps")
    ap.add_argument("--out_statepath_gz", required=True, help="Output state path")
    ap.add_argument("--out_tracts_bed_gz", required=True, help="Output tracts BED")

    # Priors
    ap.add_argument("--prior_rr", type=float, default=0.85)
    ap.add_argument("--prior_rh", type=float, default=0.05)
    ap.add_argument("--prior_hh", type=float, default=0.10)

    # HMM parameters
    ap.add_argument("--min_morgan", type=float, default=1e-8)
    ap.add_argument("--rho", type=float, default=2.0,
                    help="Transition stickiness (lower=stickier)")

    # Emission adjustments
    ap.add_argument("--eta_hh_from_rh", type=float, default=0.3,
                    help="HH rescue from RH (0-0.5)")
    ap.add_argument("--eta_rr_from_rh", type=float, default=0.0)
    ap.add_argument("--rh_penalty", type=float, default=0.5,
                    help="RH penalty in ln units")

    # RIGIDITY - the key parameter!
    ap.add_argument("--rigidity", "-R", type=int, default=50,
                    help="Rigidity: min consecutive markers to switch states (RTIGER-style). "
                         "Higher = fewer state changes. Try 20-100 for BC2S3.")

    # Method
    ap.add_argument("--method", choices=["rigidity", "majority", "standard"],
                    default="majority",
                    help="Viterbi method: rigidity (RTIGER-style), majority (simpler), standard")

    args = ap.parse_args()

    print(f"=== HMM with Rigidity ===")
    print(f"Rigidity (R): {args.rigidity}")
    print(f"Method: {args.method}")
    print(f"Priors: RR={args.prior_rr}, RH={args.prior_rh}, HH={args.prior_hh}")
    print()

    s = args.prior_rr + args.prior_rh + args.prior_hh
    pi = [args.prior_rr / s, args.prior_rh / s, args.prior_hh / s]

    # Load maps
    maps = {}
    for c in range(1, 11):
        chrom = f"chr{c}"
        maps[chrom] = read_map(f"{args.map_dir}/{chrom}.map")

    # Read GL data
    rows = []
    with gzip.open(args.in_gl_tsv_gz, "rt") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            p = line.split("\t")
            if len(p) < 7:
                continue

            chrom = p[0]
            if chrom not in maps:
                continue

            try:
                pos = int(p[1])
                ref, alt = p[2], p[3]

                if len(p) >= 8:
                    dp = int(p[4])
                    gl0, gl1, gl2 = float(p[5]), float(p[6]), float(p[7])
                else:
                    dp = -1
                    gl0, gl1, gl2 = float(p[4]), float(p[5]), float(p[6])
            except ValueError:
                continue

            rows.append((chrom, pos, ref, alt, dp, gl0, gl1, gl2))

    print(f"Loaded {len(rows)} SNPs")

    if not rows:
        with gzip.open(args.out_statepath_gz, "wt") as out:
            out.write("")
        with gzip.open(args.out_tracts_bed_gz, "wt") as out:
            out.write("")
        return

    # Sort by chromosome and position
    chr_order = {f"chr{i}": i for i in range(1, 11)}
    rows.sort(key=lambda x: (chr_order.get(x[0], 99), x[1]))

    # Build observation arrays
    chroms = []
    meta = []
    obs_ll = []
    thetas = [0.0]
    prev_chr = None
    prev_cm = None

    for i, (chrom, pos, ref, alt, dp, gl0, gl1, gl2) in enumerate(rows):
        lRR = gl0 * math.log(10.0)
        lRH = gl1 * math.log(10.0)
        lHH = gl2 * math.log(10.0)

        lRR, lRH, lHH = apply_emission_adjustments(
            lRR, lRH, lHH,
            args.eta_hh_from_rh,
            args.eta_rr_from_rh,
            args.rh_penalty
        )
        obs_ll.append([lRR, lRH, lHH])

        bp_map, cm_map = maps[chrom]
        cm_here = interp_cm(bp_map, cm_map, pos)

        if i == 0 or chrom != prev_chr:
            thetas.append(1e-6)
        else:
            d_cm = abs(cm_here - prev_cm)
            d_m = max(d_cm / 100.0, args.min_morgan)
            th = haldane_theta(d_m)
            thetas.append(max(th, 1e-8))

        prev_chr = chrom
        prev_cm = cm_here
        chroms.append(chrom)
        meta.append((chrom, pos, ref, alt, dp))

    # Run HMM
    print(f"Running {args.method} Viterbi...")

    if args.method == "majority":
        path = viterbi_majority_rigidity(obs_ll, chroms, args.rigidity)
    elif args.method == "rigidity":
        path = viterbi_rigidity(obs_ll, thetas, pi, args.rho, args.rigidity, chroms)
    else:
        path = viterbi_standard(obs_ll, thetas, pi, args.rho)

    # Compute posteriors (for output, using standard forward-backward)
    post = forward_backward(obs_ll, thetas, pi, args.rho)

    # Count states
    state_counts = [0, 0, 0]
    for s in path:
        state_counts[s] += 1
    print(f"State counts: RR={state_counts[0]}, RH={state_counts[1]}, HH={state_counts[2]}")

    # Count segments
    n_segments = 1
    for i in range(1, len(path)):
        if path[i] != path[i-1] or chroms[i] != chroms[i-1]:
            n_segments += 1
    print(f"Total segments: {n_segments}")

    # Write output
    with gzip.open(args.out_statepath_gz, "wt") as out:
        out.write("CHROM\tPOS\tREF\tALT\tDP\tSTATE\tP_RR\tP_RH\tP_HH\tTHETA\n")
        for i, (chrom, pos, ref, alt, dp) in enumerate(meta):
            st = STATES[path[i]]
            pRR, pRH, pHH = post[i]
            out.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{dp}\t{st}\t{pRR:.6g}\t{pRH:.6g}\t{pHH:.6g}\t{thetas[i+1]:.6g}\n")

    # Write tracts
    with gzip.open(args.out_tracts_bed_gz, "wt") as out:
        out.write("chrom\tstart\tend\tstate\n")
        cur_chr = meta[0][0]
        cur_state = path[0]
        cur_start = meta[0][1]
        cur_end = meta[0][1]

        for i in range(1, len(meta)):
            chrom, pos = meta[i][0], meta[i][1]
            st = path[i]

            if chrom != cur_chr or st != cur_state:
                out.write(f"{cur_chr}\t{max(cur_start-1,0)}\t{cur_end}\t{STATES[cur_state]}\n")
                cur_chr = chrom
                cur_state = st
                cur_start = pos
                cur_end = pos
            else:
                cur_end = pos

        out.write(f"{cur_chr}\t{max(cur_start-1,0)}\t{cur_end}\t{STATES[cur_state]}\n")

    print(f"\nOutput: {args.out_statepath_gz}")
    print(f"Tracts: {args.out_tracts_bed_gz}")
    print("Done!")


if __name__ == "__main__":
    main()
