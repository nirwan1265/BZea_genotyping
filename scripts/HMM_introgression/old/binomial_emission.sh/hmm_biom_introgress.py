#!/usr/bin/env python3
import gzip, math, sys, os

STATES = ["RR","RH","HH"]
S2I = {s:i for i,s in enumerate(STATES)}

def logsumexp(vals):
    m = max(vals)
    if m == -1e300:
        return m
    return m + math.log(sum(math.exp(v-m) for v in vals))

def safe_log(x):
    if x <= 0.0:
        return -1e300
    return math.log(x)

def haldane_theta(delta_cM):
    # delta_cM >= 0
    d = delta_cM / 100.0
    # theta = (1 - exp(-2d))/2
    return 0.5 * (1.0 - math.exp(-2.0 * d))

def read_map_chr(map_path):
    # Accept either 3-col (pos, 1.0, cM) or 2-col (pos, cM) or NAM_genetic_map.txt style
    xs = []
    ys = []
    with open(map_path, "r") as f:
        for line in f:
            line=line.strip()
            if not line or line.startswith("#") or line.startswith("locus"):
                continue
            parts = line.split()
            # chrN.map example: pos 1.0 cM
            if len(parts) >= 3:
                pos = int(float(parts[0]))
                cm  = float(parts[2])
            elif len(parts) == 2:
                pos = int(float(parts[0]))
                cm  = float(parts[1])
            else:
                continue
            xs.append(pos)
            ys.append(cm)
    # Must be sorted
    zipped = sorted(zip(xs, ys))
    xs = [p for p,_ in zipped]
    ys = [c for _,c in zipped]
    return xs, ys

def interp_cm(pos, xs, ys):
    # linear interpolation, clamp to edges
    if pos <= xs[0]:
        return ys[0]
    if pos >= xs[-1]:
        return ys[-1]
    # binary search
    lo, hi = 0, len(xs)-1
    while hi-lo > 1:
        mid = (lo+hi)//2
        if xs[mid] <= pos:
            lo = mid
        else:
            hi = mid
    x0,x1 = xs[lo], xs[hi]
    y0,y1 = ys[lo], ys[hi]
    if x1 == x0:
        return y0
    t = (pos - x0) / (x1 - x0)
    return y0 + t*(y1-y0)

def log_binom_pmf(k, n, p):
    # log [ C(n,k) p^k (1-p)^(n-k) ]
    # using lgamma for stability
    if p <= 0.0:
        return 0.0 if k==0 else -1e300
    if p >= 1.0:
        return 0.0 if k==n else -1e300
    return (math.lgamma(n+1) - math.lgamma(k+1) - math.lgamma(n-k+1)
            + k*safe_log(p) + (n-k)*safe_log(1.0-p))

def build_transitions(theta, pi, rho=10.0):
    # "sticky + stationary" transition using switch prob s ~ rho*theta, capped
    s = rho * theta
    if s > 0.25:
        s = 0.25
    if s < 1e-12:
        s = 1e-12
    A = [[0.0]*3 for _ in range(3)]
    for i in range(3):
        A[i][i] = 1.0 - s
        rem = s
        denom = 1.0 - pi[i]
        for j in range(3):
            if j==i: 
                continue
            A[i][j] = rem * (pi[j]/denom)
    return A

def forward_backward(obs, thetas, pi, eps_rr=0.01, eps_hh=0.05, rho=10.0):
    # obs: list of tuples (refc, altc, dp)
    # emissions: alt count modeled as binomial(dp, p_state)
    # p_RR=eps_rr, p_RH=0.5, p_HH=1-eps_hh  (eps_hh > eps_rr allows ref-bias in HH)
    p_state = [eps_rr, 0.5, 1.0-eps_hh]

    n = len(obs)
    # log emissions
    E = [[0.0]*3 for _ in range(n)]
    for t,(refc,altc,dp) in enumerate(obs):
        k = altc
        for s in range(3):
            E[t][s] = log_binom_pmf(k, dp, p_state[s])

    # forward (log-space) with scaling by logsumexp
    alpha = [[-1e300]*3 for _ in range(n)]
    c = [0.0]*n

    for s in range(3):
        alpha[0][s] = safe_log(pi[s]) + E[0][s]
    c[0] = logsumexp(alpha[0])
    for s in range(3):
        alpha[0][s] -= c[0]

    for t in range(1,n):
        A = build_transitions(thetas[t], pi, rho=rho)  # theta for step t (prev->t)
        for j in range(3):
            vals = [alpha[t-1][i] + safe_log(A[i][j]) for i in range(3)]
            alpha[t][j] = logsumexp(vals) + E[t][j]
        c[t] = logsumexp(alpha[t])
        for j in range(3):
            alpha[t][j] -= c[t]

    # backward
    beta = [[0.0]*3 for _ in range(n)]
    for s in range(3):
        beta[n-1][s] = 0.0

    for t in range(n-2, -1, -1):
        A = build_transitions(thetas[t+1], pi, rho=rho)  # step t+1
        for i in range(3):
            vals = [safe_log(A[i][j]) + E[t+1][j] + beta[t+1][j] for j in range(3)]
            beta[t][i] = logsumexp(vals) - c[t+1]

    # posterior
    post = [[0.0]*3 for _ in range(n)]
    for t in range(n):
        vals = [alpha[t][s] + beta[t][s] for s in range(3)]
        z = logsumexp(vals)
        for s in range(3):
            post[t][s] = math.exp(vals[s] - z)

    return post, E, thetas

def viterbi(obs, thetas, pi, eps_rr=0.01, eps_hh=0.05, rho=10.0):
    p_state = [eps_rr, 0.5, 1.0-eps_hh]
    n = len(obs)
    E = [[0.0]*3 for _ in range(n)]
    for t,(refc,altc,dp) in enumerate(obs):
        k=altc
        for s in range(3):
            E[t][s]=log_binom_pmf(k, dp, p_state[s])

    dpv = [[-1e300]*3 for _ in range(n)]
    ptr = [[0]*3 for _ in range(n)]

    for s in range(3):
        dpv[0][s] = safe_log(pi[s]) + E[0][s]

    for t in range(1,n):
        A = build_transitions(thetas[t], pi, rho=rho)
        for j in range(3):
            best_i = 0
            best = -1e300
            for i in range(3):
                val = dpv[t-1][i] + safe_log(A[i][j])
                if val > best:
                    best = val
                    best_i = i
            dpv[t][j] = best + E[t][j]
            ptr[t][j] = best_i

    # backtrack
    last = max(range(3), key=lambda s: dpv[n-1][s])
    path = [last]
    for t in range(n-1,0,-1):
        last = ptr[t][last]
        path.append(last)
    path.reverse()
    return path

def read_counts_tsv_gz(path):
    rows = []
    with gzip.open(path, "rt") as f:
        for line in f:
            line=line.strip()
            if not line:
                continue
            parts=line.split("\t")
            # CHROM POS REF ALT REFc ALTc DP
            chrom=parts[0]
            pos=int(parts[1])
            ref=parts[2]
            alt=parts[3]
            refc=int(parts[4])
            altc=int(parts[5])
            dp=int(parts[6])
            rows.append((chrom,pos,ref,alt,refc,altc,dp))
    return rows

def main():
    if len(sys.argv) != 5:
        print("Usage: hmm_binom_introgress.py <counts.tsv.gz> <NAM_map_dir> <out.statepath.tsv.gz> <bc_generation>", file=sys.stderr)
        print("  bc_generation: e.g. BC2S3 -> pass 'BC2S3' (used only for priors)", file=sys.stderr)
        sys.exit(2)

    counts_gz, mapdir, out_gz, gen = sys.argv[1:5]

    # Priors: BC2S3 genome-wide (tuneable)
    # donor fraction at BC2 is ~12.5%, after selfing segments go HH/RH.
    # Start with heavy RR and small donor:
    pi = [0.97, 0.02, 0.01]  # RR, RH, HH

    # Emission bias:
    # eps_rr: allow a little alt in RR
    # eps_hh: allow MORE ref in HH (reference bias), crucial for your problem
    eps_rr = 0.01
    eps_hh = 0.10   # <-- key knob; increase if HH still disappears (0.10–0.20)
    rho    = 10.0   # transition scaling

    rows = read_counts_tsv_gz(counts_gz)

    # group by chromosome
    by_chr = {}
    for r in rows:
        by_chr.setdefault(r[0], []).append(r)

    out_tmp = out_gz + ".tmp"
    with gzip.open(out_tmp, "wt") as out:
        out.write("CHROM\tPOS\tREF\tALT\tDP\tSTATE\tP_RR\tP_RH\tP_HH\tTHETA\n")

        for chrom in sorted(by_chr.keys(), key=lambda x: (len(x), x)):
            arr = sorted(by_chr[chrom], key=lambda x: x[1])

            # load map
            map_path = os.path.join(mapdir, f"{chrom}.map")
            if not os.path.exists(map_path):
                # fallback: try chrN.map with chrom already "chrN"
                raise SystemExit(f"Missing genetic map for {chrom}: {map_path}")

            xs, ys = read_map_chr(map_path)

            # Build obs + thetas
            obs = []
            meta = []
            cms = []
            for (c,pos,ref,alt,refc,altc,dp) in arr:
                cm = interp_cm(pos, xs, ys)
                cms.append(cm)
                obs.append((refc,altc,dp))
                meta.append((c,pos,ref,alt,dp))

            thetas = [1e-6]*len(obs)
            for t in range(1,len(obs)):
                dcm = max(0.0, cms[t] - cms[t-1])
                theta = haldane_theta(dcm)
                if theta < 1e-8:
                    theta = 1e-8
                thetas[t] = theta

            post, _, _ = forward_backward(obs, thetas, pi, eps_rr=eps_rr, eps_hh=eps_hh, rho=rho)
            path = viterbi(obs, thetas, pi, eps_rr=eps_rr, eps_hh=eps_hh, rho=rho)

            for t in range(len(obs)):
                c,pos,ref,alt,dp = meta[t]
                st = STATES[path[t]]
                pRR, pRH, pHH = post[t]
                out.write(f"{c}\t{pos}\t{ref}\t{alt}\t{dp}\t{st}\t{pRR:.6g}\t{pRH:.6g}\t{pHH:.6g}\t{thetas[t]:.6g}\n")

    os.replace(out_tmp, out_gz)

if __name__ == "__main__":
    main()
