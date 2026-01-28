#!/usr/bin/env python3
import argparse, gzip, math
from bisect import bisect_right

STATES = ["RR","RH","HH"]

def read_map(map_path: str):
    bp, cm = [], []
    with open(map_path, "r") as f:
        for line in f:
            line=line.strip()
            if not line or line.startswith("#") or line.startswith("locus"):
                continue
            p=line.split()
            if len(p) < 2:
                continue
            try:
                b = int(float(p[0]))
                c = float(p[-1])  # last col as cM (works for pos 1.0 cM OR pos cM)
            except ValueError:
                continue
            bp.append(b); cm.append(c)
    if len(bp) < 2:
        raise RuntimeError(f"Map too small: {map_path}")
    z = sorted(zip(bp, cm))
    bp = [x for x,_ in z]
    cm = [y for _,y in z]
    return bp, cm

def interp_cm(bp_map, cm_map, pos: int) -> float:
    i = bisect_right(bp_map, pos) - 1
    if i < 0:
        x0,x1 = bp_map[0], bp_map[1]
        y0,y1 = cm_map[0], cm_map[1]
        return y0 if x1 == x0 else y0 + (pos-x0)*(y1-y0)/(x1-x0)
    if i >= len(bp_map)-1:
        x0,x1 = bp_map[-2], bp_map[-1]
        y0,y1 = cm_map[-2], cm_map[-1]
        return y1 if x1 == x0 else y0 + (pos-x0)*(y1-y0)/(x1-x0)
    x0,x1 = bp_map[i], bp_map[i+1]
    y0,y1 = cm_map[i], cm_map[i+1]
    return y0 if x1 == x0 else y0 + (pos-x0)*(y1-y0)/(x1-x0)

def haldane_theta(d_morgan: float) -> float:
    return 0.5*(1.0 - math.exp(-2.0*d_morgan))

def logsumexp3(a,b,c):
    m = max(a,b,c)
    return m + math.log(math.exp(a-m)+math.exp(b-m)+math.exp(c-m))

def mix_log(a, b, eta):
    # log((1-eta)exp(a) + eta exp(b))
    if eta <= 0: return a
    if eta >= 1: return b
    m = max(a,b)
    return m + math.log((1-eta)*math.exp(a-m) + eta*math.exp(b-m))

def apply_emission_adjustments(lRR, lRH, lHH, eta_hh_from_rh, eta_rr_from_rh, rh_penalty):
    # HH rescue
    lHH2 = mix_log(lHH, lRH, eta_hh_from_rh)
    # optional RR tolerance (usually small)
    lRR2 = mix_log(lRR, lRH, eta_rr_from_rh)
    # RH penalty
    lRH2 = lRH - rh_penalty
    return lRR2, lRH2, lHH2

def build_A(theta, pi, rho):
    # sticky+stationary transitions
    s = rho * theta
    if s > 0.25: s = 0.25
    if s < 1e-12: s = 1e-12
    A = [[0.0]*3 for _ in range(3)]
    for i in range(3):
        A[i][i] = 1.0 - s
        denom = 1.0 - pi[i]
        for j in range(3):
            if j == i: continue
            A[i][j] = s * (pi[j] / denom)
    return A

def forward_backward(obs_ll, thetas, pi, rho):
    n = len(obs_ll)
    def slog(x): return -1e300 if x <= 0 else math.log(x)

    alpha=[[-1e300]*3 for _ in range(n)]
    c=[0.0]*n

    for s in range(3):
        alpha[0][s] = slog(pi[s]) + obs_ll[0][s]
    c[0]=logsumexp3(alpha[0][0], alpha[0][1], alpha[0][2])
    for s in range(3):
        alpha[0][s]-=c[0]

    for t in range(1,n):
        A = build_A(thetas[t], pi, rho)
        for j in range(3):
            v0 = alpha[t-1][0] + slog(A[0][j])
            v1 = alpha[t-1][1] + slog(A[1][j])
            v2 = alpha[t-1][2] + slog(A[2][j])
            alpha[t][j] = logsumexp3(v0,v1,v2) + obs_ll[t][j]
        c[t]=logsumexp3(alpha[t][0], alpha[t][1], alpha[t][2])
        for j in range(3):
            alpha[t][j]-=c[t]

    beta=[[0.0]*3 for _ in range(n)]
    for t in range(n-2, -1, -1):
        A = build_A(thetas[t+1], pi, rho)
        for i in range(3):
            v0 = math.log(max(A[i][0],1e-300)) + obs_ll[t+1][0] + beta[t+1][0]
            v1 = math.log(max(A[i][1],1e-300)) + obs_ll[t+1][1] + beta[t+1][1]
            v2 = math.log(max(A[i][2],1e-300)) + obs_ll[t+1][2] + beta[t+1][2]
            beta[t][i] = logsumexp3(v0,v1,v2) - c[t+1]

    post=[[0.0]*3 for _ in range(n)]
    for t in range(n):
        v0 = alpha[t][0]+beta[t][0]
        v1 = alpha[t][1]+beta[t][1]
        v2 = alpha[t][2]+beta[t][2]
        z = logsumexp3(v0,v1,v2)
        post[t][0]=math.exp(v0-z)
        post[t][1]=math.exp(v1-z)
        post[t][2]=math.exp(v2-z)
    return post

def viterbi(obs_ll, thetas, pi, rho):
    n=len(obs_ll)
    def slog(x): return -1e300 if x <= 0 else math.log(x)

    dp=[[-1e300]*3 for _ in range(n)]
    ptr=[[0]*3 for _ in range(n)]

    for s in range(3):
        dp[0][s]=slog(pi[s])+obs_ll[0][s]

    for t in range(1,n):
        A = build_A(thetas[t], pi, rho)
        for j in range(3):
            best=-1e300; best_i=0
            for i in range(3):
                val = dp[t-1][i] + slog(A[i][j])
                if val > best:
                    best=val; best_i=i
            dp[t][j]=best + obs_ll[t][j]
            ptr[t][j]=best_i

    last=max(range(3), key=lambda s: dp[n-1][s])
    path=[last]
    for t in range(n-1,0,-1):
        last=ptr[t][last]
        path.append(last)
    path.reverse()
    return path

def run_filter(path, chroms, positions, min_run_hh, min_run_rh):
    # Flip short runs:
    # - HH run < min_run_hh -> RH
    # - RH run < min_run_rh -> RR
    out=path[:]
    n=len(path)
    i=0
    while i<n:
        j=i
        while j+1<n and chroms[j+1]==chroms[i] and out[j+1]==out[i]:
            j+=1
        run_len=j-i+1
        st=out[i]
        if st==2 and run_len < min_run_hh:  # HH
            for k in range(i,j+1): out[k]=1
        if st==1 and run_len < min_run_rh:  # RH
            for k in range(i,j+1): out[k]=0
        i=j+1
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_gl_tsv_gz", required=True)
    ap.add_argument("--map_dir", required=True)
    ap.add_argument("--out_statepath_gz", required=True)
    ap.add_argument("--out_tracts_bed_gz", required=True)

    ap.add_argument("--prior_rr", type=float, default=0.92)
    ap.add_argument("--prior_rh", type=float, default=0.02)
    ap.add_argument("--prior_hh", type=float, default=0.06)

    ap.add_argument("--min_morgan", type=float, default=1e-8)

    # NEW knobs
    ap.add_argument("--rho", type=float, default=10.0, help="rigidity: lower = fewer switches, higher = more switches")
    ap.add_argument("--eta_hh_from_rh", type=float, default=0.20, help="HH rescue: HH borrows from RH (0-0.5 typical)")
    ap.add_argument("--eta_rr_from_rh", type=float, default=0.00, help="optional RR borrow from RH (usually 0)")
    ap.add_argument("--rh_penalty", type=float, default=1.0, help="RH penalty in NATURAL log units (try 0.5-2.0)")

    ap.add_argument("--min_run_hh", type=int, default=3, help="post-filter: HH run must be >= this many SNPs")
    ap.add_argument("--min_run_rh", type=int, default=3, help="post-filter: RH run must be >= this many SNPs")
    args = ap.parse_args()

    s = args.prior_rr + args.prior_rh + args.prior_hh
    pi = [args.prior_rr/s, args.prior_rh/s, args.prior_hh/s]

    # Load maps
    maps={}
    for c in range(1,11):
        chrom=f"chr{c}"
        maps[chrom]=read_map(f"{args.map_dir}/{chrom}.map")

    rows=[]
    with gzip.open(args.in_gl_tsv_gz, "rt") as f:
        for line in f:
            line=line.strip()
            if not line: continue
            p=line.split("\t")
            if len(p) < 8: continue
            chrom=p[0]
            if chrom not in maps: continue
            try:
                pos=int(p[1]); ref=p[2]; alt=p[3]; dp=int(p[4])
                gl0=float(p[5]); gl1=float(p[6]); gl2=float(p[7])
            except ValueError:
                continue
            rows.append((chrom,pos,ref,alt,dp,gl0,gl1,gl2))

    if not rows:
        with gzip.open(args.out_statepath_gz, "wt") as out:
            out.write("")
        with gzip.open(args.out_tracts_bed_gz, "wt") as out:
            out.write("")
        return

    chr_order={f"chr{i}":i for i in range(1,11)}
    rows.sort(key=lambda x:(chr_order.get(x[0],99), x[1]))

    chroms=[]; poss=[]; meta=[]; obs_ll=[]; thetas=[0.0]
    prev_chr=None; prev_cm=None

    for i,(chrom,pos,ref,alt,dp,gl0,gl1,gl2) in enumerate(rows):
        # convert log10 GL to ln-likelihood (relative ok)
        lRR = gl0 * math.log(10.0)
        lRH = gl1 * math.log(10.0)
        lHH = gl2 * math.log(10.0)

        lRR,lRH,lHH = apply_emission_adjustments(
            lRR,lRH,lHH,
            eta_hh_from_rh=args.eta_hh_from_rh,
            eta_rr_from_rh=args.eta_rr_from_rh,
            rh_penalty=args.rh_penalty
        )
        obs_ll.append([lRR,lRH,lHH])

        bp_map, cm_map = maps[chrom]
        cm_here = interp_cm(bp_map, cm_map, pos)

        if i==0 or chrom != prev_chr:
            thetas.append(1e-6)
        else:
            d_cm = abs(cm_here - prev_cm)
            d_m = max(d_cm/100.0, args.min_morgan)
            th = haldane_theta(d_m)
            if th < 1e-8: th = 1e-8
            thetas.append(th)

        prev_chr=chrom; prev_cm=cm_here
        chroms.append(chrom); poss.append(pos)
        meta.append((chrom,pos,ref,alt,dp))

    post = forward_backward(obs_ll, thetas, pi, rho=args.rho)
    path = viterbi(obs_ll, thetas, pi, rho=args.rho)

    # post-run filter (“r sites”)
    path2 = run_filter(path, chroms, poss, args.min_run_hh, args.min_run_rh)

    # write statepath
    with gzip.open(args.out_statepath_gz, "wt") as out:
        out.write("CHROM\tPOS\tREF\tALT\tDP\tSTATE\tP_RR\tP_RH\tP_HH\tTHETA\n")
        for i,(chrom,pos,ref,alt,dp) in enumerate(meta):
            st = STATES[path2[i]]
            pRR,pRH,pHH = post[i]
            out.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{dp}\t{st}\t{pRR:.6g}\t{pRH:.6g}\t{pHH:.6g}\t{thetas[i+1]:.6g}\n")

    # write tracts
    with gzip.open(args.out_tracts_bed_gz, "wt") as out:
        out.write("chrom\tstart\tend\tstate\n")
        cur_chr=meta[0][0]
        cur_state=path2[0]
        cur_start=meta[0][1]
        cur_end=meta[0][1]
        for i in range(1,len(meta)):
            chrom=meta[i][0]; pos=meta[i][1]; st=path2[i]
            if chrom!=cur_chr or st!=cur_state:
                out.write(f"{cur_chr}\t{max(cur_start-1,0)}\t{cur_end}\t{STATES[cur_state]}\n")
                cur_chr=chrom; cur_state=st; cur_start=pos; cur_end=pos
            else:
                cur_end=pos
        out.write(f"{cur_chr}\t{max(cur_start-1,0)}\t{cur_end}\t{STATES[cur_state]}\n")

if __name__ == "__main__":
    main()
