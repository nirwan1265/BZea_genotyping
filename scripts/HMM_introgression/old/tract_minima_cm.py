#!/usr/bin/env python3
import os, sys, gzip, math
from bisect import bisect_right
import glob

# ---------- map utils ----------
def read_map_chr(map_path):
    xs, ys = [], []
    with open(map_path) as f:
        for line in f:
            line=line.strip()
            if not line or line.startswith("#") or line.startswith("locus"):
                continue
            p=line.split()
            # chrN.map typical: bp 1.0 cM
            if len(p) >= 3:
                bp = int(float(p[0]))
                cm = float(p[2])
            elif len(p) == 2:
                bp = int(float(p[0]))
                cm = float(p[1])
            else:
                continue
            xs.append(bp); ys.append(cm)
    z = sorted(zip(xs, ys))
    xs = [a for a,_ in z]
    ys = [b for _,b in z]
    return xs, ys

def interp_cm(pos, xs, ys):
    if pos <= xs[0]: return ys[0]
    if pos >= xs[-1]: return ys[-1]
    i = bisect_right(xs, pos) - 1
    x0,x1 = xs[i], xs[i+1]
    y0,y1 = ys[i], ys[i+1]
    if x1 == x0: return y0
    t = (pos-x0)/(x1-x0)
    return y0 + t*(y1-y0)

def quantile(sorted_vals, q):
    if not sorted_vals:
        return float("nan")
    n = len(sorted_vals)
    if n == 1:
        return sorted_vals[0]
    # linear interpolation
    x = (n-1)*q
    i = int(math.floor(x))
    j = int(math.ceil(x))
    if i == j:
        return sorted_vals[i]
    return sorted_vals[i]*(j-x) + sorted_vals[j]*(x-i)

# ---------- tract extraction ----------
def runs_from_statepath(path_gz):
    # returns list of dict runs with bp_start/end and state
    runs = []
    with gzip.open(path_gz, "rt") as f:
        header = f.readline()
        if not header:
            return runs
        for line in f:
            line=line.strip()
            if not line: continue
            chrom,pos,ref,alt,dp,state = line.split("\t", 6)[:6]
            pos = int(pos)
            runs.append((chrom,pos,state))
    return runs

def make_runs(site_list):
    # site_list: [(chrom,pos,state)] sorted already
    out=[]
    if not site_list:
        return out
    cur_chr, cur_state = site_list[0][0], site_list[0][2]
    cur_start = site_list[0][1]
    cur_end   = site_list[0][1]
    n_sites   = 1
    for chrom,pos,state in site_list[1:]:
        if chrom == cur_chr and state == cur_state:
            cur_end = pos
            n_sites += 1
        else:
            out.append((cur_chr, cur_state, cur_start, cur_end, n_sites))
            cur_chr, cur_state = chrom, state
            cur_start = cur_end = pos
            n_sites = 1
    out.append((cur_chr, cur_state, cur_start, cur_end, n_sites))
    return out

def run_len_cm(run, maps):
    chrom,state,bp0,bp1,n = run
    xs, ys = maps[chrom]
    cm0 = interp_cm(bp0, xs, ys)
    cm1 = interp_cm(bp1, xs, ys)
    d = abs(cm1 - cm0)
    return d

def run_len_bp(run):
    chrom,state,bp0,bp1,n = run
    return max(0, bp1 - bp0)

# donor sets
DONOR_HH = set(["HH"])
DONOR_HH_RH = set(["HH","RH"])

def merge_runs_by_state(runs, donor_set):
    """
    Collapse runs into donor vs non-donor tracts:
    returns runs labeled DONOR or NONDONOR
    """
    out=[]
    for chrom,state,bp0,bp1,n in runs:
        lab = "DONOR" if state in donor_set else "NON"
        out.append((chrom, lab, bp0, bp1, n))
    # merge adjacent with same chrom+lab
    merged=[]
    if not out:
        return merged
    c_chr,c_lab,c0,c1,cn = out[0]
    for chrom,lab,bp0,bp1,n in out[1:]:
        if chrom==c_chr and lab==c_lab:
            c1 = bp1
            cn += n
        else:
            merged.append((c_chr,c_lab,c0,c1,cn))
            c_chr,c_lab,c0,c1,cn = chrom,lab,bp0,bp1,n
    merged.append((c_chr,c_lab,c0,c1,cn))
    return merged

# ---------- main ----------
def main():
    if len(sys.argv) != 4:
        print("Usage: tract_minima_cm.py <statepaths_dir> <NAM_map_dir> <out.tsv>", file=sys.stderr)
        sys.exit(2)

    statedir, mapdir, out_tsv = sys.argv[1:4]

    # load maps chr1..chr10
    maps={}
    for i in range(1,11):
        chrom=f"chr{i}"
        mp=os.path.join(mapdir, f"{chrom}.map")
        if not os.path.exists(mp):
            raise SystemExit(f"Missing map: {mp}")
        maps[chrom]=read_map_chr(mp)

    files = sorted(glob.glob(os.path.join(statedir, "*.statepath.tsv.gz")))
    if not files:
        raise SystemExit(f"No *.statepath.tsv.gz found in {statedir}")

    with open(out_tsv, "w") as out:
        out.write("\t".join([
            "sample",
            "n_sites",
            "min_RH_cM","q05_RH_cM","median_RH_cM",
            "min_HH_cM","q05_HH_cM","median_HH_cM",
            "min_donorHH_cM","q05_donorHH_cM","median_donorHH_cM",
            "min_donorHH_RH_cM","q05_donorHH_RH_cM","median_donorHH_RH_cM",
            "min_RH_bp","min_HH_bp","min_donorHH_bp","min_donorHH_RH_bp"
        ]) + "\n")

        for fp in files:
            base=os.path.basename(fp)
            sample=base.replace(".statepath.tsv.gz","")

            sites = runs_from_statepath(fp)
            if not sites:
                out.write(sample + "\t0\t" + "\t".join(["NA"]*15) + "\n")
                continue

            # sort sites
            chr_order={f"chr{i}":i for i in range(1,11)}
            sites.sort(key=lambda x: (chr_order.get(x[0],99), x[1]))

            runs = make_runs(sites)

            # RH runs and HH runs
            RH = [r for r in runs if r[1]=="RH"]
            HH = [r for r in runs if r[1]=="HH"]

            RH_cm = sorted([run_len_cm(r,maps) for r in RH if run_len_cm(r,maps) > 0])
            HH_cm = sorted([run_len_cm(r,maps) for r in HH if run_len_cm(r,maps) > 0])

            RH_bp = sorted([run_len_bp(r) for r in RH if run_len_bp(r) > 0])
            HH_bp = sorted([run_len_bp(r) for r in HH if run_len_bp(r) > 0])

            # donor tracts (collapse HH vs others; and HH+RH vs others)
            donorHH_runs = merge_runs_by_state(runs, DONOR_HH)
            donorHH_runs = [r for r in donorHH_runs if r[1]=="DONOR"]
            donorHH_cm = sorted([run_len_cm((r[0],"HH",r[2],r[3],r[4]),maps) for r in donorHH_runs
                                 if run_len_cm((r[0],"HH",r[2],r[3],r[4]),maps) > 0])
            donorHH_bp = sorted([max(0,r[3]-r[2]) for r in donorHH_runs if (r[3]-r[2])>0])

            donorHHRH_runs = merge_runs_by_state(runs, DONOR_HH_RH)
            donorHHRH_runs = [r for r in donorHHRH_runs if r[1]=="DONOR"]
            donorHHRH_cm = sorted([run_len_cm((r[0],"HH",r[2],r[3],r[4]),maps) for r in donorHHRH_runs
                                   if run_len_cm((r[0],"HH",r[2],r[3],r[4]),maps) > 0])
            donorHHRH_bp = sorted([max(0,r[3]-r[2]) for r in donorHHRH_runs if (r[3]-r[2])>0])

            def safe_min(v): return v[0] if v else float("nan")
            def safe_med(v): return quantile(v,0.5) if v else float("nan")
            def safe_q05(v): return quantile(v,0.05) if v else float("nan")

            out.write("\t".join(map(str,[
                sample,
                len(sites),
                safe_min(RH_cm), safe_q05(RH_cm), safe_med(RH_cm),
                safe_min(HH_cm), safe_q05(HH_cm), safe_med(HH_cm),
                safe_min(donorHH_cm), safe_q05(donorHH_cm), safe_med(donorHH_cm),
                safe_min(donorHHRH_cm), safe_q05(donorHHRH_cm), safe_med(donorHHRH_cm),
                safe_min(RH_bp), safe_min(HH_bp), safe_min(donorHH_bp), safe_min(donorHHRH_bp)
            ])) + "\n")

if __name__=="__main__":
    main()
