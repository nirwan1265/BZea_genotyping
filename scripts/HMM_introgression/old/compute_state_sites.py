#!/usr/bin/env python3
import argparse, gzip, os, sys, math

def parse_chr_len(s):
    # "chr1=308452471,chr2=243675191,..."
    d = {}
    for part in s.split(","):
        k,v = part.split("=")
        d[k.strip()] = int(v)
    return d

def midpoint(a, b):
    return (a + b) // 2  # floor

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--statepath_gz", required=True)
    ap.add_argument("--chr_len", required=True,
                    help="comma list: chr1=...,chr2=...,chr10=...")
    ap.add_argument("--out_tsv", required=True)
    args = ap.parse_args()

    chr_len = parse_chr_len(args.chr_len)

    sample = os.path.basename(args.statepath_gz).replace(".statepath.tsv.gz","")

    # site counts
    site = {"RR":0, "RH":0, "HH":0}
    # bp (genome-length) counts using midpoint interpolation
    bp = {"RR":0, "RH":0, "HH":0}

    prev_chr = None
    prev_pos = None
    prev_state = None
    seg_start = None

    with gzip.open(args.statepath_gz, "rt") as f:
        header = f.readline()
        if not header:
            raise SystemExit("Empty file")

        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            p = line.split("\t")
            chrom = p[0]
            pos = int(p[1])
            state = p[5]

            if state not in site:
                continue

            site[state] += 1

            # initialize first record
            if prev_chr is None:
                prev_chr = chrom
                prev_pos = pos
                prev_state = state
                seg_start = 1
                continue

            # chromosome change => close previous chromosome segment to chr_len
            if chrom != prev_chr:
                end_chr = chr_len.get(prev_chr)
                if end_chr is None:
                    raise SystemExit(f"Missing chr_len for {prev_chr}")
                if seg_start <= end_chr:
                    bp[prev_state] += (end_chr - seg_start + 1)

                # reset for new chromosome
                prev_chr = chrom
                prev_pos = pos
                prev_state = state
                seg_start = 1
                continue

            # same chromosome: close previous segment at midpoint(prev_pos, pos)
            end = midpoint(prev_pos, pos)
            if end < seg_start:
                end = seg_start
            bp[prev_state] += (end - seg_start + 1)

            # start next segment
            seg_start = end + 1
            prev_pos = pos
            prev_state = state

    # close last chromosome
    if prev_chr is not None:
        end_chr = chr_len.get(prev_chr)
        if end_chr is None:
            raise SystemExit(f"Missing chr_len for {prev_chr}")
        if seg_start <= end_chr:
            bp[prev_state] += (end_chr - seg_start + 1)

    total_bp = sum(chr_len[c] for c in chr_len.keys())
    total_sites = sum(site.values())

    # introgression definitions
    donor_bp_HH = bp["HH"]
    donor_bp_HH_RH = bp["HH"] + bp["RH"]

    donor_sites_HH = site["HH"]
    donor_sites_HH_RH = site["HH"] + site["RH"]

    with open(args.out_tsv, "w") as out:
        out.write(
            "sample\tRR_sites\tRH_sites\tHH_sites\ttotal_sites\t"
            "RR_bp\tRH_bp\tHH_bp\ttotal_bp\t"
            "pct_RR_sites\tpct_RH_sites\tpct_HH_sites\t"
            "pct_RR_bp\tpct_RH_bp\tpct_HH_bp\t"
            "pct_donor_bp_HH\tpct_donor_bp_HH_RH\t"
            "pct_donor_sites_HH\tpct_donor_sites_HH_RH\n"
        )
        def pct(x, tot):
            return 0.0 if tot == 0 else 100.0 * x / tot

        out.write(
            f"{sample}\t{site['RR']}\t{site['RH']}\t{site['HH']}\t{total_sites}\t"
            f"{bp['RR']}\t{bp['RH']}\t{bp['HH']}\t{total_bp}\t"
            f"{pct(site['RR'], total_sites):.6f}\t{pct(site['RH'], total_sites):.6f}\t{pct(site['HH'], total_sites):.6f}\t"
            f"{pct(bp['RR'], total_bp):.6f}\t{pct(bp['RH'], total_bp):.6f}\t{pct(bp['HH'], total_bp):.6f}\t"
            f"{pct(donor_bp_HH, total_bp):.6f}\t{pct(donor_bp_HH_RH, total_bp):.6f}\t"
            f"{pct(donor_sites_HH, total_sites):.6f}\t{pct(donor_sites_HH_RH, total_sites):.6f}\n"
        )

if __name__ == "__main__":
    main()
