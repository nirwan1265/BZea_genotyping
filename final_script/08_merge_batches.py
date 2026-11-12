#!/usr/bin/env python3
"""
08_merge_batches.py - merge ANGSD -dumpCounts 4 outputs from several BAM batches
into one .pos.gz / .counts.gz pair, as if ANGSD had been run on all BAMs at once.

Used by 08_angsd_counts_by_part.sh when the open-file limit forces the BAM list
to be split. Each batch was run on a consecutive slice of the BAM list with the
same -r and -sites, so its columns are that slice's pools in order. ANGSD only
writes sites where the batch had reads, so a site missing from a batch gets
zeros for that batch's pools. totDepth is summed over batches.

usage: 08_merge_batches.py --out PREFIX --batch PREFIX1 N1 --batch PREFIX2 N2 ...
       (N = number of BAMs in that batch, in BAM-list order)
Standard library only.
"""
import argparse
import gzip
import sys


def rows(prefix, n):
    """yield (chrom, pos, depth, counts_string_with_trailing_tab) for one batch"""
    with gzip.open(prefix + ".pos.gz", "rt") as fp, gzip.open(prefix + ".counts.gz", "rt") as fc:
        fp.readline(); fc.readline()                      # headers
        for lp, lc in zip(fp, fc):
            chrom, pos, depth = lp.rstrip("\n").split("\t")[:3]
            c = lc.rstrip("\n")
            if not c.endswith("\t"):
                c += "\t"
            if c.count("\t") != 4 * n:
                sys.exit(f"ERROR: {prefix}.counts.gz has {c.count(chr(9))} columns at {chrom}:{pos}, expected {4 * n}")
            yield chrom, int(pos), int(depth), c
        if fp.readline() or fc.readline():
            sys.exit(f"ERROR: {prefix}.pos.gz and .counts.gz have different numbers of rows")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--batch", nargs=2, action="append", required=True, metavar=("PREFIX", "N"))
    a = ap.parse_args()

    batches = [(p, int(n)) for p, n in a.batch]
    zeros = ["0\t" * (4 * n) for _, n in batches]
    its = [rows(p, n) for p, n in batches]
    cur = [next(it, None) for it in its]
    ntot = sum(n for _, n in batches)

    with gzip.open(a.out + ".pos.gz", "wt", compresslevel=3) as op, \
         gzip.open(a.out + ".counts.gz", "wt", compresslevel=3) as oc:
        op.write("chr\tpos\ttotDepth\n")
        oc.write("".join(f"ind{i}_A\tind{i}_C\tind{i}_G\tind{i}_T\t" for i in range(ntot)) + "\n")
        nrow = 0
        last = 0
        while any(r is not None for r in cur):
            pos = min(r[1] for r in cur if r is not None)
            if pos <= last:
                sys.exit(f"ERROR: positions not increasing at {pos}")
            last = pos
            chrom, depth, parts = None, 0, []
            for k, r in enumerate(cur):
                if r is not None and r[1] == pos:
                    chrom = r[0]
                    depth += r[2]
                    parts.append(r[3])
                    cur[k] = next(its[k], None)
                else:
                    parts.append(zeros[k])
            op.write(f"{chrom}\t{pos}\t{depth}\n")
            oc.write("".join(parts) + "\n")
            nrow += 1
    sys.stderr.write(f"merged {len(batches)} batches, {ntot} pools, {nrow} sites\n")


if __name__ == "__main__":
    main()
