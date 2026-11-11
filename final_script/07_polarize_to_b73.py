#!/usr/bin/env python3
"""
07_polarize_to_b73.py - put every ANGSD SNP on the B73 reference allele.

ANGSD's BCF stores alleles as MAJOR,MINOR (inferred from all pools), not as the
reference base. This script looks up the B73 base in the FASTA (via its .fai)
and writes one line per kept SNP:

    CHROM  POS  B73  ALT  SWAP
    SWAP = 1 when ANGSD's major allele was NOT the B73 base (alleles flipped here)

Dropped and counted in the log: B73 base not A/C/G/T, neither allele equals the
B73 base, non-SNP alleles, unknown contig / out-of-range position, and
unsorted or duplicated positions.

stdin : CHROM POS A1 A2  (from: bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n')
stdout: CHROM POS B73 ALT SWAP  (no header)
Standard library only.
"""
import argparse
import mmap
import sys


def read_fai(path):
    fai = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            fai[f[0]] = (int(f[1]), int(f[2]), int(f[3]), int(f[4]))  # length, offset, linebases, linewidth
    return fai


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True, help="B73 reference FASTA (needs .fai)")
    ap.add_argument("--log", help="write counts here (default: stderr)")
    args = ap.parse_args()

    fai = read_fai(args.fasta + ".fai")
    with open(args.fasta, "rb") as fh:
        mm = mmap.mmap(fh.fileno(), 0, access=mmap.ACCESS_READ)

    acgt = {"A", "C", "G", "T"}
    n = dict(read=0, kept=0, swapped=0, not_snp=0, b73_not_acgt=0,
             neither_allele_is_b73=0, bad_position=0, unsorted_or_dup=0)
    write = sys.stdout.write
    last_chrom, last_pos = None, 0

    for line in sys.stdin:
        n["read"] += 1
        f = line.rstrip("\n").split("\t")
        chrom, pos, a1, a2 = f[0], int(f[1]), f[2].upper(), f[3].upper()

        if chrom == last_chrom and pos <= last_pos:
            n["unsorted_or_dup"] += 1
            continue
        last_chrom, last_pos = chrom, pos

        if a1 not in acgt or a2 not in acgt or a1 == a2:
            n["not_snp"] += 1
            continue
        if chrom not in fai or not 1 <= pos <= fai[chrom][0]:
            n["bad_position"] += 1
            continue

        _, off, lb, lw = fai[chrom]
        i = pos - 1
        b73 = chr(mm[off + (i // lb) * lw + (i % lb)]).upper()
        if b73 not in acgt:
            n["b73_not_acgt"] += 1
            continue

        if a1 == b73:
            alt, swap = a2, 0
        elif a2 == b73:
            alt, swap = a1, 1
        else:
            n["neither_allele_is_b73"] += 1
            continue

        write(f"{chrom}\t{pos}\t{b73}\t{alt}\t{swap}\n")
        n["kept"] += 1
        n["swapped"] += swap

    report = "".join(f"{k}\t{v}\n" for k, v in n.items())
    if args.log:
        with open(args.log, "w") as fh:
            fh.write(report)
    else:
        sys.stderr.write(report)
    if n["kept"] == 0:
        sys.exit("ERROR: no SNPs kept")


if __name__ == "__main__":
    main()
