#!/usr/bin/env python3
"""
25 — Inputs for association mapping (step 26) from the k-mer introgression calls.

Per plot (pool), per 100-kb window, from the step-21 pool states:
  dosage = 0 (B73), 1 (SEG), 2 (DONOR)  -> expected teosinte allele copies per plant in the plot
           (a SEG plot holds 1/4 AA, 1/2 AB, 1/4 BB plants: on average 1 copy)
  seg    = 1 if SEG, else 0             -> for the dominance test (step 26)
Plots kept: donor pools with FLAG ok (low-depth pools dropped) that have a phenotype.
Windows kept: at least --min_carriers plots carry donor (SEG or DONOR) there.
LOCO kinship: for each chromosome c, K_c = Z Z' / p from the dosage of all windows on the OTHER
chromosomes (Z = column-centred dosage, windows with any variation), then scaled to mean diagonal 1.

IDs: the phenotype's plot name (column 1, e.g. Zv.0200_P3_P2_P1.1.1.1), linked to the pool via
its Sequencing_ID column, so step 26 can join on column 1 directly.

Outputs in <out_dir>:
  geno_dosage.csv      Line, then one column per window "chrN:midpoint" (0/1/2)
  geno_seg.csv         same layout, SEG indicator (0/1)
  window_map.tsv       Marker chr start end pos n_carriers n_seg n_donor
  pheno.csv            the phenotype rows kept, plus TAXON, F1FAMILY (accession + F1 plant) and SAMPLE
  K_matrix_chr{c}.txt  LOCO kinship, first column ID, header row of IDs
"""
import argparse
import collections
import csv
import gzip
import os

import numpy as np

CODE = {"B73": 0, "SEG": 1, "DONOR": 2}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tracts", required=True, help="step 21 tracts_all.tsv.gz")
    ap.add_argument("--pool_summary", required=True, help="step 21 pool_summary.tsv")
    ap.add_argument("--pheno", required=True, help="CSV: plot name, traits..., Sequencing_ID")
    ap.add_argument("--sid_col", default="Sequencing_ID")
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--window", type=int, default=100000)
    ap.add_argument("--min_carriers", type=int, default=20)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)

    P = {r["SAMPLE"]: r for r in csv.DictReader(open(a.pool_summary), delimiter="\t")}
    with open(a.pheno, encoding="utf-8-sig") as f:
        rd = csv.reader(f); head = next(rd); rows = [r for r in rd if r]
    isid = head.index(a.sid_col)
    keep, drop = [], collections.Counter()
    for r in rows:
        s = r[isid]
        if s not in P: drop["no k-mer calls"] += 1; continue
        if P[s]["ROLE"] != "donor": drop["not a donor pool"] += 1; continue
        if P[s]["FLAG"] != "ok": drop["low-depth (NOISY)"] += 1; continue
        keep.append(r)
    print(f"phenotyped plots: {len(rows)} | kept: {len(keep)} | dropped: {dict(drop)}")
    sids = [r[isid] for r in keep]; lines = [r[0] for r in keep]

    # window grid from the tracts (each pool's tracts cover every chromosome end to end)
    clen = collections.defaultdict(int); tr = collections.defaultdict(list); want = set(sids)
    with gzip.open(a.tracts, "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            e = int(r["END"]); clen[r["CHROM"]] = max(clen[r["CHROM"]], e)
            if r["SAMPLE"] in want: tr[r["SAMPLE"]].append((r["CHROM"], int(r["START"]), e, CODE[r["STATE"]]))
    chroms = sorted(clen, key=lambda c: int(c.replace("chr", "")))
    W = a.window
    win = [(c, w * W, min((w + 1) * W, clen[c])) for c in chroms for w in range(-(-clen[c] // W))]
    idx = {}; k = 0
    for c in chroms:
        idx[c] = k; k += -(-clen[c] // W)
    D = np.zeros((len(sids), len(win)), dtype=np.int8)
    for i, s in enumerate(sids):
        for c, b, e, st in tr[s]:
            if st == 0: continue
            w0 = idx[c] + (b + W // 2) // W                     # windows whose midpoint lies in [b, e)
            w1 = idx[c] + max((e - W // 2 + W - 1) // W, 0)
            D[i, w0:w1] = st
    seg = (D == 1); car = (D > 0)
    ncar = car.sum(0)
    ok = ncar >= a.min_carriers
    print(f"windows: {len(win)} | with >= {a.min_carriers} carrier plots: {int(ok.sum())}")

    names = [f"{c}:{(b + e) // 2}" for c, b, e in win]
    with open(f"{a.out_dir}/window_map.tsv", "w") as f:
        f.write("Marker\tchr\tstart\tend\tpos\tn_carriers\tn_seg\tn_donor\n")
        for j in np.where(ok)[0]:
            c, b, e = win[j]
            f.write(f"{names[j]}\t{c}\t{b + 1}\t{e}\t{(b + e) // 2}\t{ncar[j]}\t{seg[:, j].sum()}\t{(D[:, j] == 2).sum()}\n")
    for fn, M in (("geno_dosage.csv", D), ("geno_seg.csv", seg.astype(np.int8))):
        with open(f"{a.out_dir}/{fn}", "w") as f:
            f.write("Line," + ",".join(names[j] for j in np.where(ok)[0]) + "\n")
            Mk = M[:, ok]
            for i, l in enumerate(lines): f.write(l + "," + ",".join(map(str, Mk[i])) + "\n")
    with open(f"{a.out_dir}/pheno.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(head + ["TAXON", "F1FAMILY", "SAMPLE"])
        for r in keep: w.writerow(r + [P[r[isid]]["GROUP"], "_".join(r[0].split("_")[:2]), r[isid]])

    # LOCO kinship from the dosage of all variable windows on the other chromosomes
    var = D.std(0) > 0
    chrom_of = np.array([c for c, _, _ in win])
    for c in chroms:
        cols = var & (chrom_of != c)
        Z = D[:, cols].astype(float); Z -= Z.mean(0)
        K = Z @ Z.T / cols.sum(); K /= np.mean(np.diag(K))
        cn = c.replace("chr", "")
        with open(f"{a.out_dir}/K_matrix_chr{cn}.txt", "w") as f:
            f.write("ID," + ",".join(lines) + "\n")
            for i, l in enumerate(lines): f.write(l + "," + ",".join(f"{x:.6f}" for x in K[i]) + "\n")
    print(f"WROTE: {a.out_dir}/geno_dosage.csv geno_seg.csv window_map.tsv pheno.csv K_matrix_chr1..{len(chroms)}.txt")


if __name__ == "__main__":
    main()
