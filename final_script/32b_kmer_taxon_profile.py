#!/usr/bin/env python3
"""
32b — Teosinte k-mer profile of every pool, from sourmash sketches (32a). Reference-free check of which taxon
each pool's donor DNA comes from, and input for the PCA (32c).

Sets (sourmash FracMinHash hashes, k = 31):
  non-B73 hashes of an assembly  A* = A - B73
  taxon set                      T_t = union of A* over the taxon's assemblies
  taxon-specific set             U_t = T_t - union of the other taxa's T
Per pool (hash set H; the pool is ~90% B73, so only non-B73 hashes carry donor information):
  containment of each assembly   c_A = |H & A*| / |A*|        (input for PCA; one column per assembly)
  taxon-specific hits            n_t = |H & U_t|, f_t = n_t / |U_t|
  teosinte content               f_teo = |H & union T| / |union T|   (scales with donor share and depth)
  assigned taxon                 argmax f_t; ratio = best f_t / second best f_t
Usage: python3 32b_kmer_taxon_profile.py --refs <D32/refs> --pools <D32/pools> --assemblies NAME:TAXON ...
         --pools_tsv <D11/pools.tsv> --out <D32/kmer_taxon_profile.tsv>
"""
import argparse, csv, glob, os
import numpy as np
import sourmash


def hashes(path):
    sig = next(iter(sourmash.load_file_as_signatures(path, ksize=31)))
    return np.fromiter(sig.minhash.hashes, dtype=np.uint64)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--refs", required=True); ap.add_argument("--pools", required=True)
    ap.add_argument("--assemblies", nargs="+", required=True, help="NAME:TAXON (as in DONOR_ASSEMBLIES)")
    ap.add_argument("--pools_tsv", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    b73 = hashes(os.path.join(a.refs, "B73.sig.gz"))
    asm, tax_of = {}, {}
    for x in a.assemblies:
        name, tax = x.split(":")[:2]
        asm[name] = np.setdiff1d(hashes(os.path.join(a.refs, f"{name}.sig.gz")), b73, assume_unique=False)
        tax_of[name] = tax
    taxa = sorted(set(tax_of.values()))
    T = {t: np.unique(np.concatenate([asm[n] for n in asm if tax_of[n] == t])) for t in taxa}
    U = {t: np.setdiff1d(T[t], np.unique(np.concatenate([T[o] for o in taxa if o != t]))) for t in taxa}
    allteo = np.unique(np.concatenate(list(T.values())))
    print("B73 hashes", len(b73), "| non-B73 per assembly:", {n: len(v) for n, v in asm.items()})
    print("taxon-specific hashes:", {t: len(U[t]) for t in taxa})

    info = {}
    with open(a.pools_tsv) as f:
        for r in csv.reader(f, delimiter="\t"):
            if r and r[0] != "SAMPLE": info[r[0]] = (r[1], r[2] if len(r) > 2 else "", r[3] if len(r) > 3 else "")
    cols = ["sample", "role", "label", "accession", "n_hashes", "f_teo"] + [f"f_{t}" for t in taxa] + [f"n_{t}" for t in taxa] + \
           [f"c_{n}" for n in asm] + ["assigned", "ratio", "agrees"]
    rows = []
    for p in sorted(glob.glob(os.path.join(a.pools, "*.sig.gz"))):
        s = os.path.basename(p).replace(".sig.gz", "")
        H = hashes(p)
        n = {t: int(np.isin(U[t], H, assume_unique=True).sum()) for t in taxa}
        f = {t: n[t] / max(len(U[t]), 1) for t in taxa}
        c = {k: np.isin(v, H, assume_unique=True).sum() / max(len(v), 1) for k, v in asm.items()}
        order = sorted(taxa, key=lambda t: -f[t]); best, second = order[0], order[1]
        role, label, acc = info.get(s, ("", "", ""))
        lab_tax = label
        rows.append([s, role, label, acc, len(H), np.isin(allteo, H, assume_unique=True).sum() / len(allteo)] +
                    [f[t] for t in taxa] + [n[t] for t in taxa] + [c[k] for k in asm] +
                    [best, f[best] / f[second] if f[second] > 0 else float("inf"),
                     "" if role != "donor" else str(lab_tax == best)])
    with open(a.out, "w") as fo:
        w = csv.writer(fo, delimiter="\t"); w.writerow(cols)
        for r in rows: w.writerow([f"{v:.6g}" if isinstance(v, float) else v for v in r])
    d = [r for r in rows if r[1] == "donor"]
    agree = sum(r[-1] == "True" for r in d)
    print(f"pools: {len(rows)} | donor pools whose k-mer taxon = pedigree taxon: {agree} / {len(d)}")
    print("WROTE", a.out)


if __name__ == "__main__":
    main()
