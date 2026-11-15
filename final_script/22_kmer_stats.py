#!/usr/bin/env python3
"""
22 — Observed vs expected BC2S3 proportions, with a family bootstrap.

Pool states (k-mer HMM, step 21) -> expected plant genotypes in the pool:
  a SEG region is a 1:2:1 family, so its plants are 1/4 AA, 1/2 AB, 1/4 BB
  AA = B73 + SEG/4      AB = SEG/2      BB = DONOR + SEG/4      donor alleles = BB + AB/2

Expected for a BC2S3 family (pool state = genotype of its BC2S2 parent):
  SEG 6.25%, DONOR 9.375%, donor alleles 12.5%  ->  plants AB 3.125%, BB 10.9375%
  SEG/(SEG+DONOR) = 0.40, AB/(AB+BB) = 0.222 (share of donor-carrying plants that are het)

Test: pools of one family share their donor segments, so the unit is the family. Families are
resampled with replacement (with all their pools) NBOOT times; 95% CI = 2.5-97.5 percentiles;
two-sided p = 2 * min(P(boot <= expected), P(boot >= expected)), floor 1/NBOOT.
Only non-noisy donor pools (FLAG ok); Zm is not in the k-mer HMM.

Outputs in <out_dir>: stats_overall.tsv (all taxa and per taxon), stats_by_chrom.tsv (all taxa).
"""
import argparse
import collections
import csv
import gzip
import os

import numpy as np

EXPECTED = {"SEG": 0.0625, "DONOR": 0.09375, "DONOR_ALLELES": 0.125, "SEG_OF_INTROGRESSED": 0.40,
            "AB": 0.03125, "BB": 0.109375, "AB_OF_DONOR_PLANTS": 0.03125 / (0.03125 + 0.109375)}
LABEL = {"SEG": "SEG (pool)", "DONOR": "DONOR (pool)", "DONOR_ALLELES": "donor alleles",
         "SEG_OF_INTROGRESSED": "SEG/(SEG+DONOR)", "AB": "AB (plant)", "BB": "BB (plant)",
         "AB_OF_DONOR_PLANTS": "AB/(AB+BB)"}


def measures(pairs):
    """pairs: list of (SEG, DONOR) fractions, one per pool."""
    v = np.asarray(pairs); s, d = v[:, 0].mean(), v[:, 1].mean()
    ab, bb = s / 2, d + s / 4
    return {"SEG": s, "DONOR": d, "DONOR_ALLELES": d + s / 2, "SEG_OF_INTROGRESSED": s / (s + d),
            "AB": ab, "BB": bb, "AB_OF_DONOR_PLANTS": ab / (ab + bb)}


def bootstrap(fams, nboot, rng):
    """fams: list of lists of (SEG, DONOR). Returns {measure: (obs, lo, hi, p)}."""
    obs = measures([x for f in fams for x in f])
    boots = collections.defaultdict(list)
    for _ in range(nboot):
        pick = rng.integers(0, len(fams), len(fams))
        for k, v in measures([x for i in pick for x in fams[i]]).items(): boots[k].append(v)
    out = {}
    for k, o in obs.items():
        a = np.array(boots[k]); lo, hi = np.percentile(a, [2.5, 97.5])
        p = max(2 * min(np.mean(a <= EXPECTED[k]), np.mean(a >= EXPECTED[k])), 1.0 / nboot)
        out[k] = (o, lo, hi, min(p, 1.0))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pool_summary", required=True, help="step 21 pool_summary.tsv")
    ap.add_argument("--tracts", required=True, help="step 21 tracts_all.tsv.gz")
    ap.add_argument("--meta", required=True, help="SAMPLE_META (col 1 sample, col 2 pedigree)")
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--nboot", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    rng = np.random.default_rng(a.seed)

    fam = {}
    with open(a.meta, encoding="utf-8", errors="ignore") as f:
        next(f)
        for line in f:
            p = line.rstrip("\r\n").split("\t")
            if len(p) > 1: fam[p[0]] = "_".join(p[1].split("_")[:2])
    pools = {}
    with open(a.pool_summary) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["ROLE"] == "donor" and r["FLAG"] == "ok":
                pools[r["SAMPLE"]] = (r["GROUP"], fam.get(r["SAMPLE"], "unknown"), float(r["FRAC_SEG"]), float(r["FRAC_DONOR"]))
    print(f"non-noisy donor pools: {len(pools)} | families: {len({(g, fm) for g, fm, _, _ in pools.values()})}")

    def fam_lists(sel):
        by = collections.defaultdict(list)
        for s, (g, fm, seg, don) in pools.items():
            if sel(s, g): by[(g, fm)].append((seg, don))
        return list(by.values())

    rows = []
    for scope in ["ALL"] + sorted({g for g, *_ in pools.values()}):
        fl = fam_lists(lambda s, g: scope == "ALL" or g == scope)
        res = bootstrap(fl, a.nboot, rng)
        print(f"\n== {scope}: {len(fl)} families, {sum(map(len, fl))} pools")
        for k, (o, lo, hi, p) in res.items():
            print(f"  {LABEL[k]:18s} {100*o:6.2f}%  95% CI {100*lo:6.2f}-{100*hi:6.2f}  expected {100*EXPECTED[k]:6.2f}%  p {p:.4f}")
            rows.append((scope, len(fl), sum(map(len, fl)), k, o, lo, hi, EXPECTED[k], p))
    with open(f"{a.out_dir}/stats_overall.tsv", "w") as f:
        f.write("SCOPE\tFAMILIES\tPOOLS\tMEASURE\tOBSERVED\tCI_LOW\tCI_HIGH\tEXPECTED\tP\n")
        for r in rows: f.write("\t".join(str(x) if i < 4 else f"{x:.6f}" for i, x in enumerate(r)) + "\n")

    # per chromosome (all taxa): fractions from the tracts
    clen = collections.defaultdict(int); bp = collections.defaultdict(lambda: [0, 0])
    with gzip.open(a.tracts, "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            e = int(r["END"]); clen[r["CHROM"]] = max(clen[r["CHROM"]], e)
            if r["SAMPLE"] in pools and r["STATE"] in ("SEG", "DONOR"):
                bp[(r["SAMPLE"], r["CHROM"])][0 if r["STATE"] == "SEG" else 1] += e - int(r["START"])
    crow = []
    print("\n== per chromosome (all taxa)\n  chrom  SEG%   DONOR%  AB/(AB+BB)  95% CI        p")
    for c in sorted(clen, key=lambda x: int(x.replace("chr", "")) if x.replace("chr", "").isdigit() else 99):
        by = collections.defaultdict(list)
        for s, (g, fm, _, _) in pools.items():
            sg, dn = bp.get((s, c), (0, 0)); by[(g, fm)].append((sg / clen[c], dn / clen[c]))
        res = bootstrap(list(by.values()), a.nboot, rng)
        o, lo, hi, p = res["AB_OF_DONOR_PLANTS"]
        print(f"  {c:6s} {100*res['SEG'][0]:5.2f}  {100*res['DONOR'][0]:6.2f}  {100*o:9.1f}%  {100*lo:5.1f}-{100*hi:5.1f}  {p:.4f}")
        for k, (o, lo, hi, p) in res.items(): crow.append((c, k, o, lo, hi, EXPECTED[k], p))
    with open(f"{a.out_dir}/stats_by_chrom.tsv", "w") as f:
        f.write("CHROM\tMEASURE\tOBSERVED\tCI_LOW\tCI_HIGH\tEXPECTED\tP\n")
        for r in crow: f.write("\t".join(str(x) if i < 2 else f"{x:.6f}" for i, x in enumerate(r)) + "\n")
    print(f"\nWROTE: {a.out_dir}/stats_overall.tsv, stats_by_chrom.tsv")


if __name__ == "__main__":
    main()
