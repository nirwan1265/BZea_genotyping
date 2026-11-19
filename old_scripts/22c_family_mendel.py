#!/usr/bin/env python3
"""
22c — Mendelian test within sister-line families (README section 2.8.1).

Pools are named accession_F1_BC1_BC2.S1.S2.S3 (name_change_final.txt). Sister lines share a BC1 plant
(BC1 family) or a BC2 plant (BC2 family). At a locus a BC1 or BC2 plant is AB or AA (backcross to B73).

BC2 family, BC2 plant AB: each pool's mother (BC2S2) is two selfings below it, so
    pool state B73 : SEG : DONOR = 3/8 : 1/4 : 3/8 (Mendel).
BC1 family, BC1 plant AB: each BC2 plant is AB with probability 1/2, then as above.
If any pool of a family carries donor at a locus, its BC1 (or BC2) plant was AB there. Families whose plant
was AB but whose pools all look B73 are missed; the expected counts are conditioned on "at least one pool
carries donor" so this does not bias the test.

Outputs, per family level (BC2, BC1):
  1. Within families that carry donor: observed vs expected numbers of B73, SEG and DONOR pools.
     SEG below expectation -> calls or selection against AB; DONOR below -> selection against BB (or calls).
  2. Share of the genome where the family plant was AB, corrected for missed families,
     against Mendel (BC2 plant 25%, BC1 plant 50%) -> loss during the backcrosses.
CIs by bootstrap over families. Loci: one window every --step_mb (default 1 Mb).

Usage: python final_script/22c_family_mendel.py --tracts tracts_all.tsv.gz --pool_summary pool_summary.tsv \
          --names data/phenotypes/raw/name_change_final.txt [--step_mb 1] [--boot 500] [--out family_mendel.tsv]
Needs: numpy, pandas.
"""
import argparse
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--tracts", required=True)
ap.add_argument("--pool_summary", required=True)
ap.add_argument("--names", required=True)
ap.add_argument("--step_mb", type=float, default=1.0)
ap.add_argument("--boot", type=int, default=500)
ap.add_argument("--seed", type=int, default=1)
ap.add_argument("--out", default="family_mendel.tsv")
a = ap.parse_args()
rng = np.random.default_rng(a.seed)
P_SELF2 = np.array([3 / 8, 1 / 4, 3 / 8])          # B73, SEG, DONOR below an AB plant after two selfings

ps = pd.read_csv(a.pool_summary, sep="\t")
ok = set(ps.loc[(ps.ROLE == "donor") & (ps.FLAG == "ok"), "SAMPLE"])
nm = pd.read_csv(a.names, sep="\t", dtype=str).drop_duplicates("Sequencing_ID")
nm = nm[nm.Sequencing_ID.isin(ok)]
part = nm.new_genotype.str.split(r"[._]")
nm["BC1"] = part.str[:4].str.join("_")              # accession (2 parts) + F1 + BC1
nm["BC2"] = part.str[:5].str.join("_")              # ... + BC2
t = pd.read_csv(a.tracts, sep="\t")
t = t[t.SAMPLE.isin(nm.Sequencing_ID)]
code = {"B73": 0, "SEG": 1, "DONOR": 2}

# state matrix: pools x loci
loci = []
for c, g in t.groupby("CHROM"):
    for p in np.arange(a.step_mb / 2 * 1e6, g.END.max(), a.step_mb * 1e6):
        loci.append((c, p))
pools = nm.Sequencing_ID.tolist(); pidx = {s: i for i, s in enumerate(pools)}
S = np.full((len(pools), len(loci)), -1, np.int8)
for c, g in t.groupby("CHROM"):
    lj = [j for j, (cc, _) in enumerate(loci) if cc == c]; pos = np.array([loci[j][1] for j in lj])
    for r in g.itertuples():
        hit = (pos >= r.START) & (pos < r.END)
        S[pidx[r.SAMPLE], np.array(lj)[hit]] = code[r.STATE]
print(f"{len(pools)} pools, {len(loci)} loci (every {a.step_mb:g} Mb)")


def fam_stats(rows_by_sub):
    """one family; rows_by_sub = list of arrays of pool rows, one per BC2 plant (one entry for a BC2 family).
    Returns observed counts, expected counts (conditioned on >=1 donor pool), loci with donor, expected loci
    with donor per AB locus (detection probability), number of usable loci."""
    X = np.vstack([S[r] for r in rows_by_sub]); use = (X >= 0).all(axis=0)
    X = X[:, use]
    carry = (X > 0).any(axis=0)
    obs = np.array([(X[:, carry] == k).sum() for k in range(3)], float)
    # P(no donor pool | family plant AB) and expected counts given >=1 donor pool
    ns = [len(r) for r in rows_by_sub]
    if len(ns) == 1:                                 # BC2 family: plant AB
        n = ns[0]; p0 = P_SELF2[0] ** n
        e_all = n * P_SELF2                          # expected counts given AB
        e_none = np.array([n, 0, 0]) * p0            # contribution of the all-B73 outcome
        det = 1 - p0
    else:                                            # BC1 family: each BC2 plant AB with prob 1/2
        p0 = np.prod([0.5 + 0.5 * P_SELF2[0] ** n for n in ns])
        e_all = sum(n * (0.5 * P_SELF2 + 0.5 * np.array([1, 0, 0])) for n in ns)
        e_none = np.array([sum(ns), 0, 0]) * p0
        det = 1 - p0
    exp_per_locus = (e_all - e_none) / det           # expected counts per carrying locus
    return obs, exp_per_locus * carry.sum(), carry.sum(), det, use.sum()


rows = []
for level, keycol, ab_prior in [("BC2 family", "BC2", 0.25), ("BC1 family", "BC1", 0.50)]:
    fams = []
    for f, g in nm.groupby(keycol):
        if len(g) < 2: continue
        idx = [pidx[s] for s in g.Sequencing_ID]
        subs = [np.array(idx)] if keycol == "BC2" else [np.array([pidx[s] for s in h.Sequencing_ID]) for _, h in g.groupby("BC2")]
        o, e, ncar, det, nuse = fam_stats(subs)
        fams.append((o, e, ncar, det, nuse, len(g)))
    O = np.array([x[0] for x in fams]); E = np.array([x[1] for x in fams])
    car = np.array([x[2] for x in fams]); det = np.array([x[3] for x in fams]); nuse = np.array([x[4] for x in fams])
    ab = car / det / nuse                            # share of loci where the family plant was AB, corrected

    def summ(ix):
        o, e = O[ix].sum(0), E[ix].sum(0)
        return o / o.sum(), e / e.sum(), (o / e), np.average(ab[ix], weights=nuse[ix])
    of, ef, ratio, abm = summ(np.arange(len(fams)))
    B = [summ(rng.integers(0, len(fams), len(fams))) for _ in range(a.boot)]
    ci = lambda k, i: np.percentile([b[k][i] if np.ndim(b[k]) else b[k] for b in B], [2.5, 97.5])
    for i, st in enumerate(["B73", "SEG", "DONOR"]):
        lo, hi = ci(2, i)
        rows.append({"level": level, "families": len(fams), "pools": int(sum(x[5] for x in fams)), "measure": f"{st} pools in carrying families",
                     "observed_share": round(of[i], 4), "expected_share": round(ef[i], 4), "obs_over_exp": round(ratio[i], 3),
                     "CI95": f"{lo:.3f}-{hi:.3f}"})
    lo, hi = ci(3, 0)
    rows.append({"level": level, "families": len(fams), "pools": int(sum(x[5] for x in fams)), "measure": "genome share where family plant was AB",
                 "observed_share": round(abm, 4), "expected_share": ab_prior, "obs_over_exp": round(abm / ab_prior, 3),
                 "CI95": f"{lo / ab_prior:.3f}-{hi / ab_prior:.3f}"})
R = pd.DataFrame(rows)
R.to_csv(a.out, sep="\t", index=False)
print(R.to_string(index=False))
