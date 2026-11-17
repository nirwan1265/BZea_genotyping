#!/usr/bin/env python3
"""
Tract lengths: simulated BC2S2 mothers vs the k-mer HMM calls (README section 2.6.5).

Simulation: one donor chromosome through F1 -> BC1 -> BC2 (backcross to B73), then two selfings
(BC2S1, BC2S2). Crossovers are placed at random along the NAM genetic map (Poisson, Haldane, no
interference, no selection). The BC2S2 plant is the mother of a BC2S3 pool, so its genotype sets the
pool state. Introgressed tract = run of windows with >= 1 donor copy (SEG and DONOR joined).

Observed: tracts_all.tsv.gz (step 21), donor pools with FLAG == ok in pool_summary.tsv; adjacent SEG and
DONOR tracts are joined the same way.

Usage:
  python final_script/tract_length_sim.py --map_dir genetic_maps --tracts tracts_all.tsv.gz \
      --pool_summary pool_summary.tsv [--reps 400] [--seed 1]
  map_dir: <chr>.map per chromosome (chr1.map ...), two columns: bp, cM (no header), as used by step 20.
"""
import argparse
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--map_dir", required=True)
ap.add_argument("--tracts", required=True)
ap.add_argument("--pool_summary", required=True)
ap.add_argument("--reps", type=int, default=400, help="simulated lines per chromosome")
ap.add_argument("--step", type=float, default=20000, help="simulation grid in bp")
ap.add_argument("--seed", type=int, default=1)
a = ap.parse_args()
rng = np.random.default_rng(a.seed)
CUTS = [0.5, 1, 2, 5]


def summary(L, label):
    L = np.asarray(L)
    row = {"set": label, "n": len(L), "median_Mb": round(float(np.median(L)), 1)}
    for c in CUTS:
        row[f"<{c}Mb_%"] = round(100 * float(np.mean(L < c)), 1)
    row["bp_in_<1Mb_%"] = round(100 * L[L < 1].sum() / L.sum(), 2)
    return row


def runs(mask, pos_start, pos_end):
    """lengths (Mb) of runs of True"""
    e = np.flatnonzero(np.diff(np.r_[0, mask.astype(int), 0]))
    return [(pos_end[b - 1] - pos_start[a_]) / 1e6 for a_, b in zip(e[::2], e[1::2])]


# ---- observed ----
t = pd.read_csv(a.tracts, sep="\t")
ps = pd.read_csv(a.pool_summary, sep="\t")
ok = set(ps.loc[(ps.ROLE == "donor") & (ps.FLAG == "ok"), "SAMPLE"])
t = t[t.SAMPLE.isin(ok)].sort_values(["SAMPLE", "CHROM", "START"])
obs = []
for _, g in t.groupby(["SAMPLE", "CHROM"]):
    obs += runs((g.STATE != "B73").to_numpy(), g.START.to_numpy(), g.END.to_numpy())


# ---- simulation ----
def meiosis(h1, h2, cm, total_cM):
    """one gamete from a plant with haplotypes h1, h2 (arrays of donor copies 0/1 on the grid)"""
    xs = rng.uniform(0, total_cM, rng.poisson(total_cM / 100))
    side = (np.searchsorted(np.sort(xs), cm, side="right") + rng.integers(2)) % 2
    return np.where(side == 0, h1, h2)


sim = []
for c in range(1, 11):
    m = pd.read_csv(f"{a.map_dir}/chr{c}.map", sep=r"\s+", header=None, names=["bp", "cM"]).sort_values("bp")
    m["cM"] = np.maximum.accumulate(m.cM.to_numpy())
    start = np.arange(0, m.bp.max(), a.step); end = start + a.step
    cm = np.interp(start + a.step / 2, m.bp, m.cM); M = float(cm[-1])
    D = np.ones(len(start), int); B = np.zeros(len(start), int)
    for _ in range(a.reps):
        bc1 = meiosis(D, B, cm, M)                       # F1 gamete (+ B73 gamete) = BC1
        bc2 = meiosis(bc1, B, cm, M)                     # BC1 gamete = donor haplotype of BC2
        s1 = (meiosis(bc2, B, cm, M), meiosis(bc2, B, cm, M))
        s2 = (meiosis(*s1, cm, M), meiosis(*s1, cm, M))  # BC2S2 = mother of the BC2S3 pool
        sim += runs((s2[0] + s2[1]) > 0, start, end)

out = pd.DataFrame([summary(sim, "simulated BC2S2 mothers"), summary(obs, "k-mer HMM calls")])
print(out.to_string(index=False))
