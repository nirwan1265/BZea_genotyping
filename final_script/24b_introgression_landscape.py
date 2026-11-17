#!/usr/bin/env python3
"""
24b: where in the genome the introgressions are, and how many and how long they are (README section 2.8).

From the k-mer HMM calls (step 21, tracts_all.tsv.gz), for the donor pools with FLAG == ok (pool_summary.tsv):
  1. Each pool is put on a 100-kb grid: B73 = 0, SEG = 0.5, DONOR = 1 donor allele share
     (SEG = half the plants' alleles are donor on average: 1/4 AA, 1/2 AB, 1/4 BB).
  2. Landscape: the mean donor allele share over all pools, per bin (--bin_mb, default 1 Mb),
     drawn as a heatmap with one row per chromosome (viridis).
  3. The same per taxon, along the whole genome in coarser bins (--taxon_bin_mb, default 5 Mb).
  4. Tract statistics per taxon: number of SEG and DONOR tracts (total and per pool) and their median length;
     also introgressed tracts = adjacent SEG and DONOR joined (as in section 2.6.5).

Usage:
  python final_script/24b_introgression_landscape.py --tracts tracts_all.tsv.gz --pool_summary pool_summary.tsv \
      --out final_script/figures/results_spread [--bin_mb 1] [--taxon_bin_mb 5]
Outputs (in --out): introgression_landscape.png, introgression_landscape_by_taxon.png,
  introgression_landscape_bins.tsv, tract_stats_by_taxon.tsv
Needs: numpy, pandas, matplotlib.
"""
import argparse
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--tracts", required=True)
ap.add_argument("--pool_summary", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--bin_mb", type=float, default=1.0)
ap.add_argument("--taxon_bin_mb", type=float, default=5.0)
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)

STEP = 100_000
TAXA = ["Zd", "Zh", "Zl", "Zv", "Zx"]
NAMES = {"Zd": "diploperennis", "Zh": "huehuetenangensis", "Zl": "luxurians", "Zv": "parviglumis", "Zx": "mexicana"}
CHROMS = [f"chr{i}" for i in range(1, 11)]
VAL = {"B73": 0.0, "SEG": 0.5, "DONOR": 1.0}

t = pd.read_csv(a.tracts, sep="\t")
ps = pd.read_csv(a.pool_summary, sep="\t")
ok = set(ps.loc[(ps.ROLE == "donor") & (ps.FLAG == "ok"), "SAMPLE"])
t = t[t.SAMPLE.isin(ok) & t.GROUP.isin(TAXA)].sort_values(["SAMPLE", "CHROM", "START"])
pools = t.drop_duplicates("SAMPLE")[["SAMPLE", "GROUP"]].reset_index(drop=True)
idx = {s: i for i, s in enumerate(pools.SAMPLE)}
print(f"{len(pools)} donor pools")

# ---- 1. 100-kb grid per pool and chromosome
chrlen = t.groupby("CHROM").END.max()
grid = {}
for c in CHROMS:
    n = int(np.ceil(chrlen[c] / STEP))
    g = np.full((len(pools), n), np.nan)
    for r in t[t.CHROM == c].itertuples():
        g[idx[r.SAMPLE], int(r.START // STEP):int(np.ceil(r.END / STEP))] = VAL[r.STATE]
    grid[c] = g


def binned(g, rows, mb):
    """mean donor allele share over the given pools, in bins of mb Mb"""
    k = int(round(mb * 1e6 / STEP))
    m = np.nanmean(g[rows], axis=0)
    n = int(np.ceil(len(m) / k))
    m = np.pad(m, (0, n * k - len(m)), constant_values=np.nan).reshape(n, k)
    return np.nanmean(m, axis=1)


# ---- 2. landscape, all pools, one row per chromosome
allrows = np.arange(len(pools))
L = {c: binned(grid[c], allrows, a.bin_mb) for c in CHROMS}
nmax = max(len(v) for v in L.values())
M = np.full((10, nmax), np.nan)
for i, c in enumerate(CHROMS):
    M[i, :len(L[c])] = 100 * L[c]
cmap = plt.get_cmap("viridis").copy(); cmap.set_bad("white")
fig, ax = plt.subplots(figsize=(12, 4.2))
im = ax.imshow(np.ma.masked_invalid(M), aspect="auto", cmap=cmap, interpolation="none",
               extent=[0, nmax * a.bin_mb, 10.5, 0.5])
ax.set_yticks(range(1, 11)); ax.set_yticklabels(CHROMS)
ax.set_xlabel("Position (Mb, B73 v5)")
ax.set_title(f"Donor allele share across {len(pools)} donor pools, {a.bin_mb:g}-Mb bins (expected for BC2S3: 12.5%)",
             fontsize=11, loc="left")
cb = fig.colorbar(im, ax=ax, pad=0.01); cb.set_label("Donor alleles (%)")
fig.tight_layout(); fig.savefig(os.path.join(a.out, "introgression_landscape.png"), dpi=200); plt.close(fig)

# ---- 3. per taxon, whole genome
rowsby = {x: np.flatnonzero(pools.GROUP.values == x) for x in TAXA}
T, edges, pos = [], [0], []
for c in CHROMS:
    cols = [binned(grid[c], rowsby[x], a.taxon_bin_mb) for x in TAXA]
    T.append(np.vstack(cols)); edges.append(edges[-1] + cols[0].size)
T = 100 * np.hstack(T)
fig, ax = plt.subplots(figsize=(13, 3.2))
im = ax.imshow(np.ma.masked_invalid(T), aspect="auto", cmap=cmap, interpolation="none")
for e in edges[1:-1]:
    ax.axvline(e - 0.5, color="white", lw=1.2)
ax.set_xticks([(edges[i] + edges[i + 1]) / 2 - 0.5 for i in range(10)]); ax.set_xticklabels(CHROMS)
ax.set_yticks(range(len(TAXA)))
ax.set_yticklabels([f"{x} {NAMES[x]} ({len(rowsby[x])})" for x in TAXA])
ax.set_title(f"Donor allele share per taxon, {a.taxon_bin_mb:g}-Mb bins (number of pools in brackets)", fontsize=11, loc="left")
cb = fig.colorbar(im, ax=ax, pad=0.01); cb.set_label("Donor alleles (%)")
fig.tight_layout(); fig.savefig(os.path.join(a.out, "introgression_landscape_by_taxon.png"), dpi=200); plt.close(fig)

# bins table (bin_mb), all pools and per taxon
rows = []
for c in CHROMS:
    per = {x: binned(grid[c], rowsby[x], a.bin_mb) for x in TAXA}
    for i, v in enumerate(L[c]):
        rows.append({"chr": c, "bin_start_Mb": i * a.bin_mb, "all_pct": round(100 * v, 2),
                     **{f"{x}_pct": round(100 * per[x][i], 2) for x in TAXA}})
B = pd.DataFrame(rows)
B.to_csv(os.path.join(a.out, "introgression_landscape_bins.tsv"), sep="\t", index=False)
q = B.all_pct.quantile([0.05, 0.5, 0.95]).round(1).to_dict()
print("bin donor allele share (%): 5th / median / 95th percentile", q)
print("highest bins:\n", B.sort_values("all_pct", ascending=False).head(8)[["chr", "bin_start_Mb", "all_pct"]].to_string(index=False))
print("lowest bins:\n", B.sort_values("all_pct").head(8)[["chr", "bin_start_Mb", "all_pct"]].to_string(index=False))

# ---- 4. tract statistics per taxon
t["L"] = (t.END - t.START) / 1e6
t["intro"] = t.STATE != "B73"
t["run"] = (t.intro != t.groupby(["SAMPLE", "CHROM"]).intro.shift()).cumsum()
joined = t[t.intro].groupby(["SAMPLE", "CHROM", "run", "GROUP"]).agg(start=("START", "min"), end=("END", "max")).reset_index()
joined["L"] = (joined.end - joined.start) / 1e6


def stats(tt, jj, n):
    s, d = tt[tt.STATE == "SEG"], tt[tt.STATE == "DONOR"]
    return {"pools": n,
            "SEG_tracts": len(s), "SEG_per_pool": round(len(s) / n, 1), "SEG_median_Mb": round(s.L.median(), 1),
            "DONOR_tracts": len(d), "DONOR_per_pool": round(len(d) / n, 1), "DONOR_median_Mb": round(d.L.median(), 1),
            "introgressed_tracts_per_pool": round(len(jj) / n, 1), "introgressed_median_Mb": round(jj.L.median(), 1)}


S = [{"taxon": x, **stats(t[t.GROUP == x], joined[joined.GROUP == x], len(rowsby[x]))} for x in TAXA]
S.append({"taxon": "All", **stats(t, joined, len(pools))})
S = pd.DataFrame(S)
S.to_csv(os.path.join(a.out, "tract_stats_by_taxon.tsv"), sep="\t", index=False)
print(S.to_string(index=False))
