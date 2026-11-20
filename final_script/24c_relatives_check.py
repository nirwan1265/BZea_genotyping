#!/usr/bin/env python3
"""
24c — Relatives check: are introgression tracts inherited?
For every SEG or DONOR tract of a pool, take the other pools with a given relationship to it (read from the
pedigree name, e.g. Zv.0200_P3_P2_P1.1.1.1 = accession Zv.0200, F1 plant P3, BC1 plant P2, BC2 plant P1) and
count the share of them that also carry an introgression (SEG or DONOR) somewhere inside the tract (100-kb grid). Averaged over
all tracts, separately for short (< 10 Mb) and long (>= 10 Mb) tracts. Donor pools with FLAG ok only.
Usage: python final_script/24c_relatives_check.py --base data/introgression_fix1 \
          --names data/phenotypes/raw/name_change_final.txt [--out <tsv>]
Needs: numpy, pandas.
"""
import argparse
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--base", required=True, help="folder with pool_summary.tsv and tracts_all.tsv.gz")
ap.add_argument("--names", required=True, help="Sequencing_ID -> new_genotype (pedigree name)")
ap.add_argument("--short", type=float, default=10.0, help="short/long cut, Mb")
ap.add_argument("--out", default="")
a = ap.parse_args()

ps = pd.read_csv(f"{a.base}/pool_summary.tsv", sep="\t")
ok = set(ps.loc[(ps.ROLE == "donor") & (ps.FLAG == "ok"), "SAMPLE"])
nm = pd.read_csv(a.names, sep="\t", dtype=str).drop_duplicates("Sequencing_ID")
nm = nm[nm.Sequencing_ID.isin(ok)].reset_index(drop=True)
p = nm.new_genotype.str.split(r"[._]")
for k, n in [("TAX", 1), ("ACC", 2), ("F1", 3), ("BC1", 4), ("BC2", 5)]:
    nm[k] = p.str[:n].str.join("_")
t = pd.read_csv(f"{a.base}/tracts_all.tsv.gz", sep="\t")
t = t[t.SAMPLE.isin(nm.Sequencing_ID)]
intro = t[t.STATE != "B73"]
pools = nm.Sequencing_ID.tolist()
pid = {s: i for i, s in enumerate(pools)}
# carrier matrix on a 100-kb grid: C[pool, bin] = pool has SEG or DONOR there
BIN = 100_000
off, o = {}, 0
for c, g in t.groupby("CHROM"):
    off[c] = o
    o += int(g.END.max() // BIN) + 1
C = np.zeros((len(pools), o), bool)
for r in intro.itertuples():
    C[pid[r.SAMPLE], off[r.CHROM] + int(r.START // BIN): off[r.CHROM] + int(np.ceil(r.END / BIN))] = True
lab = {k: nm[k].to_numpy() for k in ("TAX", "ACC", "F1", "BC1", "BC2")}

# groups of other pools, per relationship (each other pool counts once, in its closest relationship)
REL = [("Same taxon, other accession (unrelated: the baseline)", "TAX", "ACC"),
       ("Same accession, other F1 plant", "ACC", "F1"),
       ("Same F1 plant, other BC1 plant", "F1", "BC1"),
       ("Same BC1 plant, other BC2 plant", "BC1", "BC2"),
       ("Same BC2 plant", "BC2", None)]
rows = []
for r in intro.itertuples():
    i = pid[r.SAMPLE]
    b0, b1 = off[r.CHROM] + int(r.START // BIN), off[r.CHROM] + int(np.ceil(r.END / BIN))
    L = (r.END - r.START) / 1e6
    for name, same, diff in REL:
        m = lab[same] == lab[same][i]
        m = m & (lab[diff] != lab[diff][i]) if diff else m
        m[i] = False
        if m.any():
            rows.append((name, "short" if L < a.short else "long", C[m, b0:b1].any(1).mean()))
D = pd.DataFrame(rows, columns=["relationship", "length", "share"])
T = D.pivot_table(index="relationship", columns="length", values="share", aggfunc="mean").reindex([x[0] for x in REL])
N = D.pivot_table(index="relationship", columns="length", values="share", aggfunc="size").reindex([x[0] for x in REL])
out = pd.DataFrame({"short_pct": (100 * T["short"]).round(1), "long_pct": (100 * T["long"]).round(1),
                    "short_tracts": N["short"], "long_tracts": N["long"]})
print(f"{len(nm)} pools, {len(intro)} SEG/DONOR tracts")
print(out.to_string())
if a.out:
    import os; os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    out.to_csv(a.out, sep="\t")
