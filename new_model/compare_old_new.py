#!/usr/bin/env python3
"""
new_model / compare — current calls (assembly DONOR levels) against the new model (levels from all pools, 20c).

1. Genome shares per taxon: SEG and DONOR (% of the genome), old vs new.
2. SEG tracts between two DONOR tracts (the error found with the targeted sequencing): % of the genome, old vs new,
   and the hotspot bins (pools whose sandwiched SEG covers the 1-Mb bin; chr6 21-23, chr1 297, chr7 54-65, chr3 75-78 Mb).
3. With --targseq_old/--targseq_new (34d tracts_targseq.tsv of each call set): the validation table
   (targeted-sequencing genotype by our call) and the BB share of sandwiched SEG tracts, old vs new.
Usage: python new_model/compare_old_new.py --old $D21/tracts_all.tsv.gz --new $DNEW/introgression/tracts_all.tsv.gz \
          --pool_summary $D21/pool_summary.tsv [--targseq_old A/tracts_targseq.tsv --targseq_new B/tracts_targseq.tsv]
Needs: numpy, pandas.
"""
import argparse
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--old", required=True); ap.add_argument("--new", required=True)
ap.add_argument("--pool_summary", required=True)
ap.add_argument("--targseq_old", default=""); ap.add_argument("--targseq_new", default="")
a = ap.parse_args()
ps = pd.read_csv(a.pool_summary, sep="\t")
ok = set(ps.loc[(ps.ROLE == "donor") & (ps.FLAG == "ok"), "SAMPLE"])
BAD = {"Zv0480_P3_P2_P1211", "Zx0090_P2_P5_P1111", "Zx0160_P2_P4_P1111"}
HOT = [("chr6", 21, 23), ("chr1", 297, 297), ("chr7", 54, 65), ("chr3", 75, 78)]


def load(p):
    t = pd.read_csv(p, sep="\t"); t = t[t.SAMPLE.isin(ok)].sort_values(["SAMPLE", "CHROM", "START"]).copy()
    t["Mb"] = (t.END - t.START) / 1e6
    t["prev"] = t.groupby(["SAMPLE", "CHROM"]).STATE.shift(1); t["next"] = t.groupby(["SAMPLE", "CHROM"]).STATE.shift(-1)
    t["sand"] = (t.STATE == "SEG") & (t.prev == "DONOR") & (t.next == "DONOR")
    return t


O, N = load(a.old), load(a.new)
rows = []
for tx in sorted(set(O.GROUP) | {"All"}):
    r = {"taxon": tx}
    for lab, T in [("old", O), ("new", N)]:
        x = T if tx == "All" else T[T.GROUP == tx]; tot = x.Mb.sum()
        r[f"SEG_{lab}"] = round(100 * x.Mb[x.STATE == "SEG"].sum() / tot, 2)
        r[f"DONOR_{lab}"] = round(100 * x.Mb[x.STATE == "DONOR"].sum() / tot, 2)
        r[f"sandwichedSEG_{lab}"] = round(100 * x.Mb[x.sand].sum() / tot, 2)
    rows.append(r)
print("1-2. Genome shares (%), old = assembly levels, new = levels from all pools:")
print(pd.DataFrame(rows).to_string(index=False))
print("\nHotspot bins: pools with a SEG tract between two DONOR tracts covering the bin (old -> new)")
for c, s, e in HOT:
    def n(T):
        x = T[T.sand & (T.CHROM == c)]
        return max(int(((x.START <= b * 1e6 + 5e5) & (x.END > b * 1e6 + 5e5)).sum()) for b in range(s, e + 1))
    print(f"  {c} {s}-{e} Mb: {n(O)} -> {n(N)}")

if a.targseq_old and a.targseq_new:
    print("\n3. Targeted sequencing (34d), share of clearly called tracts:")
    for lab, p, T in [("old", a.targseq_old, O), ("new", a.targseq_new, N)]:
        R = pd.read_csv(p, sep="\t"); R = R[~R["sample"].isin(BAD) & R.plant_genotype.isin(["AA", "AB", "BB"])]
        tab = pd.crosstab(R.our_call, R.plant_genotype, normalize="index").reindex(["B73", "SEG", "DONOR"]).round(3)
        tab["tracts"] = R.our_call.value_counts().reindex(tab.index)
        m = R[R.our_call == "SEG"].merge(T[["SAMPLE", "CHROM", "START", "sand"]], left_on=["pool", "chr", "start"],
                                         right_on=["SAMPLE", "CHROM", "START"], how="left")
        sw = m[m.sand == True]
        print(f"--- {lab}\n{tab.to_string()}\n  SEG between two DONOR tracts: {len(sw)} tracts, BB {100 * (sw.plant_genotype == 'BB').mean():.0f}%"
              if len(sw) else f"--- {lab}\n{tab.to_string()}\n  SEG between two DONOR tracts: none")
