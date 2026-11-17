#!/usr/bin/env python3
"""
ra2 x tsh4 check for plant height (README section 3.6.4, PH Q3 and Q14).

RA2 represses tsh4 (Xie et al. 2023) and tsh4 acts against the ramosa genes (Xiao et al. 2022), so we ask
whether the teosinte segments at the two genes interact for PH. Simple check: ordinary least squares on
the combined PH BLUEs, PH ~ taxon + ra2 + tsh4 + ra2:tsh4, with the donor allele share of each pool at
each gene (B73 = 0, SEG = 0.5, DONOR = 1) from the k-mer HMM calls. No kinship or spatial term.

Usage:
  python final_script/ra2_tsh4_interaction.py --tracts tracts_all.tsv.gz --pool_summary pool_summary.tsv \
      --pheno BZea_spats_BLUEs_combined.csv [--trait PH]
Needs: pandas, statsmodels.
"""
import argparse
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

GENES = {"ra2": ("chr3", 12830057), "tsh4": ("chr7", 138540580)}  # Zm00001eb123060, Zm00001eb316740 (v5 start)

ap = argparse.ArgumentParser()
ap.add_argument("--tracts", required=True)
ap.add_argument("--pool_summary", required=True)
ap.add_argument("--pheno", required=True)
ap.add_argument("--trait", default="PH")
a = ap.parse_args()

t = pd.read_csv(a.tracts, sep="\t")
ps = pd.read_csv(a.pool_summary, sep="\t")
ok = set(ps.loc[ps.FLAG == "ok", "SAMPLE"])
d = pd.read_csv(a.pheno)[["Sequencing_ID", a.trait]].rename(columns={"Sequencing_ID": "SAMPLE", a.trait: "y"})
for g, (ch, pos) in GENES.items():
    x = t[(t.CHROM == ch) & (t.START <= pos) & (t.END > pos)]
    d[g] = d.SAMPLE.map(x.set_index("SAMPLE").STATE.map({"B73": 0, "SEG": 0.5, "DONOR": 1}))
d["taxon"] = d.SAMPLE.map(t.drop_duplicates("SAMPLE").set_index("SAMPLE").GROUP)
d = d[d.SAMPLE.isin(ok)].dropna()

m = smf.ols("y ~ C(taxon) + ra2 * tsh4", d).fit()
print(f"{a.trait}: {len(d)} pools")
print(m.summary().tables[1])
print("per copy (coefficient / 2):", {k: round(m.params[k] / 2, 2) for k in ["ra2", "tsh4", "ra2:tsh4"]})
d["group"] = np.select([(d.ra2 == 0) & (d.tsh4 == 0), (d.ra2 > 0) & (d.tsh4 == 0), (d.ra2 == 0) & (d.tsh4 > 0)],
                       ["neither", "ra2 only", "tsh4 only"], "both")
print(d.groupby("group").y.agg(["mean", "count"]).round(1))
