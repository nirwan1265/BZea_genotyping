#!/usr/bin/env python3
"""
NBR Q8 (chr10 98.2 Mb) and flowering time (README section 3.6.4).

NBR Q8 sits 1-2 Mb from ZmCCT10, so we ask whether its effect on nodes with brace roots runs through flowering
time. Simple check: ordinary least squares on the combined BLUEs, NBR ~ taxon + g [+ DTS (+ DTA + PH)], where g is the
donor allele share of each pool at Zm00001eb418930 (B73 = 0, SEG = 0.5, DONOR = 1) from the k-mer HMM calls.
No kinship or spatial term.

Usage:
  python final_script/nbr_q8_flowering_check.py --tracts tracts_all.tsv.gz --pool_summary pool_summary.tsv \
      --pheno BZea_spats_BLUEs_combined.csv
Needs: pandas, statsmodels.
"""
import argparse
import pandas as pd
import statsmodels.formula.api as smf

GENE = ("chr10", 97405806)    # Zm00001eb418930 (v5 start)
CCT10 = ("chr10", 96180000)   # ZmCCT10, Zm00001eb418700

ap = argparse.ArgumentParser()
ap.add_argument("--tracts", required=True)
ap.add_argument("--pool_summary", required=True)
ap.add_argument("--pheno", required=True)
a = ap.parse_args()

t = pd.read_csv(a.tracts, sep="\t")
ps = pd.read_csv(a.pool_summary, sep="\t")
ok = set(ps.loc[ps.FLAG == "ok", "SAMPLE"])
d = pd.read_csv(a.pheno).rename(columns={"Sequencing_ID": "SAMPLE"})


def dose(ch, pos):
    x = t[(t.CHROM == ch) & (t.START <= pos) & (t.END > pos)]
    return x.set_index("SAMPLE").STATE.map({"B73": 0, "SEG": 0.5, "DONOR": 1})


d["g"] = d.SAMPLE.map(dose(*GENE))
d["cct10"] = d.SAMPLE.map(dose(*CCT10))
d["taxon"] = d.SAMPLE.map(t.drop_duplicates("SAMPLE").set_index("SAMPLE").GROUP)
d = d[d.SAMPLE.isin(ok)].dropna(subset=["NBR", "DTS", "DTA", "PH", "g", "taxon"])
print(f"{len(d)} pools; {(d.g > 0).sum()} carry a donor segment at the gene; "
      f"{((d.g > 0) != (d.cct10 > 0)).sum()} differ between the gene and ZmCCT10")
print("correlation of NBR with:", d[["DTS", "DTA", "PH"]].corrwith(d.NBR).round(2).to_dict())
for f in ["NBR ~ C(taxon) + g", "NBR ~ C(taxon) + g + DTS", "NBR ~ C(taxon) + g + DTS + DTA + PH"]:
    m = smf.ols(f, d).fit()
    print(f"{f:40s} effect per copy {m.params['g'] / 2:+.3f}  p = {m.pvalues['g']:.2g}")
