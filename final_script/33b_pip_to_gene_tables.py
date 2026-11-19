#!/usr/bin/env python3
"""
33b — Add the SuSiE fine-mapping result (step 33) to the per-trait gene tables (qtl_genes/<TRAIT>_qtl_genes.tsv).

New columns, after gwas_mlog10p:
  window_pip        posterior probability (SuSiE PIP) of the 100-kb window holding the gene (highest if it spans two)
  in_credible_set   yes / no: the window is in the QTL's 95% credible set
  in_top_windows    yes / no: the window is among the fewest windows that hold half of the QTL's probability
Run again after the gene tables are rebuilt; existing columns of the same name are replaced.

Usage: python final_script/33b_pip_to_gene_tables.py --susie data/introgression_spread/finemap_susie --genes qtl_genes
Needs: pandas.
"""
import argparse
import glob
import os
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--susie", required=True, help="step 33 output folder (susie_<TRAIT>_genes.tsv)")
ap.add_argument("--genes", required=True, help="folder with <TRAIT>_qtl_genes.tsv")
a = ap.parse_args()
NEW = ["window_pip", "in_credible_set", "in_top_windows"]

for f in sorted(glob.glob(os.path.join(a.genes, "*_qtl_genes.tsv"))):
    t = os.path.basename(f).replace("_qtl_genes.tsv", "")
    s = os.path.join(a.susie, f"susie_{t}_genes.tsv")
    if not os.path.exists(s):
        print(f"{t}: no {s}, skipped"); continue
    g = pd.read_csv(f, sep="\t").drop(columns=NEW, errors="ignore")
    p = pd.read_csv(s, sep="\t")[["QTL_window", "gene"] + NEW]
    for c in NEW[1:]:
        p[c] = p[c].map({True: "yes", False: "no"})
    g = g.merge(p, on=["QTL_window", "gene"], how="left")
    cols = [c for c in g.columns if c not in NEW]
    i = cols.index("gwas_mlog10p") + 1
    g = g[cols[:i] + NEW + cols[i:]]
    g.to_csv(f, sep="\t", index=False)
    print(f"{t}: {len(g)} genes, {g.window_pip.notna().sum()} with a PIP, {(g.in_top_windows == 'yes').sum()} in top windows")
