#!/usr/bin/env python3
"""
33a — Per-trait table of every gene inside each joint-model QTL support interval (qtl_genes/<TRAIT>_qtl_genes.tsv),
the input of step 33 (genes table) and 33b (adds window_pip, in_credible_set, in_top_windows).

QTL: step 29's qtl2_<TRAIT>_qtl.tsv, numbered Q1, Q2, ... in file order; holds_joint_model = joint -log10 p above
the Bonferroni cut over the trait's QTL (-log10(0.05 / number of QTL)), as in step 33 and the README.
Genes: step 29's qtl2_<TRAIT>_genes.tsv (every B73 v5 gene in each support interval). gwas_mlog10p = highest
-log10 p of the step-26 scan (scan_<TRAIT>.tsv, additive) over the windows the gene overlaps. Name, family and GO
biological process: step 30's genes_GO_all.tsv (all traits; taxon and subtaxon runs), GO ids dropped, cut at 200 characters.
Usage: python final_script/33a_qtl_gene_tables.py --qtl_dir <joint_qtl/noperm> --scan_dir <assoc_gridlmm_spats> \
          --go_all <GO/genes_GO_all.tsv> <GO_subtaxon/genes_GO_all.tsv> --out qtl_genes [--traits "DTS DTA"]
Needs: pandas.
"""
import argparse
import math
import os
import re
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--qtl_dir", required=True, help="step 29 output: qtl2_<TRAIT>_qtl.tsv, qtl2_<TRAIT>_genes.tsv")
ap.add_argument("--scan_dir", required=True, help="step 26 output: scan_<TRAIT>.tsv")
ap.add_argument("--go_all", required=True, nargs="+", help="step 30 output(s): genes_GO_all.tsv (taxon and subtaxon runs)")
ap.add_argument("--out", required=True)
ap.add_argument("--traits", default="DTS DTA ASI NBR SPAD2 PH EH EN Prolif SL BL BW")
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
GO = pd.concat([pd.read_csv(f, sep="\t") for f in a.go_all if os.path.exists(f)]).drop_duplicates("GENE").set_index("GENE")
clean = lambda v: "" if pd.isna(v) else str(v)
for T in a.traits.split():
    qf, gf, sf = (f"{a.qtl_dir}/qtl2_{T}_qtl.tsv", f"{a.qtl_dir}/qtl2_{T}_genes.tsv", f"{a.scan_dir}/scan_{T}.tsv")
    miss = [f for f in (qf, gf, sf) if not os.path.exists(f)]
    if miss: print(f"{T}: missing {miss[0]}, skipped"); continue
    Q = pd.read_csv(qf, sep="\t"); G = pd.read_csv(gf, sep="\t"); S = pd.read_csv(sf, sep="\t")
    cut = -math.log10(0.05 / len(Q))
    rows = []
    for i, q in enumerate(Q.itertuples(), 1):
        for r in G[G.qtl == q.qtl].sort_values("gene_start").itertuples():
            w = S[(S.chr == r.chr) & (S.start <= r.gene_end) & (S.end >= r.gene_start)].mlog10p_add
            an = GO.loc[r.gene] if r.gene in GO.index else None
            fam = re.sub(r"\s*\(PTHR[^)]*\)", "", clean(an.Family_Subfamily)) if an is not None else ""
            gobp = re.sub(r"\(GO:\d+\)", "", clean(an.GO_BP))[:200] if an is not None else ""
            rows.append({"QTL": f"Q{i}", "QTL_window": q.qtl, "holds_joint_model": "yes" if q.mlog10p_joint > cut else "no",
                         "joint_mlog10p": q.mlog10p_joint, "effect_per_copy_joint": q.effect_per_copy_joint,
                         "gene": r.gene, "chr": r.chr, "start": r.gene_start, "end": r.gene_end,
                         "gwas_mlog10p": round(float(w.max()), 2) if len(w) else "",
                         "gene_name": clean(an.GeneName) if an is not None else "", "family": fam, "GO_biological_process": gobp})
    out = pd.DataFrame(rows)
    out.to_csv(f"{a.out}/{T}_qtl_genes.tsv", sep="\t", index=False)
    print(f"{T}: {len(Q)} QTL ({sum(Q.mlog10p_joint > cut)} hold), {len(out)} gene rows -> {a.out}/{T}_qtl_genes.tsv")
