#!/usr/bin/env python3
"""
35 — README tables of the k-mer ancestry GWAS + joint QTL model + fine-mapping, built from the result files
(so the numbers in the README come from the data, not by hand).

Writes to --out (one markdown file per trait, plus overview.md):
  overview.md          one row per trait: plots, threshold, QTL, QTL that hold, strongest scan window
  <TRAIT>.md           the per-QTL table (position, interval, scan and joint -log10 p, effect, taxa differ,
                       per-taxon effects, genes in the interval, top windows, best-window PIP), the Bonferroni cut,
                       and, for every gene ID given with --genes, which QTL interval holds it and its PIP
Usage: python final_script/35_trait_tables.py --base data/introgression_fix1 --out <dir> [--genes Zm00001eb418700,...]
Needs: pandas.
"""
import argparse
import math
import os
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--base", required=True, help="data/introgression_<run>: assoc_gridlmm_spats, joint_qtl/noperm, finemap_susie, qtl_genes")
ap.add_argument("--out", required=True)
ap.add_argument("--traits", default="DTS DTA ASI NBR SPAD2 PH EH EN Prolif SL BL BW SPAD Kinki LAE StPu StPi")
ap.add_argument("--sub", default="", help="'_subtaxon' for the races run")
ap.add_argument("--genes", default="", help="comma-separated gene IDs to locate (named candidates)")
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
UNIT = {"DTS": " d", "DTA": " d", "ASI": " d", "PH": " cm", "EH": " cm"}
B = a.base
link = lambda T: f"{B}/qtl_genes/{T}_qtl_genes.tsv"
ov = []
for T in a.traits.split():
    sf, qf = f"{B}/assoc_gridlmm_spats{a.sub}/summary_{T}.tsv", f"{B}/joint_qtl/noperm{a.sub}/qtl2_{T}_qtl.tsv"
    if not (os.path.exists(sf) and os.path.exists(qf)):
        continue
    S = pd.read_csv(sf, sep="\t").iloc[0]
    Q = pd.read_csv(qf, sep="\t")
    cut = -math.log10(0.05 / len(Q)) if len(Q) else float("nan")
    Q["holds"] = Q.mlog10p_joint > cut
    u = UNIT.get(T, "")
    c, p = S.top_window.split(":")
    ov.append(f"| {T} | {S.plots} | {S.threshold:.2f} | {len(Q)} | {int(Q.holds.sum())} | {c} {int(p) / 1e6:.1f} Mb "
              f"(−log10 p {S.top_mlog10p:.1f}, {S.top_beta:+.2f}{u} per copy) |")
    BT = pd.read_csv(f"{B}/joint_qtl/noperm{a.sub}/qtl2_{T}_by_taxon.tsv", sep="\t")
    SU = f"{B}/finemap_susie/susie_{T}_qtl.tsv"
    SU = pd.read_csv(SU, sep="\t").set_index("qtl") if os.path.exists(SU) and not a.sub else None
    G = pd.read_csv(link(T), sep="\t") if os.path.exists(link(T)) and not a.sub else None
    L = [f"Bonferroni over the {len(Q)} QTL: −log10 p > {cut:.2f} (bold).", "",
         "| QTL | Position | Interval (Mb) | Scan −log10 p | Effect per copy (joint) | Joint −log10 p | Taxa differ −log10 p | "
         "Per-taxon effects (carrier plots) | Genes in interval | Top windows, Mb (genes) | Best window PIP |",
         "|---|---|---|---|---|---|---|---|---|---|---|"]
    for i, q in enumerate(Q.itertuples(), 1):
        bt = BT[(BT.qtl == q.qtl) & (BT.carriers >= 5) & BT.effect_per_copy.notna()]
        pt = ", ".join(f"{r.taxon} {r.effect_per_copy:+.2f} ({r.carriers})" for r in bt.itertuples())
        jp = f"**{q.mlog10p_joint:.1f}**" if q.holds else f"{q.mlog10p_joint:.1f}"
        td = q.mlog10p_taxa_differ
        tds = "" if pd.isna(td) else (f"**{td:.1f}**" if td > cut else f"{td:.1f}")
        ng = int((G.QTL == f"Q{i}").sum()) if G is not None else ""
        gl = f"[{ng} genes]({link(T)})" if G is not None else ""
        tw, bp = "", ""
        if SU is not None and q.qtl in SU.index:
            s = SU.loc[q.qtl]
            if pd.notna(s.top_windows_span):
                t0, t1 = (float(v) for v in str(s.top_windows_span).split("-"))
                tw = f"{t0:.1f}–{t1:.1f} ({int(s.genes_top_windows)})"
            else:
                tw = "–"
            bp = f"{s.best_window_pip:.2f}" if pd.notna(s.best_window_pip) else "–"
        tri = " (cond.)" if str(q.found_by).startswith("conditional") else ""
        L.append(f"| Q{i}{tri} | {q.chr} {q.pos_Mb:.1f} Mb | {q.interval_start_Mb:.1f}–{q.interval_end_Mb:.1f} | {q.gridlmm_mlog10p:.1f} | "
                 f"{q.effect_per_copy_joint:+.2f}{u} ± {q.se:.2f} | {jp} | {tds} | {pt} | {gl} | {tw} | {bp} |")
    if a.genes and G is not None:
        L += ["", "Named genes:"]
        for g in a.genes.split(","):
            x = G[G.gene == g]
            for r in x.itertuples():
                L.append(f"- {g} ({r.gene_name}): {r.QTL} {r.QTL_window}, holds {r.holds_joint_model}, window PIP {r.window_pip}, "
                         f"credible set {r.in_credible_set}, top windows {r.in_top_windows}")
    open(f"{a.out}/{T}.md", "w").write("\n".join(L) + "\n")
open(f"{a.out}/overview.md", "w").write("| Trait | Plots | 5% threshold (−log10 p) | QTL | Hold in joint model | Strongest QTL |\n"
                                        "|---|---|---|---|---|---|\n" + "\n".join(ov) + "\n")
print("wrote", a.out)
