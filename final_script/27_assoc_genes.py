#!/usr/bin/env python3
"""
27 — Peaks from the step-26 scans and the B73 v5 genes inside them.

Peaks: windows above the threshold, joined when they are less than --gap (5 Mb) apart on a chromosome;
a peak inside a stronger peak's interval is dropped, so each region is reported once.
Interval: around each peak's top window, the windows within --drop (-log10 p units) of its top
value (a support interval, like a 1.5-LOD drop). Genes: every 'gene' feature of the B73 GFF3
overlapping the interval, with its ID, Name and description when the GFF has them.
--find REGEX also lists every gene (genome-wide) whose ID/Name/description matches, with the
scan's best -log10 p within --near bp of it (e.g. --find 'CCT' to locate ZmCCT10).

Outputs in <out_dir>: peaks_<TAG>.tsv, peak_genes_<TAG>.tsv, find_<TAG>.tsv (with --find)
"""
import argparse
import collections
import csv
import glob
import os
import re
from urllib.parse import unquote


def read_gff_genes(gff_glob):
    genes = []
    for f in sorted(glob.glob(gff_glob)):
        with open(f) as fh:
            for line in fh:
                if line.startswith("#"): continue
                p = line.rstrip("\n").split("\t")
                if len(p) < 9 or p[2] != "gene": continue
                a = dict(kv.split("=", 1) for kv in p[8].split(";") if "=" in kv)
                c = p[0] if p[0].startswith("chr") else "chr" + p[0]
                genes.append(dict(chr=c, start=int(p[3]), end=int(p[4]),
                                  id=a.get("gene_id", a.get("ID", "").replace("gene:", "")),
                                  name=unquote(a.get("Name", "")), desc=unquote(a.get("description", ""))))
    return genes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scan", required=True, help="step-26 scan_<TRAIT>[_<GROUP>].tsv")
    ap.add_argument("--threshold", type=float, required=True, help="-log10 p threshold (step-26 summary)")
    ap.add_argument("--gff", required=True, help="glob of B73 v5 GFF3 file(s), quoted")
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--tag", default=None)
    ap.add_argument("--gap", type=int, default=5_000_000)
    ap.add_argument("--drop", type=float, default=1.5)
    ap.add_argument("--find", default=None)
    ap.add_argument("--near", type=int, default=5_000_000)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    tag = a.tag or os.path.basename(a.scan).replace("scan_", "").replace(".tsv", "")

    by = collections.defaultdict(list)
    with open(a.scan) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            by[r["chr"]].append((int(r["start"]), int(r["end"]), float(r["mlog10p_add"]), float(r["beta_add"]), r["Marker"]))
    for c in by: by[c].sort()
    genes = read_gff_genes(a.gff)
    gby = collections.defaultdict(list)
    for g in genes: gby[g["chr"]].append(g)
    print(f"genes in GFF: {len(genes)}")

    peaks = []
    for c, ws in by.items():
        sig = [i for i, w in enumerate(ws) if w[2] > a.threshold]
        groups = []
        for i in sig:
            if groups and ws[i][0] - ws[groups[-1][-1]][1] <= a.gap: groups[-1].append(i)
            else: groups.append([i])
        for gi in groups:
            top = max(gi, key=lambda i: ws[i][2]); m = ws[top][2]
            lo = top
            while lo > 0 and ws[lo - 1][2] >= m - a.drop: lo -= 1
            hi = top
            while hi < len(ws) - 1 and ws[hi + 1][2] >= m - a.drop: hi += 1
            peaks.append(dict(chr=c, sig_start=ws[gi[0]][0], sig_end=ws[gi[-1]][1], top=ws[top][4], top_mlog10p=m,
                              beta=ws[top][3], int_start=ws[lo][0], int_end=ws[hi][1]))
    peaks.sort(key=lambda p: -p["top_mlog10p"])
    kept = []                                   # one peak per region: drop peaks inside a stronger peak's interval
    for p in peaks:
        if not any(q["chr"] == p["chr"] and p["int_start"] <= q["int_end"] and p["int_end"] >= q["int_start"] for q in kept):
            kept.append(p)
    peaks = kept
    with open(f"{a.out_dir}/peaks_{tag}.tsv", "w") as f, open(f"{a.out_dir}/peak_genes_{tag}.tsv", "w") as g:
        f.write("PEAK\tCHR\tSIG_START\tSIG_END\tTOP_WINDOW\tTOP_MLOG10P\tBETA_DAYS_PER_COPY\tINTERVAL_START\tINTERVAL_END\tN_GENES\n")
        g.write("PEAK\tCHR\tSTART\tEND\tGENE\tNAME\tDESCRIPTION\n")
        for k, p in enumerate(peaks, 1):
            gl = [x for x in gby[p["chr"]] if x["end"] >= p["int_start"] and x["start"] <= p["int_end"]]
            f.write(f"{k}\t{p['chr']}\t{p['sig_start']}\t{p['sig_end']}\t{p['top']}\t{p['top_mlog10p']:.2f}\t{p['beta']:.3f}\t"
                    f"{p['int_start']}\t{p['int_end']}\t{len(gl)}\n")
            for x in gl: g.write(f"{k}\t{x['chr']}\t{x['start']}\t{x['end']}\t{x['id']}\t{x['name']}\t{x['desc']}\n")
            print(f"peak {k}: {p['chr']}:{p['int_start']/1e6:.1f}-{p['int_end']/1e6:.1f} Mb, top {p['top']} "
                  f"-log10p {p['top_mlog10p']:.2f}, {p['beta']:+.2f} days/copy, {len(gl)} genes")
    if a.find:
        rx = re.compile(a.find, re.I)
        with open(f"{a.out_dir}/find_{tag}.tsv", "w") as f:
            f.write("GENE\tNAME\tCHR\tSTART\tEND\tBEST_MLOG10P_WITHIN_NEAR\tDESCRIPTION\n")
            for x in genes:
                if not (rx.search(x["id"]) or rx.search(x["name"]) or rx.search(x["desc"])): continue
                near = [w[2] for w in by.get(x["chr"], []) if w[1] >= x["start"] - a.near and w[0] <= x["end"] + a.near]
                best = max(near) if near else float("nan")
                f.write(f"{x['id']}\t{x['name']}\t{x['chr']}\t{x['start']}\t{x['end']}\t{best:.2f}\t{x['desc']}\n")
                print(f"find: {x['id']} {x['name']} {x['chr']}:{x['start']}-{x['end']} best -log10p within "
                      f"{a.near/1e6:.0f} Mb: {best:.2f}  {x['desc'][:80]}")
    print(f"WROTE: {a.out_dir}/peaks_{tag}.tsv peak_genes_{tag}.tsv" + (f" find_{tag}.tsv" if a.find else ""))


if __name__ == "__main__":
    main()
