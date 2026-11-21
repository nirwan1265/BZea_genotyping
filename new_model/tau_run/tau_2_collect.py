#!/usr/bin/env python
"""
tau_run / 2 — current model vs the tau-floor rerun, on the pools of tau_samples.tsv (run on the cluster, after tau_1).

Both runs are read at the current model's rho (the best total of ${FX}/with_levels/introgression/rho_loglik.tsv,
as in step 21), states_R1 (the calls of step 21).
Outputs (--out):
  pools.tsv           per pool: taxon, scale a (old/new), SEG and DONOR share (old/new), % windows that change,
                      windows per change (B73>SEG, SEG>DONOR, ...), tracts >= 1 Mb (old/new)
  changed_regions.tsv runs of windows whose state changed: pool, chr, start, end, Mb, old state, new state
  tracts_old.tsv.gz / tracts_new.tsv.gz   the two runs' tracts of these pools, in the format of tracts_all.tsv.gz
                      (for the targeted-sequencing check on the laptop, tau_3_targseq.sh)
Usage (cluster, from new_model/tau_run/, python env active):
  source ../../final_script/config.sh
  python tau_2_collect.py --samples tau_samples.tsv --old_hmm /share/maize/ntanduk/kmer/new_model/fix1/with_levels/hmm \
    --new_hmm /share/maize/ntanduk/kmer/new_model/tau_run/zv_floor/hmm \
    --rho_tsv /share/maize/ntanduk/kmer/new_model/fix1/with_levels/introgression/rho_loglik.tsv \
    --ref_dir $D19 --out /share/maize/ntanduk/kmer/new_model/tau_run/zv_floor/compare
Needs: numpy, pandas.
"""
import argparse
import gzip
import importlib.util
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("c21", os.path.join(HERE, "..", "..", "final_script", "21_kmer_collect.py"))
C21 = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(C21)
ST = ("B73", "SEG", "DONOR")

ap = argparse.ArgumentParser()
ap.add_argument("--samples", required=True)
ap.add_argument("--old_hmm", required=True)
ap.add_argument("--new_hmm", required=True)
ap.add_argument("--rho_tsv", required=True)
ap.add_argument("--ref_dir", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--window", type=int, default=100000)
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)

chrom, start = C21.read_windows(f"{a.ref_dir}/B73.tsv.gz")
R = pd.read_csv(a.rho_tsv, sep="\t")
rho = float(R.RHO[R.ALL.idxmax()])
print(f"current model's rho: {rho:g}")
S = pd.read_csv(a.samples, sep="\t")


def load(f):
    with np.load(f, allow_pickle=True) as z:
        return {k: z[k] for k in z.files}


rows, chg, tr = [], [], {"old": [], "new": []}
for s in S.itertuples():
    fo, fn = f"{a.old_hmm}/{s.SAMPLE}.npz", f"{a.new_hmm}/{s.SAMPLE}.npz"
    if not (os.path.exists(fo) and os.path.exists(fn)):
        print(f"{s.SAMPLE}: missing {'old' if not os.path.exists(fo) else 'new'} run, skipped"); continue
    O, N = load(fo), load(fn)
    jo = int(np.where(np.isclose(O["rho_grid"], rho))[0][0]); jn = int(np.where(np.isclose(N["rho_grid"], rho))[0][0])
    so, sn = O["states_R1"][jo], N["states_R1"][jn]
    m = (so >= 0) & (sn >= 0)
    r = {"SAMPLE": s.SAMPLE, "TAXON": s.TAXON, "TARGSEQ_SAMPLE": s.TARGSEQ_SAMPLE,
         "scale_old": float(O["par_best"][jo]), "scale_new": float(N["par_best"][jn])}
    for k, st in ((1, "SEG"), (2, "DONOR")):
        r[f"{st}_old_%"] = round(100 * np.mean(so[m] == k), 2); r[f"{st}_new_%"] = round(100 * np.mean(sn[m] == k), 2)
    r["windows_changed_%"] = round(100 * np.mean(so[m] != sn[m]), 2)
    for i in range(3):
        for k in range(3):
            if i != k: r[f"{ST[i]}>{ST[k]}"] = int(np.sum(m & (so == i) & (sn == k)))
    for name, st in (("old", so), ("new", sn)):
        bl = C21.runs(st, chrom, start, a.window)
        r[f"tracts_1Mb_{name}"] = sum(1 for b in bl if b[3] > 0 and b[2] - b[1] >= 1e6)
        tr[name] += [(s.SAMPLE, s.TAXON, c, b0, b1, ST[k]) for c, b0, b1, k in bl]
    rows.append(r)
    d = (so != sn) & m                                              # runs of changed windows
    for c in dict.fromkeys(chrom):
        idx = np.where((chrom == c) & d)[0]
        if not len(idx): continue
        brk = np.where(np.diff(idx) > 1)[0]
        for g in np.split(idx, brk + 1):
            chg.append({"SAMPLE": s.SAMPLE, "TAXON": s.TAXON, "CHROM": c, "START": int(start[g[0]]),
                        "END": int(start[g[-1]]) + a.window, "Mb": round(len(g) * a.window / 1e6, 1),
                        "OLD": ST[int(np.bincount(so[g]).argmax())], "NEW": ST[int(np.bincount(sn[g]).argmax())]})

P = pd.DataFrame(rows)
P.to_csv(f"{a.out}/pools.tsv", sep="\t", index=False)
pd.DataFrame(chg).to_csv(f"{a.out}/changed_regions.tsv", sep="\t", index=False)
for name in ("old", "new"):
    with gzip.open(f"{a.out}/tracts_{name}.tsv.gz", "wt") as f:
        f.write("SAMPLE\tGROUP\tCHROM\tSTART\tEND\tSTATE\n")
        for t in tr[name]: f.write("\t".join(map(str, t)) + "\n")
pd.set_option("display.width", 250)
print(P.to_string(index=False))
print("\nper taxon: mean % windows changed, SEG and DONOR share old -> new")
print(P.groupby("TAXON")[["windows_changed_%", "SEG_old_%", "SEG_new_%", "DONOR_old_%", "DONOR_new_%"]].mean().round(2).to_string())
print(f"\nwrote {a.out}/pools.tsv, changed_regions.tsv, tracts_old.tsv.gz, tracts_new.tsv.gz")
