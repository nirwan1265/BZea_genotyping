#!/usr/bin/env python
"""
tau_run / 3 — laptop: the targeted-sequencing check (step 34d) on the current calls and on the tau-floor calls of the
same pools, side by side per taxon (three bad samples of README 2.9.1 left out, tracts >= 1 Mb, clearly called).
Inputs: the tau_2 outputs copied to --run_dir (tracts_old.tsv.gz, tracts_new.tsv.gz, pools.tsv), and the current
model's targeted-sequencing sites (data/introgression_fix1/targseq/sites_used.tsv.gz).
Outputs in --run_dir: targseq_old/ and targseq_new/ (34d outputs), targseq_by_taxon.tsv (the side-by-side table).
Usage (from the repo root):
  python new_model/tau_run/tau_3_targseq.py --run_dir data/tau_run/zv_floor
Needs: numpy, pandas, scipy.
"""
import argparse
import os
import subprocess
import sys
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--run_dir", required=True)
ap.add_argument("--sites", default="data/introgression_fix1/targseq/sites_used.tsv.gz")
a = ap.parse_args()
BAD = ["Zv0480_P3_P2_P1211", "Zx0090_P2_P5_P1111", "Zx0160_P2_P4_P1111"]
EXP = {"B73": "AA", "SEG": "AB", "DONOR": "BB"}
s34d = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "final_script", "34d_targseq_tracts.py")

out = []
for name in ("old", "new"):
    od = f"{a.run_dir}/targseq_{name}"
    subprocess.run([sys.executable, s34d, "--sites", a.sites, "--tracts", f"{a.run_dir}/tracts_{name}.tsv.gz", "--out", od],
                   check=True, stdout=subprocess.DEVNULL)
    T = pd.read_csv(f"{od}/tracts_targseq.tsv", sep="\t")
    T = T[(T.Mb >= 1) & ~T["sample"].isin(BAD) & T.plant_genotype.isin(["AA", "AB", "BB"])]
    T["TAXON"] = T["sample"].str[:2]
    for tax, g in list(T.groupby("TAXON")) + [("all", T)]:
        r = {"run": name, "TAXON": tax, "samples": g["sample"].nunique()}
        for call, gt in EXP.items():
            x = g[g.our_call == call]
            r[f"{call}_clear"] = len(x)
            r[f"{call}_agree_%"] = round(100 * (x.plant_genotype == gt).mean(), 1) if len(x) else float("nan")
        r["SEG_look_BB"] = int(((g.our_call == "SEG") & (g.plant_genotype == "BB")).sum())
        r["DONOR_look_AB"] = int(((g.our_call == "DONOR") & (g.plant_genotype == "AB")).sum())
        out.append(r)
R = pd.DataFrame(out)
R = R.sort_values(["TAXON", "run"], key=lambda c: c.map({"old": 0, "new": 1}) if c.name == "run" else c)
R.to_csv(f"{a.run_dir}/targseq_by_taxon.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
print("targeted sequencing, current calls (old) vs tau floor (new), same pools:")
print(R.to_string(index=False))
P = pd.read_csv(f"{a.run_dir}/pools.tsv", sep="\t")
print("\nwindows that change state, per taxon (mean % of windows):")
print(P.groupby("TAXON")[["windows_changed_%", "SEG_old_%", "SEG_new_%", "DONOR_old_%", "DONOR_new_%"]].mean().round(2).to_string())
print(f"\nwrote {a.run_dir}/targseq_by_taxon.tsv")
