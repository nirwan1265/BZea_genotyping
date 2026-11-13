#!/usr/bin/env python3
"""
12_signal_check.py - is the teosinte-panel site set strong enough on its own, per taxon?

For every pool, step 11 summed B73 and ALT reads over PASS sites that are
  (a) in the teosinte panel, and (b) not in the panel.
Here, for each site class:
  background  = ALT fraction of the clean B73 controls (errors, leftover artifacts)
  excess      = pool ALT fraction - background
  d           = excess / 0.125   (a BC2-derived line is ~12.5% donor on average, so this
                is roughly the ALT fraction a fully-donor region adds above background)
Then, per 1 Mb window, for a typical pool of each taxon:
  reads/Mb          reads that pool has on sites of this class
  ALT/Mb in B73     reads/Mb x background
  ALT/Mb in donor   reads/Mb x (background + d)
  z_fixed           (ALT/Mb donor - ALT/Mb B73) / sqrt(ALT/Mb B73 + 1)   ~ how clearly a
                    fixed-donor window stands out; segregating windows give about half.
Rule of thumb: z_fixed >= 10 on panel sites -> panel sites alone are enough for that taxon.

usage: 12_signal_check.py --d11 DIR --genome-mb MB --out DIR
"""
import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--d11", required=True)
    ap.add_argument("--genome-mb", type=float, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--donor-share", type=float, default=0.125)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    files = sorted(glob.glob(os.path.join(a.d11, "*.pool_class.tsv")))
    if not files:
        sys.exit(f"ERROR: no *.pool_class.tsv in {a.d11}")
    pc = pd.concat([pd.read_csv(f, sep="\t", dtype={"GROUP": str, "ACCESSION": str}) for f in files])
    num = ["PANEL_REF", "PANEL_ALT", "NONPANEL_REF", "NONPANEL_ALT"]
    p = pc.groupby(["SAMPLE", "ROLE", "GROUP", "ACCESSION"], dropna=False, as_index=False)[num].sum()

    summ = pd.concat([pd.read_csv(f, sep="\t") for f in sorted(glob.glob(os.path.join(a.d11, "*.summary.tsv")))])
    summ.to_csv(os.path.join(a.out, "site_filter_summary.tsv"), sep="\t", index=False)

    bg = {}
    for cls in ("PANEL", "NONPANEL"):
        c = p[p.ROLE == "control"]
        r, x = c[f"{cls}_REF"].sum(), c[f"{cls}_ALT"].sum()
        bg[cls] = x / max(r + x, 1)
        p[f"{cls}_READS"] = p[f"{cls}_REF"] + p[f"{cls}_ALT"]
        p[f"{cls}_ALT_FRAC"] = p[f"{cls}_ALT"] / p[f"{cls}_READS"].clip(lower=1)
        p[f"{cls}_EXCESS"] = p[f"{cls}_ALT_FRAC"] - bg[cls]
    p.to_csv(os.path.join(a.out, "pool_signal.tsv"), sep="\t", index=False)
    ctrl = p[p.ROLE == "control"]

    rows = []
    for g, d in p[p.ROLE == "donor"].groupby("GROUP"):
        row = {"TAXON": g, "N_POOLS": len(d)}
        for cls in ("PANEL", "NONPANEL"):
            rpm = d[f"{cls}_READS"].median() / a.genome_mb
            exc = d[f"{cls}_EXCESS"].median()
            dd = max(exc, 0) / a.donor_share
            alt_b73 = rpm * bg[cls]
            alt_don = rpm * (bg[cls] + dd)
            row.update({f"{cls}_READS_PER_MB": round(rpm, 1),
                        f"{cls}_MEDIAN_EXCESS": round(exc, 5),
                        f"{cls}_ALT_PER_MB_B73": round(alt_b73, 2),
                        f"{cls}_ALT_PER_MB_DONOR": round(alt_don, 1),
                        f"{cls}_Z_FIXED": round((alt_don - alt_b73) / np.sqrt(alt_b73 + 1), 1)})
        row["N_POOLS_NO_SIGNAL_PANEL"] = int((d.PANEL_EXCESS <= 2 * ctrl.PANEL_ALT_FRAC.std(ddof=0)).sum()) if len(ctrl) > 1 else -1
        row["PANEL_ENOUGH"] = "yes" if row["PANEL_Z_FIXED"] >= 10 else "NO - add non-panel sites"
        rows.append(row)
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(a.out, "taxon_signal.tsv"), sep="\t", index=False)

    tot = summ[["sites", "PASS", "PASS_in_panel", "PASS_not_in_panel", "ctrl_alt", "dp_high", "dp_low", "oth", "panel_alleles"]].sum()
    print("Sites (all chromosomes):")
    print(tot.to_string())
    print(f"\nB73 background ALT fraction: panel sites {bg['PANEL']:.5f} | non-panel sites {bg['NONPANEL']:.5f}")
    print(f"pools: {p.ROLE.value_counts().to_dict()}\n")
    pd.set_option("display.width", 250)
    print(t[["TAXON", "N_POOLS", "PANEL_READS_PER_MB", "PANEL_ALT_PER_MB_B73", "PANEL_ALT_PER_MB_DONOR", "PANEL_Z_FIXED",
             "NONPANEL_ALT_PER_MB_B73", "NONPANEL_ALT_PER_MB_DONOR", "NONPANEL_Z_FIXED", "PANEL_ENOUGH"]].to_string(index=False))
    print("\nWrote:", *[os.path.join(a.out, x) for x in ("taxon_signal.tsv", "pool_signal.tsv", "site_filter_summary.tsv")], sep="\n  ")


if __name__ == "__main__":
    main()
