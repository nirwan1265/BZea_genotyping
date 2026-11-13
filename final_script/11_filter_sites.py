#!/usr/bin/env python3
"""
11_filter_sites.py - site filters, teosinte-panel annotation and donor-group read
counts for one chromosome. Called by 11_filter_sites.sh (two sub-commands).

  pools   : build the pool table (role / taxon group / accession for every pool)
  filter  : stream the AD counts from bcftools query and write the outputs

Site filters (a site must pass all of them):
  ctrl_alt   ALT reads in the clean B73 controls > --max-ctrl-alt
             (a real B73 plant cannot carry ALT: repeat, paralog or reference error)
  dp_high    total depth > --max-dp-mult x chromosome median   (collapsed repeats)
  dp_low     total depth < --min-dp-mult x chromosome median   (hard-to-map region)
  oth        third-base reads / depth > --max-oth               (messy site)
  panel_alleles  site is in the teosinte panel but with different alleles

Outputs (prefix = --out):
  .sites.tsv.gz      one row per site: depth, control reads, filter result, teosinte
                     panel AC/AN, and B73/ALT read sums per donor taxon (included pools)
  .acc_counts.npz    PASS sites only: B73/ALT read sums per donor accession
                     (used later for "other pools of the same accession" evidence)
  .pool_class.tsv    per pool: B73/ALT reads on PASS sites in the panel / not in the panel.
                     For the B73 controls, each control is counted on sites that pass
                     the filters using the OTHER controls (leave-one-out), so their ALT
                     fraction is a fair background estimate.
  .summary.tsv       how many sites each filter removed
Needs numpy and pandas.
"""
import argparse
import gzip
import sys

import numpy as np
import pandas as pd

BASES = {"A": 0, "C": 1, "G": 2, "T": 3}


def read_list(path):
    with open(path) as fh:
        return {l.strip() for l in fh if l.strip()}


# --------------------------------------------------------------------------- pools
def cmd_pools(a):
    samples = [l.strip() for l in open(a.samples) if l.strip()]
    meta = pd.read_csv(a.metadata, sep="\t", dtype=str)
    meta = meta.apply(lambda s: s.str.replace("\r", "", regex=False).str.strip())
    meta.columns = ["SAMPLE", "GENOTYPE"] + list(meta.columns[2:])
    geno = meta.groupby("SAMPLE")["GENOTYPE"].agg(lambda x: "|".join(sorted(set(x)))).to_dict()
    controls = read_list(a.controls)
    exclude = read_list(a.exclude)

    rows = []
    for s in samples:
        g = geno.get(s, "")
        if s in controls:
            role, group, acc = "control", "B73", "B73"
        elif s in exclude:
            role, group, acc = "excluded", "", ""
        elif len(g) > 3 and g[0] == "Z" and g[2] == "." and "|" not in g:
            role, group, acc = "donor", g.split(".")[0], g.split("_")[0]
        else:
            role, group, acc = "other", g, ""
        rows.append((s, role, group, acc, g))
    t = pd.DataFrame(rows, columns=["SAMPLE", "ROLE", "GROUP", "ACCESSION", "GENOTYPE"])
    t.to_csv(a.out, sep="\t", index=False)
    keep = t[t.ROLE.isin(["control", "donor"])]
    keep.SAMPLE.to_csv(a.keep, index=False, header=False)
    sys.stderr.write(t.groupby("ROLE").size().to_string() + "\n")
    sys.stderr.write("donor pools per group: " + keep[keep.ROLE == "donor"].groupby("GROUP").size().to_dict().__str__() + "\n")
    missing = controls - set(samples)
    if missing:
        sys.stderr.write(f"WARNING: controls not in samples.txt: {sorted(missing)}\n")


# --------------------------------------------------------------------------- filter
def load_panel(path):
    p = pd.read_csv(path, sep="\t", header=None, names=["POS", "REF", "ALT", "AC", "AN"],
                    dtype={"POS": np.int64, "REF": str, "ALT": str, "AC": str, "AN": str})
    p["AC"] = pd.to_numeric(p.AC, errors="coerce").fillna(-1).astype(np.int32)
    p["AN"] = pd.to_numeric(p.AN, errors="coerce").fillna(-1).astype(np.int32)
    p = p.drop_duplicates("POS", keep=False).sort_values("POS")
    return p.POS.to_numpy(), p.REF.to_numpy(), p.ALT.to_numpy(), p.AC.to_numpy(), p.AN.to_numpy()


def cmd_filter(a):
    pools = pd.read_csv(a.pools, sep="\t", dtype=str).fillna("")
    keep = [l.strip() for l in open(a.keep) if l.strip()]
    pools = pools.set_index("SAMPLE").loc[keep].reset_index()      # query column order
    P = len(pools)
    is_ctrl = (pools.ROLE == "control").to_numpy()
    groups = sorted(pools.loc[pools.ROLE == "donor", "GROUP"].unique())
    gcols = {g: np.where(((pools.ROLE == "donor") & (pools.GROUP == g)).to_numpy())[0] for g in groups}
    accs = sorted(pools.loc[pools.ROLE == "donor", "ACCESSION"].unique())
    acc_idx = [np.where(((pools.ROLE == "donor") & (pools.ACCESSION == x)).to_numpy())[0] for x in accs]
    acc_perm = np.concatenate(acc_idx)
    acc_starts = np.cumsum([0] + [len(i) for i in acc_idx[:-1]])
    sys.stderr.write(f"pools: {P} ({is_ctrl.sum()} controls) | groups: {groups} | accessions: {len(accs)}\n")

    dp_all = pd.read_csv(a.sitesum, sep="\t", usecols=["DP"]).DP.to_numpy()
    med = float(np.median(dp_all))
    hi, lo = a.max_dp_mult * med, a.min_dp_mult * med
    sys.stderr.write(f"median site depth {med:.0f}: keep {lo:.0f} <= DP <= {hi:.0f}\n")

    ppos, pref, palt, pac, pan = load_panel(a.panel)
    sys.stderr.write(f"teosinte panel SNPs: {len(ppos)}\n")

    dtypes = {0: np.int64, 1: str, 2: str, 3: np.int64, 4: np.int64}
    dtypes.update({i: np.int32 for i in range(5, 5 + 2 * P)})
    reader = pd.read_csv(sys.stdin, sep="\t", header=None, dtype=dtypes, chunksize=a.chunk,
                         na_values=["."], keep_default_na=False)

    cls_ref = np.zeros((2, P), np.int64)       # class 0 = PASS & in panel, 1 = PASS & not in panel
    cls_alt = np.zeros((2, P), np.int64)
    ctrl_idx = np.where(is_ctrl)[0]
    hold_ref = np.zeros((2, P), np.int64)      # controls: leave-one-out version of the above
    hold_alt = np.zeros((2, P), np.int64)
    acc_pos, acc_ref, acc_alt = [], [], []
    counts = dict(sites=0, ctrl_alt=0, dp_high=0, dp_low=0, oth=0, panel_alleles=0,
                  PASS=0, PASS_in_panel=0, PASS_not_in_panel=0, in_panel=0)

    out = gzip.open(a.out + ".sites.tsv.gz", "wt", compresslevel=4)
    header = ["CHROM", "POS", "B73", "ALT", "DP", "OTH", "CTRL_REF", "CTRL_ALT", "FILTER",
              "IN_PANEL", "PANEL_AC", "PANEL_AN"] + [f"{g}_{x}" for g in groups for x in ("REF", "ALT")]
    out.write("\t".join(header) + "\n")

    for ch in reader:
        if ch.shape[1] != 5 + 2 * P:
            sys.exit(f"ERROR: got {ch.shape[1]} columns, expected {5 + 2 * P}")
        pos = ch[0].to_numpy(); ref = ch[1].to_numpy(); alt = ch[2].to_numpy()
        dp = ch[3].to_numpy(); oth = ch[4].to_numpy()
        M = ch.iloc[:, 5:].to_numpy(dtype=np.int32, na_value=0)
        R, A = M[:, 0::2], M[:, 1::2]
        n = len(pos)

        cref = R[:, is_ctrl].sum(1); calt = A[:, is_ctrl].sum(1)
        i = np.searchsorted(ppos, pos); i[i >= len(ppos)] = len(ppos) - 1
        inp = ppos[i] == pos
        same = inp & (pref[i] == ref) & (palt[i] == alt)
        f_ctrl = calt > a.max_ctrl_alt
        f_hi = dp > hi; f_lo = dp < lo
        f_oth = (dp > 0) & (oth > a.max_oth * dp)
        f_pan = inp & ~same
        ok = ~(f_ctrl | f_hi | f_lo | f_oth | f_pan)

        filt = np.full(n, "", dtype=object)
        for f, name in ((f_ctrl, "ctrl_alt"), (f_hi, "dp_high"), (f_lo, "dp_low"), (f_oth, "oth"), (f_pan, "panel_alleles")):
            filt[f] = filt[f] + name + ";"
            counts[name] += int(f.sum())
        filt = np.where(ok, "PASS", np.char.rstrip(filt.astype(str), ";"))

        counts["sites"] += n; counts["PASS"] += int(ok.sum()); counts["in_panel"] += int(same.sum())
        for c, m in enumerate((ok & same, ok & ~same)):
            cls_ref[c] += R[m].sum(0); cls_alt[c] += A[m].sum(0)
        # controls: score each one on sites filtered by the OTHER controls only
        # (otherwise its ALT count on PASS sites is 0 by construction)
        ok_other = ~(f_hi | f_lo | f_oth | f_pan)
        for k in ctrl_idx:
            mk = ok_other & ((calt - A[:, k]) <= a.max_ctrl_alt)
            for c, m in enumerate((mk & same, mk & ~same)):
                hold_ref[c, k] += R[m, k].sum(); hold_alt[c, k] += A[m, k].sum()
        counts["PASS_in_panel"] += int((ok & same).sum()); counts["PASS_not_in_panel"] += int((ok & ~same).sum())

        cols = {"CHROM": a.chrom, "POS": pos, "B73": ref, "ALT": alt, "DP": dp, "OTH": oth,
                "CTRL_REF": cref, "CTRL_ALT": calt, "FILTER": filt, "IN_PANEL": same.astype(np.int8),
                "PANEL_AC": np.where(same, pac[i], -1), "PANEL_AN": np.where(same, pan[i], -1)}
        for g in groups:
            cols[f"{g}_REF"] = R[:, gcols[g]].sum(1); cols[f"{g}_ALT"] = A[:, gcols[g]].sum(1)
        pd.DataFrame(cols, columns=header).to_csv(out, sep="\t", header=False, index=False)

        if ok.any():
            Rp = R[ok][:, acc_perm]; Ap = A[ok][:, acc_perm]
            acc_pos.append(pos[ok])
            acc_ref.append(np.minimum(np.add.reduceat(Rp, acc_starts, axis=1), 65535).astype(np.uint16))
            acc_alt.append(np.minimum(np.add.reduceat(Ap, acc_starts, axis=1), 65535).astype(np.uint16))
    out.close()

    np.savez_compressed(a.out + ".acc_counts.npz", pos=np.concatenate(acc_pos) if acc_pos else np.zeros(0, np.int64),
                        accessions=np.array(accs),
                        ref=np.vstack(acc_ref) if acc_ref else np.zeros((0, len(accs)), np.uint16),
                        alt=np.vstack(acc_alt) if acc_alt else np.zeros((0, len(accs)), np.uint16))

    pc = pools[["SAMPLE", "ROLE", "GROUP", "ACCESSION"]].copy()
    pc["PANEL_REF"], pc["PANEL_ALT"] = cls_ref[0], cls_alt[0]
    pc["NONPANEL_REF"], pc["NONPANEL_ALT"] = cls_ref[1], cls_alt[1]
    for c, cls in enumerate(("PANEL", "NONPANEL")):       # controls get their leave-one-out counts
        pc.loc[is_ctrl, f"{cls}_REF"] = hold_ref[c, is_ctrl]
        pc.loc[is_ctrl, f"{cls}_ALT"] = hold_alt[c, is_ctrl]
    pc.to_csv(a.out + ".pool_class.tsv", sep="\t", index=False)

    s = pd.DataFrame([{"CHROM": a.chrom, **counts, "MEDIAN_DP": med}])
    s.to_csv(a.out + ".summary.tsv", sep="\t", index=False)
    sys.stderr.write(s.T.to_string(header=False) + "\n")
    if counts["sites"] == 0:
        sys.exit("ERROR: no sites read from stdin")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("pools")
    p.add_argument("--samples", required=True); p.add_argument("--metadata", required=True)
    p.add_argument("--controls", required=True); p.add_argument("--exclude", required=True)
    p.add_argument("--out", required=True); p.add_argument("--keep", required=True)
    f = sub.add_parser("filter")
    f.add_argument("--pools", required=True); f.add_argument("--keep", required=True)
    f.add_argument("--sitesum", required=True); f.add_argument("--panel", required=True)
    f.add_argument("--chrom", required=True); f.add_argument("--out", required=True)
    f.add_argument("--max-dp-mult", type=float, default=2.0); f.add_argument("--min-dp-mult", type=float, default=0.25)
    f.add_argument("--max-oth", type=float, default=0.05); f.add_argument("--max-ctrl-alt", type=int, default=0)
    f.add_argument("--chunk", type=int, default=20000)
    a = ap.parse_args()
    cmd_pools(a) if a.cmd == "pools" else cmd_filter(a)


if __name__ == "__main__":
    main()
