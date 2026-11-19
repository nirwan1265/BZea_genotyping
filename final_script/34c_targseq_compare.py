#!/usr/bin/env python3
"""
34c — Check the k-mer introgression calls against high-depth targeted sequencing of the same plots (steps 34, 34b).

The targeted-sequencing (targseq) samples are named e.g. Zd0010_P3_P3_P3211 = plot Zd.0010_P3_P3_P3.2.1.1, linked to
our pool through name_change_final.txt. At every SNP of the panel, a sample's donor allele fraction is
AF = alt reads / (ref + alt reads); the alternative allele is the non-B73 one (the reference is B73 v5).
No donor genotype is needed. Where the donor carries the B73 allele, AF is 0 in every state, so the tests use
"alt-present" sites (AF >= 0.1) where needed.

Expected AF by our call, if the targseq DNA is the same 6-plant pool / one plant of the plot:
  B73    0                       / 0
  SEG    around 0.5 (0.46, 22b)   / 0, 0.5 or 1 (the plant is AA, AB or BB)
  DONOR  1 at alt-present sites  / 1
The AF distribution inside our SEG windows tells which of the two the targseq samples are.

Site filters (they do not use our calls):
  - depth >= --min_dp in at least --min_samples samples;
  - not alt in most samples: sites with AF >= 0.1 in more than --max_alt_share of the covered samples are dropped
    (paralogs, or places where our B73 differs from the B73 v5 reference); at any one place most pools are B73.
Outputs (--out):
  af_by_state.tsv          AF classes by our call: all sites, and alt-present sites
  per_sample.tsv           per sample: sites, hidden donor in B73 windows, AF >= 0.9 share in DONOR and SEG windows
  window_concordance.tsv   1-Mb bins (one state of ours over the whole bin): our call vs the targseq call
  af_by_state.png          AF histograms by our call
  sites_used.tsv.gz        every sample x site used (for further checks)
Targseq call of a bin: DONOR if >= 2 sites with AF >= 0.9 and more of those than sites with 0.1 <= AF < 0.9;
SEG if >= 2 sites with 0.1 <= AF < 0.9 and more of those; B73 if >= --min_bin_sites covered sites and at most one
alt-present site; else unclear. Where the donor shares most B73 alleles, a DONOR bin can look B73 (few alt sites).

Usage: python final_script/34c_targseq_compare.py --vcf_dir $D34/vcf --tracts $D21/tracts_all.tsv.gz \
         --pool_summary $D21/pool_summary.tsv --names data/phenotypes/raw/name_change_final.txt --out $D34/compare
Needs: numpy, pandas, matplotlib.
"""
import argparse
import glob
import gzip
import os
import re
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--vcf_dir", required=True)
ap.add_argument("--tracts", required=True)
ap.add_argument("--pool_summary", required=True)
ap.add_argument("--names", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--min_dp", type=int, default=10)
ap.add_argument("--min_samples", type=int, default=20)
ap.add_argument("--max_alt_share", type=float, default=0.5)
ap.add_argument("--bin_mb", type=float, default=1.0)
ap.add_argument("--min_bin_sites", type=int, default=5)
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
CODE = {"B73": 0, "SEG": 1, "DONOR": 2}; ST = ["B73", "SEG", "DONOR"]


def plot_name(n):
    """Zd0010_P3_P3_P3211 -> Zd.0010_P3_P3_P3.2.1.1"""
    m = re.match(r"(Z[a-z])(\d{4})_(P\d+)_(P\d+)_P(\d)(\d)(\d)(\d)$", n)
    if not m: return None
    t, acc, f1, bc1, bc2, s1, s2, s3 = m.groups()
    return f"{t}.{acc}_{f1}_{bc1}_P{bc2}.{s1}.{s2}.{s3}"


# ---- VCF: allele depths
rows, samples = [], None
for f in sorted(glob.glob(os.path.join(a.vcf_dir, "targseq_*.vcf.gz"))):
    with gzip.open(f, "rt") as fh:
        for line in fh:
            if line.startswith("##"): continue
            p = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                samples = [os.path.basename(s).replace(".bam", "") for s in p[9:]]; continue
            fmt = p[8].split(":"); i = fmt.index("AD")
            ad = [x.split(":")[i] if len(x.split(":")) > i else "." for x in p[9:]]
            ref = np.array([int(x.split(",")[0]) if x not in (".", "") and "," in x else 0 for x in ad])
            alt = np.array([int(x.split(",")[1]) if x not in (".", "") and "," in x else 0 for x in ad])
            rows.append((p[0], int(p[1]), ref, alt))
print(f"{len(rows):,} SNPs in the VCFs, {len(samples)} samples")
CH = np.array([r[0] for r in rows]); POS = np.array([r[1] for r in rows])
REF = np.vstack([r[2] for r in rows]); ALT = np.vstack([r[3] for r in rows]); DP = REF + ALT
with np.errstate(invalid="ignore", divide="ignore"):
    AF = np.where(DP > 0, ALT / np.maximum(DP, 1), np.nan)
cov = DP >= a.min_dp
alt_share = np.nansum((AF >= 0.1) & cov, axis=1) / np.maximum(cov.sum(1), 1)
keep = (cov.sum(1) >= a.min_samples) & (alt_share <= a.max_alt_share)
print(f"sites with depth >= {a.min_dp} in >= {a.min_samples} samples: {(cov.sum(1) >= a.min_samples).sum():,}; "
      f"after dropping sites alt in > {a.max_alt_share:.0%} of samples: {keep.sum():,}")
CH, POS, AF, DP, cov = CH[keep], POS[keep], AF[keep], DP[keep], cov[keep]

# ---- link samples to pools, look up our state
nm = pd.read_csv(a.names, sep="\t", dtype=str).drop_duplicates("new_genotype").set_index("new_genotype").Sequencing_ID
ps = pd.read_csv(a.pool_summary, sep="\t")
ok = set(ps.loc[(ps.ROLE == "donor") & (ps.FLAG == "ok"), "SAMPLE"])
link = {s: nm.get(plot_name(s)) for s in samples}
use = [j for j, s in enumerate(samples) if link[s] in ok]
print(f"{len(use)} of {len(samples)} samples linked to a pool with k-mer calls")
t = pd.read_csv(a.tracts, sep="\t")
t = t[t.SAMPLE.isin({link[samples[j]] for j in use})]
TR = {k: (g.START.values, g.END.values, g.STATE.map(CODE).values) for k, g in t.sort_values("START").groupby(["SAMPLE", "CHROM"])}


def state_at(sample, ch, pos):
    out = np.full(len(pos), -1)
    if (sample, ch) not in TR: return out
    s, e, c = TR[(sample, ch)]
    i = np.searchsorted(s, pos, side="right") - 1
    okk = (i >= 0) & (pos < e[np.clip(i, 0, None)])
    out[okk] = c[i[okk]]
    return out


recs = []
for j in use:
    sm = samples[j]; pool = link[sm]
    for ch in np.unique(CH):
        m = (CH == ch) & cov[:, j]
        if not m.any(): continue
        st = state_at(pool, ch, POS[m])
        recs.append(pd.DataFrame({"sample": sm, "pool": pool, "chr": ch, "pos": POS[m], "dp": DP[m, j], "af": AF[m, j], "state": st}))
D = pd.concat(recs); D = D[D.state >= 0]
D["call"] = D.state.map(dict(enumerate(ST)))
D.to_csv(os.path.join(a.out, "sites_used.tsv.gz"), sep="\t", index=False, compression="gzip")
CLS = pd.cut(D.af, [0, 0.1, 0.35, 0.65, 0.9, 1.0001], right=False, labels=["<0.1", "0.1-0.35", "0.35-0.65", "0.65-0.9", ">=0.9"])

# ---- 1. AF classes by our call
tab = []
for sub, lab in [(D.index == D.index, "all sites"), (D.af >= 0.1, "alt-present sites (AF >= 0.1)")]:
    x = D[sub]; c = CLS[sub]
    for k in ST:
        y = c[x.call == k]
        tab.append({"sites": lab, "our_call": k, "n": len(y), **(y.value_counts(normalize=True).reindex(c.cat.categories).round(4).to_dict()),
                    "mean_af": round(x.af[x.call == k].mean(), 3)})
T = pd.DataFrame(tab); T.to_csv(os.path.join(a.out, "af_by_state.tsv"), sep="\t", index=False)
print("\nAF classes by our call (share of sample x site):"); print(T.to_string(index=False))

# ---- 2. per sample
ps_rows = []
for (sm, pool), g in D.groupby(["sample", "pool"]):
    b, s, d = g[g.call == "B73"], g[(g.call == "SEG") & (g.af >= 0.1)], g[(g.call == "DONOR") & (g.af >= 0.1)]
    ps_rows.append({"sample": sm, "pool": pool, "sites": len(g), "B73_sites": len(b),
                    "B73_af_ge_0.9": round((b.af >= 0.9).mean(), 4) if len(b) else np.nan,
                    "B73_af_0.1_0.9": round(((b.af >= 0.1) & (b.af < 0.9)).mean(), 4) if len(b) else np.nan,
                    "SEG_altpresent": len(s), "SEG_af_ge_0.9": round((s.af >= 0.9).mean(), 3) if len(s) else np.nan,
                    "SEG_mean_af": round(s.af.mean(), 3) if len(s) else np.nan,
                    "DONOR_altpresent": len(d), "DONOR_af_ge_0.9": round((d.af >= 0.9).mean(), 3) if len(d) else np.nan,
                    "DONOR_mean_af": round(d.af.mean(), 3) if len(d) else np.nan})
P = pd.DataFrame(ps_rows); P.to_csv(os.path.join(a.out, "per_sample.tsv"), sep="\t", index=False)
print(f"\nper sample (medians over {len(P)} samples):")
print(P.drop(columns=["sample", "pool"]).median(numeric_only=True).round(4).to_string())

# ---- 3. bins: our call vs targseq call (bins with one state of ours throughout)
D["bin"] = (D.pos // (a.bin_mb * 1e6)).astype(int)
B = D.groupby(["sample", "chr", "bin"]).agg(n=("af", "size"), ours=("state", "nunique"), state=("state", "first"),
                                            hom=("af", lambda x: (x >= 0.9).sum()), het=("af", lambda x: ((x >= 0.1) & (x < 0.9)).sum())).reset_index()
B = B[B.ours == 1]
B["targseq"] = np.select([(B.hom >= 2) & (B.hom > B.het), (B.het >= 2) & (B.het >= B.hom), (B.n >= a.min_bin_sites) & (B.hom + B.het <= 1)],
                         ["DONOR", "SEG", "B73"], "unclear")
C = pd.crosstab(B.state.map(dict(enumerate(ST))), B.targseq).reindex(index=ST, columns=["B73", "SEG", "DONOR", "unclear"], fill_value=0)
C.to_csv(os.path.join(a.out, "window_concordance.tsv"), sep="\t")
print(f"\n{a.bin_mb:g}-Mb bins, rows = our call, columns = targseq call:"); print(C.to_string())

# ---- figure
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 3, figsize=(12, 3.4), sharey=False)
for k, axx in zip(ST, ax):
    x = D.af[(D.call == k) & (D.af >= 0.1)]
    axx.hist(x, bins=np.linspace(0.1, 1, 37), color="#3b6ea5")
    axx.set_title(f"our call {k}: {len(x):,} alt-present sites\n({(D.call == k).sum():,} sites in all)", fontsize=10, loc="left")
    axx.set_xlabel("donor allele fraction (targseq)")
    axx.spines[["top", "right"]].set_visible(False)
ax[0].set_ylabel("sample x site")
fig.tight_layout(); fig.savefig(os.path.join(a.out, "af_by_state.png"), dpi=200)
print("wrote", a.out)
