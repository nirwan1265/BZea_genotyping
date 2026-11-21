"""
Step 37 vs step 29: for each joint-model QTL of a trait, does the stepwise mixed model (37_stepwise_mlmm.R) have it?
  yes            a final stepwise QTL lies in the QTL's support interval (+-2 Mb) or within 10.5 Mb of its lead window
  moved (chrN x) the stepwise model added the QTL there but its refine step moved it along the same chromosome to x Mb
  no             neither
Used by 35_trait_tables.py (Stepwise column). stepwise_dir holds stepwise_<TRAIT>_qtl.tsv and _steps.tsv.
"""
import os
import pandas as pd


def _near(chr_, pos, q):
    return chr_ == q.chr and (q.interval_start_Mb - 2 <= pos <= q.interval_end_Mb + 2 or abs(q.pos_Mb - pos) <= 10.5)


def stepwise_status(Q, stepwise_dir, trait):
    """Q: step-29 qtl2_<TRAIT>_qtl.tsv as a DataFrame (file order = Q1, Q2, ...). Returns a list of strings, or None."""
    fq, fs = f"{stepwise_dir}/stepwise_{trait}_qtl.tsv", f"{stepwise_dir}/stepwise_{trait}_steps.tsv"
    if not os.path.exists(fq):
        return None
    S = pd.read_csv(fq, sep="\t") if os.path.getsize(fq) else pd.DataFrame(columns=["chr", "pos_Mb", "Marker"])
    St = pd.read_csv(fs, sep="\t") if os.path.exists(fs) and os.path.getsize(fs) else pd.DataFrame(columns=["action", "Marker"])
    origin = {}                                  # final marker -> where the forward step added it
    for r in St.itertuples():
        if str(r.action).startswith("move from"):
            src = str(r.action).split("move from ")[1]
            origin[r.Marker] = origin.get(src, src)
    out = []
    for q in Q.itertuples():
        st = "no"
        if any(_near(s.chr, s.pos_Mb, q) for s in S.itertuples()):
            st = "yes"
        else:
            for s in S.itertuples():
                o = origin.get(s.Marker)
                if o and _near(o.split(":")[0], int(o.split(":")[1]) / 1e6, q):
                    st = f"moved ({s.chr} {s.pos_Mb:.1f})"
                    break
        out.append(st)
    return out
