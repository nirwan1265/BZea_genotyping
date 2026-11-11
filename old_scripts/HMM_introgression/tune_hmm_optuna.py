#!/usr/bin/env python3
"""
Optuna-based HMM parameter tuning for introgression calling.

This script tunes HMM parameters (R, rho, eta) by maximizing stability
under marker thinning perturbations.

Usage:
    # Set environment variables
    export HMM_PY="/path/to/hmm_viterbi_bc2s3.py"
    export MAPDIR="/path/to/genetic_maps"  # folder with chr1.map, chr2.map, ...
    export POST_R="/path/to/batch_introgression_analysis.R"
    export OUTROOT="optuna_tuning_runs"
    export MAX_GAP_KB=10000
    export MIN_BLOCK_KB=5000
    export DROP_FRAC=0.10
    export K_REPS=3

    python tune_hmm_optuna.py calib_samples.txt
"""
import os, sys, gzip, math, random, shutil, subprocess, traceback
from dataclasses import dataclass
from typing import List, Tuple, Optional

import optuna

# ----------------------------
# USER SETTINGS (edit these)
# ----------------------------
# HMM_PY: path to the HMM script (hmm_viterbi_bc2s3.py)
HMM_PY = os.environ.get("HMM_PY", "/path/to/hmm_viterbi_bc2s3.py")
# MAPDIR: directory containing chr1.map, chr2.map, ..., chr10.map
MAPDIR = os.environ.get("MAPDIR", "/path/to/maps")
# POST_R: path to the batch R analysis script (batch_introgression_analysis.R)
POST_R = os.environ.get("POST_R", "/path/to/batch_introgression_analysis.R")

OUTROOT = os.environ.get("OUTROOT", "optuna_tuning_runs")

# Freeze post-analysis knobs during HMM tuning
MAX_GAP_KB = int(os.environ.get("MAX_GAP_KB", "10000"))
MIN_BLOCK_KB = int(os.environ.get("MIN_BLOCK_KB", "5000"))

# Stability experiment
DROP_FRAC = float(os.environ.get("DROP_FRAC", "0.10"))
K_REPS = int(os.environ.get("K_REPS", "3"))  # use 3 for tuning; validate later with 5-10

# HMM fixed knobs (keep constant while tuning)
PRIOR_RR = 0.85
PRIOR_RH = 0.05
PRIOR_HH = 0.10
MIN_MORGAN = 1e-12
RH_PENALTY = 0.0
MIN_RUN_HH = 3
MIN_RUN_RH = 5
DECODE = "posterior_hysteresis"
ETA_RR_FROM_RH = 0.0

# ----------------------------
# Helpers
# ----------------------------
def run(cmd: List[str], cwd: Optional[str] = None):
    p = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"Command failed:\n{' '.join(cmd)}\nSTDOUT:\n{p.stdout}\nSTDERR:\n{p.stderr}")

def open_text(path: str, mode: str = "rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)

def thin_gl(infile: str, outfile: str, drop_frac: float, seed: int):
    rng = random.Random(seed)
    with open_text(infile, "rt") as f, open_text(outfile, "wt") as out:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if rng.random() < drop_frac:
                continue
            out.write(line + "\n")

def sample_name_from_gl(path: str) -> str:
    """Extract sample name from GL file path."""
    b = os.path.basename(path)
    # Handle various naming conventions
    for suf in [".chr1-10.GL.filtered.tsv.gz", ".chr1-10.filtered.GL.tsv.gz",
                ".GL.filtered.tsv.gz", ".filtered.GL.tsv.gz",
                ".GL.tsv.gz", ".tsv.gz", ".gz"]:
        if b.endswith(suf):
            return b[: -len(suf)]
    return b

def bed_read(path: str):
    rows = []
    with open(path, "rt") as f:
        for line in f:
            line = line.strip()
            if not line: 
                continue
            chrom, start, end, geno = line.split("\t")[:4]
            rows.append((chrom, int(start), int(end), geno))
    rows.sort(key=lambda x: (x[0], x[1], x[2]))
    return rows

def donor_intervals(bed_rows):
    # AB/BB are donor
    out = []
    for chrom, start, end, geno in bed_rows:
        if geno in ("AB", "BB"):
            if end > start:
                out.append((chrom, start, end))
    return out

def intervals_total_len(iv):
    return sum(max(0, e - s) for _, s, e in iv)

def intervals_intersection_len(a, b):
    # two-pointer per chrom
    inter = 0
    ia = 0
    ib = 0
    while ia < len(a) and ib < len(b):
        ca, sa, ea = a[ia]
        cb, sb, eb = b[ib]
        if ca < cb:
            ia += 1
            continue
        if cb < ca:
            ib += 1
            continue
        # same chrom
        s = max(sa, sb)
        e = min(ea, eb)
        if e > s:
            inter += (e - s)
        if ea < eb:
            ia += 1
        else:
            ib += 1
    return inter

def donor_jaccard(base_bed: str, rep_bed: str) -> float:
    A = donor_intervals(bed_read(base_bed))
    B = donor_intervals(bed_read(rep_bed))
    if not A and not B:
        return 1.0
    if not A or not B:
        return 0.0
    inter = intervals_intersection_len(A, B)
    uni = intervals_total_len(A) + intervals_total_len(B) - inter
    return 1.0 if uni <= 0 else inter / uni

def state_concordance(base_statepath_gz: str, rep_statepath_gz: str) -> float:
    # merge-join two sorted statepath files by (CHROM, POS)
    def key(chrom, pos):
        # chroms are chr1..chr10 so this sorts fine lexicographically if consistent
        return (chrom, int(pos))

    with gzip.open(base_statepath_gz, "rt") as fb, gzip.open(rep_statepath_gz, "rt") as fr:
        hb = next(fb, None)
        hr = next(fr, None)
        if hb is None or hr is None:
            return float("nan")

        lb = next(fb, None)
        lr = next(fr, None)

        same = 0
        total = 0

        while lb is not None and lr is not None:
            pb = lb.rstrip("\n").split("\t")
            pr = lr.rstrip("\n").split("\t")

            kb = key(pb[0], pb[1])
            kr = key(pr[0], pr[1])

            if kb == kr:
                sb = pb[5]  # STATE
                sr = pr[5]
                total += 1
                if sb == sr:
                    same += 1
                lb = next(fb, None)
                lr = next(fr, None)
            elif kb < kr:
                lb = next(fb, None)
            else:
                lr = next(fr, None)

        return 0.0 if total == 0 else same / total

def breakpoints_per_mb(complete_bed: str) -> float:
    # count blocks in bed (AA/AB/BB) -> breakpoints = blocks - 1
    rows = bed_read(complete_bed)
    if not rows:
        return 0.0
    blocks = len(rows)
    # genome Mb: use B73v5-ish sum, but constant cancels in tuning; approximate:
    genome_mb = 0.0
    # estimate from bed coverage (sum lengths)
    genome_mb = sum(max(0, e - s) for _, s, e, _ in rows) / 1e6
    if genome_mb <= 0:
        genome_mb = 2000.0
    return max(0.0, (blocks - 1) / genome_mb)

def read_summary_het_pct(introgression_summary_total_tsv: str, sample: str) -> float:
    # file has header: sample Ref_pct Het_pct Teo_pct Donor_pct
    with open(introgression_summary_total_tsv, "rt") as f:
        header = next(f).strip().split("\t")
        col = {h:i for i,h in enumerate(header)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p[col["sample"]] == sample:
                return float(p[col["Het_pct"]])
    return float("nan")

def run_pipeline(gl_path: str, outdir: str, R: int, rho: float, eta: float):
    """
    Run HMM + post-analysis pipeline for a single sample.

    Returns: (sample_name, statepath_file, complete_bed_file, intro_summary_file)
    Raises: RuntimeError if any step fails
    """
    os.makedirs(outdir, exist_ok=True)
    sample = sample_name_from_gl(gl_path)

    statepath = os.path.join(outdir, f"{sample}.statepath.tsv.gz")
    tracts = os.path.join(outdir, f"{sample}.tracts.bed.gz")

    # Step 1: Run HMM
    run([
        "python3", HMM_PY,
        "--in_gl_tsv_gz", gl_path,
        "--map_dir", MAPDIR,
        "--out_statepath_gz", statepath,
        "--out_tracts_bed_gz", tracts,
        "--prior_rr", str(PRIOR_RR),
        "--prior_rh", str(PRIOR_RH),
        "--prior_hh", str(PRIOR_HH),
        "--min_morgan", str(MIN_MORGAN),
        "--eta_hh_from_rh", str(eta),
        "--eta_rr_from_rh", str(ETA_RR_FROM_RH),
        "--rh_penalty", str(RH_PENALTY),
        "--decode", DECODE,
        "--rho", str(rho),
        "--rigidity", str(R),
        "--min_run_hh", str(MIN_RUN_HH),
        "--min_run_rh", str(MIN_RUN_RH),
    ])

    # Verify statepath was created
    if not os.path.exists(statepath):
        raise RuntimeError(f"HMM did not produce statepath: {statepath}")

    # Step 2: Run R post-analysis on that folder (expects .statepath.tsv.gz files)
    # POST_R should be the path to batch_introgression_analysis.R
    run([
        "Rscript", POST_R,
        outdir,              # input folder with .statepath.tsv.gz files
        str(MAX_GAP_KB),     # max gap in kb
        str(MIN_BLOCK_KB),   # min block size in kb
        outdir,              # output folder
    ])

    complete_bed = os.path.join(outdir, f"{sample}.complete_blocks.bed")
    intro_total = os.path.join(outdir, "introgression_summary_total.tsv")

    # Verify outputs exist
    if not os.path.exists(complete_bed):
        raise RuntimeError(f"R script did not produce BED: {complete_bed}")
    if not os.path.exists(intro_total):
        raise RuntimeError(f"R script did not produce summary: {intro_total}")

    return sample, statepath, complete_bed, intro_total

# ----------------------------
# Optuna objective
# ----------------------------
def objective(trial: optuna.Trial, samples: List[str]) -> float:
    """
    Objective function for Optuna optimization.

    Maximizes stability (Jaccard similarity + state concordance) under marker thinning.
    Penalizes excessive fragmentation and low heterozygosity.

    Returns: float score (higher is better)
    Raises: optuna.TrialPruned if trial fails
    """
    # Wide search space (expand if optimum hits boundary)
    R = trial.suggest_int("R", 2, 300, log=True)      # log helps explore small/large
    rho = trial.suggest_float("rho", 0.02, 1.0, log=True)
    eta = trial.suggest_float("eta", 0.0, 0.5)

    print(f"\n[Trial {trial.number}] R={R}, rho={rho:.4f}, eta={eta:.4f}")

    # Evaluate on a subset each trial (for speed)
    use_samples = samples[: min(5, len(samples))]

    scores = []
    n_failed = 0

    for gl in use_samples:
        sname = sample_name_from_gl(gl)
        base_dir = os.path.join(OUTROOT, "trial_tmp", f"{trial.number}", sname, "baseline")
        rep_root = os.path.join(OUTROOT, "trial_tmp", f"{trial.number}", sname, "rep")
        thin_root = os.path.join(OUTROOT, "thin_cache", sname)

        try:
            # clean per-sample trial dirs
            if os.path.exists(os.path.dirname(base_dir)):
                shutil.rmtree(os.path.dirname(base_dir), ignore_errors=True)
            os.makedirs(rep_root, exist_ok=True)
            os.makedirs(thin_root, exist_ok=True)

            # baseline run
            sample, base_state, base_bed, base_intro = run_pipeline(gl, base_dir, R, rho, eta)

            # get het% (penalty guardrail)
            het_pct = read_summary_het_pct(base_intro, sample)
            het_penalty = 0.0
            if math.isfinite(het_pct) and het_pct < 0.05:   # adjust threshold if needed
                het_penalty = (0.05 - het_pct) * 2.0        # mild penalty

            # replicates
            jac = []
            conc = []
            bpmb = []

            for k in range(1, K_REPS + 1):
                thin_gl_path = os.path.join(thin_root, f"drop{int(DROP_FRAC*100)}_seed{k}.GL.tsv.gz")
                if not os.path.exists(thin_gl_path):
                    thin_gl(gl, thin_gl_path, DROP_FRAC, k)

                rep_dir = os.path.join(rep_root, f"rep{k}")
                _, rep_state, rep_bed, _ = run_pipeline(thin_gl_path, rep_dir, R, rho, eta)

                jac.append(donor_jaccard(base_bed, rep_bed))
                conc.append(state_concordance(base_state, rep_state))
                bpmb.append(breakpoints_per_mb(rep_bed))

            m_jac = sum(jac) / len(jac) if jac else 0.0
            m_conc = sum(conc) / len(conc) if conc else 0.0
            m_bpmb = sum(bpmb) / len(bpmb) if bpmb else 0.0

            # Final score: maximize Jaccard + concordance, penalize fragmentation and low het
            score = m_jac + 0.5 * m_conc - 0.02 * m_bpmb - het_penalty
            scores.append(score)
            print(f"  {sname}: Jac={m_jac:.3f}, Conc={m_conc:.3f}, BP/Mb={m_bpmb:.2f}, Het={het_pct:.1f}%, Score={score:.3f}")

        except Exception as e:
            print(f"  {sname}: FAILED - {e}")
            traceback.print_exc()
            n_failed += 1
            continue

    # If all samples failed, prune this trial
    if not scores:
        print(f"[Trial {trial.number}] All samples failed, pruning trial")
        raise optuna.TrialPruned(f"All {n_failed} samples failed")

    # If more than half failed, penalize heavily
    if n_failed > len(use_samples) / 2:
        print(f"[Trial {trial.number}] Many failures ({n_failed}/{len(use_samples)}), penalizing")
        return sum(scores) / len(scores) * 0.5

    final_score = sum(scores) / len(scores)
    print(f"[Trial {trial.number}] Final score: {final_score:.4f}")
    return final_score

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 tune_hmm_optuna.py calib_samples.txt [n_trials]", file=sys.stderr)
        print("\nEnvironment variables:", file=sys.stderr)
        print(f"  HMM_PY={HMM_PY}", file=sys.stderr)
        print(f"  MAPDIR={MAPDIR}", file=sys.stderr)
        print(f"  POST_R={POST_R}", file=sys.stderr)
        print(f"  OUTROOT={OUTROOT}", file=sys.stderr)
        print(f"  MAX_GAP_KB={MAX_GAP_KB}", file=sys.stderr)
        print(f"  MIN_BLOCK_KB={MIN_BLOCK_KB}", file=sys.stderr)
        print(f"  DROP_FRAC={DROP_FRAC}", file=sys.stderr)
        print(f"  K_REPS={K_REPS}", file=sys.stderr)
        sys.exit(2)

    with open(sys.argv[1], "rt") as f:
        samples = [line.strip() for line in f if line.strip() and not line.startswith("#")]

    n_trials = int(sys.argv[2]) if len(sys.argv) > 2 else 50

    # Validate environment
    print("=== HMM Parameter Tuning with Optuna ===")
    print(f"HMM script: {HMM_PY}")
    print(f"Map directory: {MAPDIR}")
    print(f"R post-analysis: {POST_R}")
    print(f"Output root: {OUTROOT}")
    print(f"Samples: {len(samples)}")
    print(f"Trials: {n_trials}")
    print()

    # Check if required files exist
    if not os.path.exists(HMM_PY):
        print(f"ERROR: HMM script not found: {HMM_PY}", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(POST_R):
        print(f"ERROR: R script not found: {POST_R}", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(MAPDIR):
        print(f"ERROR: Map directory not found: {MAPDIR}", file=sys.stderr)
        sys.exit(1)

    # Check for genetic maps
    for c in range(1, 11):
        map_file = os.path.join(MAPDIR, f"chr{c}.map")
        if not os.path.exists(map_file):
            print(f"ERROR: Genetic map not found: {map_file}", file=sys.stderr)
            print(f"  Run: python split_genetic_map.py NAM_genetic_map.txt {MAPDIR}", file=sys.stderr)
            sys.exit(1)

    # Check sample files
    for s in samples:
        if not os.path.exists(s):
            print(f"WARNING: Sample file not found: {s}", file=sys.stderr)

    os.makedirs(OUTROOT, exist_ok=True)

    # Create study with SQLite storage for persistence
    storage_path = os.path.join(OUTROOT, "optuna_study.db")
    storage = f"sqlite:///{storage_path}"

    study = optuna.create_study(
        study_name="hmm_tuning",
        direction="maximize",
        storage=storage,
        load_if_exists=True,
    )

    print(f"Study storage: {storage_path}")
    print(f"Previous trials: {len(study.trials)}")
    print()

    # Run optimization
    study.optimize(
        lambda t: objective(t, samples),
        n_trials=n_trials,
        catch=(Exception,),  # Catch all exceptions, don't crash
    )

    # Print results
    print("\n" + "=" * 60)
    print("OPTIMIZATION COMPLETE")
    print("=" * 60)
    print(f"Best parameters: {study.best_params}")
    print(f"Best score: {study.best_value:.4f}")

    # Save best params to file
    best_params_file = os.path.join(OUTROOT, "best_params.txt")
    with open(best_params_file, "w") as f:
        f.write(f"# Best HMM parameters from Optuna tuning\n")
        f.write(f"# Score: {study.best_value:.4f}\n")
        f.write(f"# Trials: {len(study.trials)}\n\n")
        for k, v in study.best_params.items():
            f.write(f"{k}={v}\n")
    print(f"\nBest params saved to: {best_params_file}")

    # Boundary check hints
    bp = study.best_params
    if bp["R"] >= 280:
        print("\nNOTE: best R near upper bound; consider rerunning with larger max R (e.g., 600).")
    if bp["rho"] >= 0.9:
        print("\nNOTE: best rho near upper bound; consider expanding rho range or rethinking transitions.")
    if bp["eta"] >= 0.48:
        print("\nNOTE: best eta near upper bound; consider expanding eta range.")

    # Print top 5 trials
    print("\nTop 5 trials:")
    try:
        trials_df = study.trials_dataframe()
        if len(trials_df) > 0:
            trials_df = trials_df.sort_values("value", ascending=False).head(5)
            print(trials_df[["number", "value", "params_R", "params_rho", "params_eta"]].to_string())
    except ImportError:
        # pandas not installed
        completed = [t for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE]
        completed.sort(key=lambda t: t.value or 0, reverse=True)
        for t in completed[:5]:
            print(f"  Trial {t.number}: score={t.value:.4f}, R={t.params['R']}, rho={t.params['rho']:.4f}, eta={t.params['eta']:.4f}")

if __name__ == "__main__":
    main()
