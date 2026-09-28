#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
15_stratified_gate.py -- calibrated gate for the STRATIFIED age estimator.

Motivation
----------
The production pooled fit on the weighted Kepler-only FGK control gives
beta_age/IQR = -0.0296 +- 0.1061. Splitting the sample on host brightness
and refitting gives -0.136 (bright half) and -0.170 (faint half); splitting
on distance gives -0.097 (near) and -0.219 (far). Both halves of both
partitions lie further from zero than the pooled fit, which is the signature
of a nuisance term that is removed by stratification. High-v_tan hosts in
this sample are more distant (686 -> 856 pc), fainter (Kp 14.42 -> 14.66)
and have blurrier radii (rel. err 0.149 -> 0.177, MWU p = 6e-6), so radius
measurement error is correlated with the age axis by construction.

That pattern was found post hoc by refitting six slices. It therefore
CANNOT be reported as a result on its own -- it has to be run through the
same calibrated gate as the pooled estimator. This script does that.

Estimator under test
--------------------
  1. split the sample at the median Kepler magnitude
  2. fit the production mixture model separately in each half, with the
     mixture components PINNED to the full-sample values (so the two halves
     use an identical model and the split cannot move the components)
  3. inverse-variance pool the two per-IQR coefficients

The pooled (production) estimator is fitted on the SAME replicate so the two
estimators are compared on identical synthetic data rather than across runs.

Everything else is production-identical to 13_calibration_fgk.py: same
weighting (Kepler-only, inverse-detection, trimmed at the 95th percentile),
same synth() generative process, same alpha calibration to a 0.408 target
sub-Neptune fraction, same B_INJ = -0.14 per SD of v_tan, same fixed_su
working-independence fit with cluster-robust sandwich SEs.

synth() is reimplemented with an index-based host gather instead of the
original per-host concat. This is a pure speed change: the RNG call order
and sizes are unchanged, so the random stream is identical to production.
A --verify flag asserts this against the original implementation.

Usage
-----
  python 15_stratified_gate.py                  # 400 reps/arm (default)
  python 15_stratified_gate.py --reps 1000      # full production
  python 15_stratified_gate.py --reps 20 --quick-check
  python 15_stratified_gate.py --verify         # RNG-equivalence test only
  python 15_stratified_gate.py --strat dist     # sensitivity: distance split
  python 15_stratified_gate.py --strat relerr   # sensitivity: precision split

Checkpoints every 25 replicates to OUT/checkpoints/, resumable. Writes
raw replicate-level CSV (the ground truth) plus a summary JSON.
"""
import argparse
import functools
import importlib.util
import multiprocessing as mp
import os
import json
import pathlib
import sys
import time

import numpy as np
import pandas as pd

# ----------------------------------------------------------------- paths ---
HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent if (HERE.parent / "data").is_dir() else HERE
DATA = ROOT / "data"
OUT = ROOT / "results" / "stratified_gate"
CKPT = OUT / "checkpoints"


def _set_paths(sample):
    global OUT, CKPT
    OUT = ROOT / "results" / "stratified_gate" / sample
    CKPT = OUT / "checkpoints"

_spec = importlib.util.spec_from_file_location(
    "h9", (ROOT / "code" / "09_hierarchical_model.py"))
m9 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m9)

# ---------------------------------------------------- production constants --
B_INJ = -0.14          # per SD of v_tan; see derivation note below
SIGMA_U = 1.5          # injected host random intercept, logits
TRIM = 0.95
SEED = 20260826        # same seed as 13_calibration_fgk.py

# Per-sample configuration. gamma_logP and valley are taken from each sample's
# own archived hierarchical fit rather than shared, because the period
# dependence of the sub-Neptune fraction is ~2x steeper around M dwarfs
# (+4.39) than around FGK hosts (+2.3). Using the FGK value for the M sample
# would understate the period mixing the injections have to reproduce.
SAMPLES = {
    "FGK": dict(suffix="FGK", valley=1.88, gamma_logp=2.3, target_frac=0.408,
                kepler_only=True, strat_default="mag", magcol="sy_kepmag"),
    # M: gamma_logP = 4.391 +- 0.561 and valley = 1.85 from
    # results/hierarchical_M_vtan.json. target_frac is the completeness
    # weighted sub-Neptune fraction of the real M sample (0.234), which is
    # far below its unweighted value (0.408) -- see the ESS warning below.
    "M":   dict(suffix="M", valley=1.85, gamma_logp=4.391, target_frac=0.234,
                kepler_only=False, strat_default="vmag", magcol="sy_vmag"),
}
CFG = SAMPLES["FGK"]          # rebound in main()
TARGET_FRAC = CFG["target_frac"]
VALLEY = CFG["valley"]
GAMMA_LOGP = CFG["gamma_logp"]

# B_INJ derivation (documented here because it belongs in the manuscript):
#   Kamulali et al. 2026 report a super-Earth/sub-Neptune number ratio of
#   0.51 below 3 Gyr and 0.64 above it. As sub-Neptune fractions those are
#   1/(1+0.51) = 0.662 and 1/(1+0.64) = 0.610, a logit difference of
#   0.674 - 0.447 = -0.227. Their 3 Gyr cut splits 136/687, i.e. the 17th
#   percentile; for a standard normal the two group means are then separated
#   by 1.79 SD, giving -0.127 per SD. The production value -0.14 assumes a
#   median split (1.60 SD separation). Both are carried below as B_INJ and
#   B_INJ_ALT so the choice is visible rather than buried.
B_INJ_ALT = -0.127


# ------------------------------------------------------------------ data ---
def load_control(cfg=None, dr25_only=False):
    """Completeness-weighted control sample for the configured host class."""
    cfg = cfg or CFG
    sfx = cfg["suffix"]
    pl = pd.read_csv(DATA / f"planet_sample_{sfx}.csv")
    kin = pd.read_csv(DATA / f"hosts_kinematics_{sfx}.csv")
    wts = pd.read_csv(DATA / f"completeness_weights_{sfx}.csv")
    mag = pd.read_csv(DATA / f"mags_{sfx}.csv")

    df = (pl.merge(kin[["hostname", "vtan", "W"]], on="hostname", how="inner")
            .merge(wts[[c for c in ("pl_name", "w", "p_det", "src")
                        if c in wts.columns]], on="pl_name", how="left")
            .merge(mag[["pl_name", "sy_kepmag", "sy_vmag"]], on="pl_name",
                   how="left"))
    df["wgt"] = df.pop("w").fillna(1.0)
    if cfg["kepler_only"]:
        df = df[df.hostname.str.startswith("Kepler")].copy()
    if dr25_only:
        # Section 3.1 restricts the FGK primary to hosts whose weights derive
        # from real DR25 CDPP rather than the IDEM proxy. That rule was never
        # applied to the M sample, 70% of which is proxy weighted. This flag
        # applies it. Cost on the M sample: N 417 -> 126, ESS 140 -> 117.
        if "src" not in df.columns:
            sys.exit("--dr25-only needs a 'src' column in the weights table")
        df = df[df.src == "DR25_cdpp"].copy()
    df = df.copy()
    df["wgt"] = np.minimum(df.wgt, np.quantile(df.wgt, TRIM))
    df["relerr"] = np.clip(df.pl_radeerr1.abs() / df.pl_rade, 0.01, 0.30)

    # median-impute the stratification variables so no replicate silently
    # drops rows and changes N between arms
    for c in ("sy_kepmag", "sy_vmag", "sy_dist"):
        df[c] = df[c].fillna(df[c].median())
    return df.reset_index(drop=True)


STRAT_COL = {"mag": "sy_kepmag", "vmag": "sy_vmag",
             "dist": "sy_dist", "relerr": "relerr"}


# ------------------------------------------------------------ generative ---
def calibrate_alpha(df, target_frac):
    """Bisection for the intercept giving target_frac mean sub-Neptunes.
    Deterministic; takes no rng (matches production, which passed an unused
    rng argument)."""
    lp = np.log10(df.pl_orbper)
    zp = (lp - lp.mean()) / lp.std()
    lo, hi = -6.0, 6.0
    for _ in range(28):
        mid = 0.5 * (lo + hi)
        frac = np.mean(1.0 / (1.0 + np.exp(-(mid + GAMMA_LOGP * zp))))
        if frac < target_frac:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def synth(base, beta_inj, alpha, rng, _idx_cache={}, independent_hosts=False):
    """Host-bootstrap + inject. RNG-stream-identical to production synth().
    independent_hosts=True: one draw per bootstrap OCCURRENCE (fixed form).
    False (default): legacy shared-effect form for comparability."""
    key = id(base)
    if key not in _idx_cache:
        hosts_seen = pd.unique(base.hostname)
        codes = pd.Categorical(base.hostname, categories=hosts_seen).codes
        order = np.argsort(codes, kind="stable")
        starts = np.searchsorted(codes[order], np.arange(codes.max() + 1))
        ends = np.r_[starts[1:], len(order)]
        _idx_cache[key] = (hosts_seen, order, starts, ends)
    hosts, order, starts, ends = _idx_cache[key]

    # --- RNG call 1: host bootstrap (same call, same size as production) ---
    picked = rng.choice(base.hostname.unique(), size=len(hosts), replace=True)

    hpos = {h: i for i, h in enumerate(hosts)}
    sel = [hpos[h] for h in picked]
    take = np.concatenate([order[starts[i]:ends[i]] for i in sel])
    d = base.iloc[take].reset_index(drop=True)
    n = len(d)
    # Host-effect assignment. Legacy production pos={h:i...} maps duplicates
    # to LAST index (shared effect, wasted draws) — preserved by default for
    # comparability; independent_hosts=True uses per-occurrence form.
    uh = rng.normal(0, SIGMA_U, len(picked))
    if independent_hosts:
        sizes = np.array([ends[i] - starts[i] for i in sel])
        eta_add = np.repeat(uh, sizes)
    else:
        pos = {h: i for i, h in enumerate(picked)}
        rep_of = d.hostname.map(pos).values
        eta_add = uh[rep_of]

    zv = (d.vtan.values - base.vtan.mean()) / base.vtan.std()
    lp = np.log10(base.pl_orbper)
    zp = (np.log10(d.pl_orbper.values) - lp.mean()) / lp.std()
    eta = alpha + beta_inj * zv + GAMMA_LOGP * zp

    # --- RNG calls 2-5: identical order and sizes to production ---
    eta += eta_add
    pi = 1.0 / (1.0 + np.exp(-np.clip(eta, -30, 30)))
    sn = rng.uniform(size=n) < pi
    r_true = np.where(sn, rng.normal(2.10, 0.35, n), rng.normal(1.20, 0.18, n))

    e = d.pl_radeerr1.abs().values
    e = np.where(np.isfinite(e) & (e > 0.01), e, 0.05 * d.pl_rade.values)
    e = np.clip(e, 0.01 * d.pl_rade.values, 0.30 * d.pl_rade.values)
    d["pl_rade"] = r_true + rng.normal(0, 1, n) * e
    return d


# ------------------------------------------------------------ estimators ---
def _fit_beta(d, comp=None):
    """Production fit -> (beta_per_IQR, sandwich_se_per_IQR). None on failure."""
    try:
        mdl = m9.build_design(d, "vtan", valley=VALLEY, fe=False,
                              weights=d.wgt.values)
        if comp is not None:
            mdl["init_comp"] = comp
        th, _, _, _ = m9.fit(mdl, fixed_su=True, fix_comp=True)
        se = m9.sandwich_se(mdl, th, fixed_su=True)
        sd = d.vtan.std()
        iqr = d.vtan.quantile(0.75) - d.vtan.quantile(0.25)
        if not np.isfinite(sd) or sd <= 0:
            return None
        k = iqr / sd
        b, s = th[1] * k, se[1] * k
        if not (np.isfinite(b) and np.isfinite(s) and s > 0):
            return None
        return float(b), float(s)
    except (np.linalg.LinAlgError, ValueError, FloatingPointError):
        return None


def pooled_estimator(d, comp):
    return _fit_beta(d, comp)


def stratified_estimator(d, comp, strat_col, cut):
    """Fit each half at a FIXED cut, then inverse-variance pool."""
    lo = d[d[strat_col] <= cut]
    hi = d[d[strat_col] > cut]
    if min(len(lo), len(hi)) < 200:
        return None
    a, b = _fit_beta(lo, comp), _fit_beta(hi, comp)
    if a is None or b is None:
        return None
    est = np.array([a[0], b[0]])
    se = np.array([a[1], b[1]])
    w = 1.0 / se ** 2
    pooled = float((est * w).sum() / w.sum())
    pooled_se = float(1.0 / np.sqrt(w.sum()))
    return pooled, pooled_se, a[0], a[1], b[0], b[1]


# ----------------------------------------------------------------- gate ----
def percentile_of(value, draws):
    """Empirical percentile of value within draws, plus a binomial 95% CI."""
    draws = np.asarray(draws, float)
    draws = draws[np.isfinite(draws)]
    n = len(draws)
    if n == 0:
        return float("nan"), (float("nan"), float("nan")), 0
    p = float(np.mean(draws <= value))
    se = np.sqrt(max(p * (1 - p), 1e-12) / n)
    return p, (max(0.0, p - 1.96 * se), min(1.0, p + 1.96 * se)), n


_W = {}


def _init_worker(base, comp, strat_col, cut, beta_inj, alpha, independent_hosts=False):
    """Runs once per worker. Pin BLAS to one thread: each worker is already a
    process, and letting OpenBLAS spawn its own threads inside every worker
    oversubscribes the CPU and makes the pool SLOWER than serial."""
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ[v] = "1"
    _W.update(base=base, comp=comp, strat_col=strat_col, cut=cut,
              beta_inj=beta_inj, alpha=alpha,
              independent_hosts=independent_hosts)


def _one_rep(task):
    """One replicate. `task` carries its own SeedSequence, so the result for
    replicate i is identical no matter which worker runs it, how many workers
    there are, or where a resume picked up."""
    rep, seedseq = task
    rng = np.random.default_rng(seedseq)
    d = synth(_W["base"], _W["beta_inj"], _W["alpha"], rng,
              independent_hosts=_W.get("independent_hosts", False))
    pol = pooled_estimator(d, _W["comp"])
    strat = stratified_estimator(d, _W["comp"], _W["strat_col"], _W["cut"])
    return {"rep": rep,
            "beta_pooled": pol[0] if pol else np.nan,
            "se_pooled": pol[1] if pol else np.nan,
            "beta_strat": strat[0] if strat else np.nan,
            "se_strat": strat[1] if strat else np.nan,
            "beta_lo": strat[2] if strat else np.nan,
            "beta_hi": strat[4] if strat else np.nan,
            "n": len(d)}


def run_arm(base, comp, strat_col, cut, beta_inj, reps, seed, tag,
            workers=1, legacy_stream=False, independent_hosts=False):
    """One injection arm.

    Two RNG modes:
      legacy_stream=True  -- one sequential Generator, bit-identical to
                             13_calibration_fgk.py. Serial only; resuming
                             replays the stream, which costs O(done) synths.
      legacy_stream=False -- SeedSequence(seed).spawn(reps), one independent
                             substream per replicate (default). Parallel-safe,
                             worker-count-independent, and resume is free.
    """
    ck = CKPT / f"{tag}.json"
    done, rows = 0, []
    if ck.exists():
        blob = json.loads(ck.read_text())
        rows, done = blob["rows"], blob["done"]
        if blob.get("legacy_stream", False) != legacy_stream:
            sys.exit(f"checkpoint {tag} used a different RNG mode; "
                     f"delete it or match the flag")
        print(f"  [{tag}] resuming from replicate {done}")

    alpha = calibrate_alpha(base, TARGET_FRAC)
    t0 = time.time()

    def flush(k):
        CKPT.mkdir(parents=True, exist_ok=True)
        ck.write_text(json.dumps({"done": k, "rows": rows,
                                  "legacy_stream": legacy_stream}))
        rate = (time.time() - t0) / max(k - done, 1)
        print(f"  [{tag}] {k}/{reps}  {rate:.2f}s/rep  "
              f"ETA {rate * (reps - k) / 60:.0f} min", flush=True)

    if legacy_stream:
        rng = np.random.default_rng(seed)
        for _ in range(done):
            synth(base, beta_inj, alpha, rng,
                  independent_hosts=independent_hosts)
        for rep in range(done, reps):
            d = synth(base, beta_inj, alpha, rng,
                      independent_hosts=independent_hosts)
            pol = pooled_estimator(d, comp)
            st = stratified_estimator(d, comp, strat_col, cut)
            rows.append({"rep": rep,
                         "beta_pooled": pol[0] if pol else np.nan,
                         "se_pooled": pol[1] if pol else np.nan,
                         "beta_strat": st[0] if st else np.nan,
                         "se_strat": st[1] if st else np.nan,
                         "beta_lo": st[2] if st else np.nan,
                         "beta_hi": st[4] if st else np.nan,
                         "n": len(d)})
            if (rep + 1) % 25 == 0 or rep + 1 == reps:
                flush(rep + 1)
        return pd.DataFrame(rows)

    children = np.random.SeedSequence(seed).spawn(reps)
    todo = [(i, children[i]) for i in range(done, reps)]
    if not todo:
        return pd.DataFrame(rows)

    init = (base, comp, strat_col, cut, beta_inj, alpha, independent_hosts)
    if workers <= 1:
        _init_worker(*init)
        for k, task in enumerate(todo, start=done + 1):
            rows.append(_one_rep(task))
            if k % 25 == 0 or k == reps:
                flush(k)
    else:
        chunk = max(1, min(8, len(todo) // (workers * 4) or 1))
        method = ("fork" if "fork" in mp.get_all_start_methods()
                  else "spawn")  # fork is POSIX-only; Windows must spawn
        with mp.get_context(method).Pool(
                workers, initializer=_init_worker, initargs=init) as pool:
            for k, row in enumerate(
                    pool.imap_unordered(_one_rep, todo, chunksize=chunk),
                    start=done + 1):
                rows.append(row)
                if k % 25 == 0 or k == reps:
                    rows.sort(key=lambda r: r["rep"])
                    flush(k)
    rows.sort(key=lambda r: r["rep"])
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=400)
    ap.add_argument("--sample", choices=list(SAMPLES), default="FGK")
    ap.add_argument("--strat", choices=list(STRAT_COL), default=None)
    ap.add_argument("--dr25-only", action="store_true",
                    help="keep only planets whose completeness weight comes "
                         "from real DR25 CDPP (drops the IDEM proxy planets)")
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--binj", type=float, default=B_INJ)
    ap.add_argument("--quick-check", action="store_true")
    ap.add_argument("--verify", action="store_true")
    ap.add_argument("--independent-hosts", action="store_true",
                    help="per-occurrence host effects (fixed form); default "
                         "preserves legacy shared-effect stream")
    ap.add_argument("--workers", type=int, default=0,
                    help="0 = auto (cpu_count-1); 1 = serial")
    ap.add_argument("--legacy-stream", action="store_true",
                    help="sequential RNG, bit-identical to production; serial only")
    args = ap.parse_args()

    global CFG, TARGET_FRAC, VALLEY, GAMMA_LOGP
    CFG = SAMPLES[args.sample]
    TARGET_FRAC, VALLEY = CFG["target_frac"], CFG["valley"]
    GAMMA_LOGP = CFG["gamma_logp"]
    if args.strat is None:
        args.strat = CFG["strat_default"]
    _set_paths(args.sample + ("_dr25" if args.dr25_only else ""))

    if args.workers == 0:
        args.workers = max(1, (os.cpu_count() or 1) - 1)
    if args.legacy_stream and args.workers > 1:
        print("--legacy-stream is serial-only; forcing --workers 1")
        args.workers = 1
    if args.workers == 1:
        for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
            os.environ.setdefault(v, "1")
    print(f"workers={args.workers}  rng="
          f"{'legacy sequential' if args.legacy_stream else 'spawned per-replicate'}")
    OUT.mkdir(parents=True, exist_ok=True)
    CKPT.mkdir(parents=True, exist_ok=True)

    base = load_control(CFG, dr25_only=args.dr25_only)
    if args.dr25_only:
        # the configured target_frac describes the full weighted sample; the
        # DR25 subset has a different weighted sub-Neptune fraction, so the
        # synthetic intercept must be recalibrated to the sample actually used
        _sn = (base.pl_rade > VALLEY).astype(float)
        TARGET_FRAC = float(np.average(_sn, weights=base.wgt))
        print(f"--dr25-only: target_frac recomputed from this subset "
              f"= {TARGET_FRAC:.4f} (config value {CFG['target_frac']})")
    strat_col = STRAT_COL[args.strat]
    cov = base[strat_col].notna().mean()
    if cov < 0.9:
        print(f"WARNING: {strat_col} covers only {cov*100:.0f}% of "
              f"{args.sample}; median-imputed rows will not stratify cleanly")
    cut = float(base[strat_col].median())
    ess = base.wgt.sum() ** 2 / (base.wgt ** 2).sum()
    ret = ess / len(base)
    print(f"control [{args.sample}]: N={len(base)}  "
          f"hosts={base.hostname.nunique()}  ESS={ess:.0f} "
          f"({ret*100:.0f}% retained, SE inflation {1/np.sqrt(ret):.2f}x)")
    print(f"valley={VALLEY}  gamma_logP={GAMMA_LOGP}  "
          f"target_frac={TARGET_FRAC}")
    if ret < 0.5:
        print("  !! ESS retention below 50%: the injected ensembles operate "
              "in a much thinner information regime than the FGK control")
    print(f"stratify on {strat_col} at median {cut:.4f}")

    full = m9.build_design(base, "vtan", valley=VALLEY, fe=False,
                           weights=base.wgt.values)
    comp = full["init_comp"]
    print(f"components pinned to full-sample: "
          f"{[round(c, 4) for c in comp]}")

    if args.verify:
        return verify_stream(base)

    # ------------------------------------------------------- real values ---
    real_pooled = pooled_estimator(base, comp)
    real_strat = stratified_estimator(base, comp, strat_col, cut)
    print(f"\nREAL pooled      beta/IQR = {real_pooled[0]:+.4f} "
          f"+- {real_pooled[1]:.4f}")
    if real_strat is None:
        n_lo = int((base[strat_col] <= cut).sum())
        print(f"REAL stratified  NOT ESTIMABLE: halves are {n_lo} and "
              f"{len(base)-n_lo}, below the 200/half floor. Reporting the "
              f"pooled gate only; stratification is not meaningful at this N.")
        estimators = ["pooled"]
    else:
        print(f"REAL stratified  beta/IQR = {real_strat[0]:+.4f} "
              f"+- {real_strat[1]:.4f}   "
              f"(low {real_strat[2]:+.4f}, high {real_strat[4]:+.4f})")
        estimators = ["pooled", "strat"]

    # ------------------------------------------------------------- arms ---
    arms = {}
    for name, binj in [("null", 0.0), ("signal", args.binj)]:
        print(f"\n=== arm: {name} (beta_inj = {binj:+.3f}/SD) ===")
        arms[name] = run_arm(base, comp, strat_col, cut, binj, args.reps,
                             args.seed + (0 if name == "null" else 1),
                             f"{args.sample}_{args.strat}_{name}_r{args.reps}",
                             workers=args.workers,
                             legacy_stream=args.legacy_stream,
                             independent_hosts=args.independent_hosts)
        arms[name].to_csv(
            OUT / f"strat_gate_{args.sample}_{args.strat}_{name}"
            f"_r{args.reps}_raw.csv",
            index=False)

    # ------------------------------------------------------- summarise ---
    summary = {"sample": args.sample, "dr25_only": args.dr25_only,
               "strat_on": strat_col,
               "valley": VALLEY, "gamma_logP": GAMMA_LOGP, "cut": cut, "reps": args.reps,
               "b_inj": args.binj, "seed": args.seed, "N": int(len(base)),
               "ESS": float(ess),
               "real": {"pooled": real_pooled[0], "pooled_se": real_pooled[1],
                        "strat": real_strat[0] if real_strat else None,
                        "strat_se": real_strat[1] if real_strat else None,
                        "strat_low_half": real_strat[2] if real_strat else None,
                        "strat_high_half": real_strat[4] if real_strat else None}}

    print("\n" + "=" * 72)
    _reals = {"pooled": real_pooled[0],
              "strat": real_strat[0] if real_strat else None}
    for est in estimators:
        real = _reals[est]
        col = f"beta_{est}"
        nul, sig = arms["null"][col].values, arms["signal"][col].values
        pn, cin, nn = percentile_of(real, nul)
        ps, cis, ns = percentile_of(real, sig)
        sep = (np.nanmedian(sig) - np.nanmedian(nul)) / np.nanstd(nul)
        verdict = "PASS" if (pn <= 0.05 and 0.025 <= ps <= 0.975) else "FAIL"
        summary[est] = {
            "real": real,
            "null_median": float(np.nanmedian(nul)),
            "null_sd": float(np.nanstd(nul)),
            "null_95": [float(np.nanpercentile(nul, 2.5)),
                        float(np.nanpercentile(nul, 97.5))],
            "signal_median": float(np.nanmedian(sig)),
            "signal_sd": float(np.nanstd(sig)),
            "signal_95": [float(np.nanpercentile(sig, 2.5)),
                          float(np.nanpercentile(sig, 97.5))],
            "pct_in_null": pn, "pct_in_null_ci": list(cin),
            "pct_in_signal": ps, "pct_in_signal_ci": list(cis),
            "arm_separation_sd": float(sep),
            "n_null": nn, "n_signal": ns,
            "recovery": float(np.nanmedian(sig) / (args.binj * 1.236)),
            "verdict": verdict,
        }
        print(f"\n--- estimator: {est} ---")
        print(f"  real                 {real:+.4f}")
        print(f"  null   median {np.nanmedian(nul):+.4f}  sd {np.nanstd(nul):.4f}"
              f"  95% [{np.nanpercentile(nul,2.5):+.3f}, "
              f"{np.nanpercentile(nul,97.5):+.3f}]")
        print(f"  signal median {np.nanmedian(sig):+.4f}  sd {np.nanstd(sig):.4f}"
              f"  95% [{np.nanpercentile(sig,2.5):+.3f}, "
              f"{np.nanpercentile(sig,97.5):+.3f}]")
        print(f"  real percentile in null   {pn*100:5.1f}% "
              f"[{cin[0]*100:.1f}, {cin[1]*100:.1f}]")
        print(f"  real percentile in signal {ps*100:5.1f}% "
              f"[{cis[0]*100:.1f}, {cis[1]*100:.1f}]")
        print(f"  arm separation {sep:.2f} null-SD   recovery "
              f"{summary[est]['recovery']*100:.0f}%   -> {verdict}")

    (OUT / f"strat_gate_{args.sample}_{args.strat}_r{args.reps}"
     "_summary.json").write_text(
        json.dumps(summary, indent=2))
    print(f"\nwrote {OUT}")


def verify_stream(base):
    """Assert the fast synth() is RNG-identical to the production version."""
    src = (ROOT / "code" / "13_calibration_fgk.py")
    s = importlib.util.spec_from_file_location("c13", src)
    c13 = importlib.util.module_from_spec(s)
    sys.modules["c13"] = c13
    s.loader.exec_module(c13)
    alpha = calibrate_alpha(base, TARGET_FRAC)
    ok = True
    for i in range(3):
        a = synth(base, -0.14, alpha, np.random.default_rng(1000 + i))
        b = c13.synth(base, -0.14, alpha, np.random.default_rng(1000 + i))
        same = (len(a) == len(b)
                and np.allclose(np.sort(a.pl_rade.values),
                                np.sort(b.pl_rade.values), atol=1e-12))
        print(f"  replicate {i}: n {len(a)} vs {len(b)}  radii match: {same}")
        ok &= same
    print("RNG-equivalent" if ok else "STREAM MISMATCH -- do not use")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main() or 0)