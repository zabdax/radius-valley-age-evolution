#!/usr/bin/env python3
"""appendix_A_arm_delta.py — independent-host sensitivity + arm-delta appendix.

Runs POOLED-only calibrated arms (null beta=0, signal beta=-0.14/SD) for the
M dwarf (200 reps/arm) and FGK Kepler-only (100 reps/arm) samples, each under
legacy shared-effect (independent_hosts=False, archived behavior) and fixed
per-occurrence (independent_hosts=True) host-effect assignment.

Outputs (results/appendix/):
  arm_delta_M_raw.csv, arm_delta_FGK_raw.csv (one row per replicate:
    mode, arm, rep, beta_inj, beta_hat, se)
  arm_delta_summary.json (per sample/mode/arm: n, median, sd, q025, q975,
    real percentile vs published real beta, analytical detection rate,
    arm separation (signal-null medians / null SD))

Published real values (finalised draft Table gate):
  M real +0.0458 +- 0.3210 per IQR; FGK real -0.0296 +- 0.1061 per IQR.
"""
import json
import multiprocessing as mp
import os
import pathlib
import sys
import time

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
import importlib.util as _ilu
_spec = _ilu.spec_from_file_location("s15", HERE / "15_stratified_gate.py")
s15 = _ilu.module_from_spec(_spec)
_spec.loader.exec_module(s15)

OUT = ROOT / "results" / "appendix"
OUT.mkdir(parents=True, exist_ok=True)

REAL = {"M": 0.0458, "FGK": -0.0296}
REPS = {"M": 200, "FGK": 100}
B_SIGNAL = -0.14

_W = {}


def _init(base, comp, beta_fixed, alpha, independent_hosts):
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "NUMEXPR_NUM_THREADS"):
        os.environ[v] = "1"
    _W.update(base=base, comp=comp, beta_fixed=beta_fixed, alpha=alpha,
              independent_hosts=independent_hosts)


def _one(task):
    rep, seedseq = task
    rng = np.random.default_rng(seedseq)
    d = s15.synth(_W["base"], _W["beta_fixed"], _W["alpha"], rng,
                  independent_hosts=_W["independent_hosts"])
    r = s15.pooled_estimator(d, _W["comp"])
    if r is None:
        return {"rep": rep, "beta_hat": float("nan"), "se": float("nan"),
                "n": len(d)}
    return {"rep": rep, "beta_hat": r[0], "se": r[1], "n": len(d)}


def run_one(sample, mode, arm, beta_fixed, reps, seed, workers):
    cfg = s15.SAMPLES[sample]
    s15.CFG = cfg
    s15.TARGET_FRAC, s15.VALLEY = cfg["target_frac"], cfg["valley"]
    s15.GAMMA_LOGP = cfg["gamma_logp"]
    base = s15.load_control(cfg)
    full = s15.m9.build_design(base, "vtan", valley=cfg["valley"], fe=False,
                               weights=base.wgt.values)
    comp = full["init_comp"]
    alpha = s15.calibrate_alpha(base, cfg["target_frac"])
    indep = (mode == "independent")
    tag = f"{sample}_{mode}_{arm}"
    print(f"[{tag}] N={len(base)} hosts={base.hostname.nunique()} "
          f"beta_inj={beta_fixed:+.3f} reps={reps} indep={indep}", flush=True)
    children = np.random.SeedSequence(seed).spawn(reps)
    todo = [(i, children[i]) for i in range(reps)]
    t0 = time.time()
    init = (base, comp, beta_fixed, alpha, indep)
    rows = []
    if workers <= 1:
        _init(*init)
        for k, task in enumerate(todo, 1):
            rows.append(_one(task))
            if k % 25 == 0 or k == reps:
                print(f"  [{tag}] {k}/{reps} "
                      f"{(time.time()-t0)/k:.2f}s/rep", flush=True)
    else:
        chunk = max(1, min(8, len(todo) // (workers * 4) or 1))
        method = ("fork" if "fork" in mp.get_all_start_methods() else "spawn")
        with mp.get_context(method).Pool(
                workers, initializer=_init, initargs=init) as pool:
            for k, row in enumerate(
                    pool.imap_unordered(_one, todo, chunksize=chunk), 1):
                rows.append(row)
                if k % 25 == 0 or k == reps:
                    print(f"  [{tag}] {k}/{reps} "
                          f"{(time.time()-t0)/k:.2f}s/rep", flush=True)
    df = pd.DataFrame(rows).sort_values("rep")
    df.insert(0, "arm", arm)
    df.insert(0, "mode", mode)
    df.insert(0, "sample", sample)
    df["beta_inj"] = beta_fixed
    return df


def summarize(df, sample):
    real = REAL[sample]
    out = {}
    sub = df[df["sample"] == sample]
    for mode in ("legacy", "independent"):
        m = {}
        arms = {}
        for arm in ("null", "signal"):
            b = sub[(sub["mode"] == mode) & (sub["arm"] == arm)]["beta_hat"].values
            b = b[np.isfinite(b)]
            se = sub[(sub["mode"] == mode) & (sub["arm"] == arm)]["se"].values
            se = se[np.isfinite(se)]
            det = float(np.mean((b - 1.96 * se) * (b + 1.96 * se) > 0))
            arms[arm] = {"n": int(len(b)), "median": float(np.median(b)),
                         "sd": float(np.std(b, ddof=1)) if len(b) > 1 else 0.0,
                         "q025": float(np.quantile(b, 0.025)),
                         "q975": float(np.quantile(b, 0.975)),
                         "pct_rank_of_real": float(np.mean(b <= real)),
                         "analytical_det_rate": det}
        armsep = ((arms["signal"]["median"] - arms["null"]["median"])
                  / arms["null"]["sd"]) if arms["null"]["sd"] > 0 else nan_(0)
        m["arms"] = arms
        m["arm_separation_nullSD"] = float(armsep)
        m["delta_vs_legacy"] = None
        out[mode] = m
    for arm in ("null", "signal"):
        for stat in ("median", "sd"):
            out["independent"]["delta_vs_legacy"] = out["independent"].get(
                "delta_vs_legacy") or {}
            l = out["legacy"]["arms"][arm][stat]
            f = out["independent"]["arms"][arm][stat]
            out["independent"]["delta_vs_legacy"][f"{arm}_{stat}"] = float(f - l)
    return out


def nan_(x):
    return float("nan")


def main():
    workers = max(1, (os.cpu_count() or 1) - 1)
    print(f"workers={workers}")
    alldfs = []
    # seeds fixed and distinct per cell
    seeds = {("M", "legacy", "null"): 701, ("M", "legacy", "signal"): 702,
             ("M", "independent", "null"): 703,
             ("M", "independent", "signal"): 704,
             ("FGK", "legacy", "null"): 711, ("FGK", "legacy", "signal"): 712,
             ("FGK", "independent", "null"): 713,
             ("FGK", "independent", "signal"): 714}
    for sample in ("M", "FGK"):
        for mode in ("legacy", "independent"):
            for arm, beta in (("null", 0.0), ("signal", B_SIGNAL)):
                df = run_one(sample, mode, arm, beta, REPS[sample],
                             seeds[(sample, mode, arm)], workers)
                alldfs.append(df)
                df.to_csv(OUT / f"arm_delta_{sample}_{mode}_{arm}_raw.csv",
                          index=False)
    full = pd.concat(alldfs, ignore_index=True)
    full.to_csv(OUT / "arm_delta_all_raw.csv", index=False)
    summary = {s: summarize(full, s) for s in ("M", "FGK")}
    summary["meta"] = {"reps": REPS, "B_SIGNAL": B_SIGNAL, "REAL": REAL,
                       "note": "pooled estimator only; components pinned to "
                               "full-sample; TRIM=0.95 via load_control"}
    with open(OUT / "arm_delta_summary.json", "w") as fh:
        json.dump(summary, fh, indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
