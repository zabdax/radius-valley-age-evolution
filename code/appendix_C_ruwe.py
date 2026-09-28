#!/usr/bin/env python3
"""appendix_C_ruwe.py — RUWE + parallax-S/N sensitivity gate.

Requires data/gaia_quality_DR3.csv (see fetch_gaia_quality.py).
Cleaning rule (pre-specified here, before seeing outcomes):
  keep RUWE <= 1.4 AND parallax/parallax_error > 5.
Falls back gracefully if quality file is absent (reports counts only).

For each sample (M, FGK Kepler-only):
  1. report hosts/planets dropped and their vtan distribution vs kept
     (are dropped hosts kinematically hotter?);
  2. refit the real pooled beta on the cleaned base;
  3. run fresh pooled-only arms (null 0, signal -0.14), 200 reps/arm,
     legacy host-effect mode (comparable to production), and report
     arm separation + real percentile + verdict under the operational rule
     (inside-signal-95% AND outside-null-95%).

Outputs results/appendix/: ruwe_clean_summary.json, ruwe_gate_all_raw.csv.
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
DATA = ROOT / "data"

REPS = {"M": 200, "FGK": 200}
SEEDS = {("M", "null"): 7201, ("M", "signal"): 7202,
         ("FGK", "null"): 7211, ("FGK", "signal"): 7212}

_W = {}


def load_cleaned(sample):
    cfg = s15.SAMPLES[sample]
    s15.CFG = cfg
    s15.TARGET_FRAC, s15.VALLEY = cfg["target_frac"], cfg["valley"]
    s15.GAMMA_LOGP = cfg["gamma_logp"]
    base = s15.load_control(cfg)
    n0, h0 = len(base), base.hostname.nunique()
    q = pd.read_csv(DATA / "gaia_quality_DR3.csv")
    g = pd.read_csv(DATA / f"gaia_hosts_{cfg['suffix']}.csv",
                    usecols=["hostname", "source_id"])
    base = base.merge(g, on="hostname", how="left")
    n_nomatch = int(base.source_id.isna().sum())
    base = base.merge(q[["source_id", "ruwe", "parallax", "parallax_error"]],
                      on="source_id", how="left")
    base["psnr"] = base["parallax"] / base["parallax_error"]
    qstats = {"planets_pre": int(n0), "hosts_pre": int(h0),
              "planets_no_source_id": n_nomatch,
              "planets_no_quality_row": int(base["ruwe"].isna().sum()),
              "ruwe_gt1.4": int((base["ruwe"] > 1.4).sum()),
              "psnr_le5": int((base["psnr"] <= 5).sum()),
              "ruwe_median": float(base["ruwe"].median()),
              "psnr_median": float(base["psnr"].median())}
    keep = base[(base["ruwe"] <= 1.4) & (base["psnr"] > 5)].copy()
    drop = base[~base.index.isin(keep.index)].copy()
    qstats["planets_post"] = int(len(keep))
    qstats["hosts_post"] = int(keep.hostname.nunique())
    qstats["planets_dropped"] = int(len(drop))
    qstats["hosts_dropped"] = int(drop.hostname.nunique())
    if len(drop):
        qstats["vtan_median_kept"] = float(keep.vtan.median())
        qstats["vtan_median_dropped"] = float(drop.vtan.median())
    keep = keep.reset_index(drop=True)
    return keep, qstats


def _init(base, comp, beta_fixed, alpha):
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "NUMEXPR_NUM_THREADS"):
        os.environ[v] = "1"
    _W.update(base=base, comp=comp, beta_fixed=beta_fixed, alpha=alpha)


def _one(task):
    rep, seedseq = task
    rng = np.random.default_rng(seedseq)
    d = s15.synth(_W["base"], _W["beta_fixed"], _W["alpha"], rng)
    r = s15.pooled_estimator(d, _W["comp"])
    if r is None:
        return {"rep": rep, "beta_hat": float("nan"), "se": float("nan"),
                "n": len(d)}
    return {"rep": rep, "beta_hat": r[0], "se": r[1], "n": len(d)}


def run_cells(base, comp, alpha, beta_fixed, reps, seed, workers):
    children = np.random.SeedSequence(seed).spawn(reps)
    todo = [(i, children[i]) for i in range(reps)]
    init = (base, comp, beta_fixed, alpha)
    rows = []
    t0 = time.time()
    if workers <= 1:
        _init(*init)
        for k, task in enumerate(todo, 1):
            rows.append(_one(task))
    else:
        chunk = max(1, min(8, len(todo) // (workers * 4) or 1))
        method = ("fork" if "fork" in mp.get_all_start_methods() else "spawn")
        with mp.get_context(method).Pool(
                workers, initializer=_init, initargs=init) as pool:
            for k, row in enumerate(
                    pool.imap_unordered(_one, todo, chunksize=chunk), 1):
                rows.append(row)
                if k % 25 == 0 or k == reps:
                    print(f"    {k}/{reps} {(time.time()-t0)/k:.2f}s/rep",
                          flush=True)
    return pd.DataFrame(rows).sort_values("rep")


def main():
    workers = max(1, (os.cpu_count() or 1) - 1)
    print(f"workers={workers}", flush=True)
    summary, frames = {}, []
    for sample in ("M", "FGK"):
        cfg = s15.SAMPLES[sample]
        base, qs = load_cleaned(sample)
        print(f"[{sample}] {qs}", flush=True)
        full = s15.m9.build_design(base, "vtan", valley=cfg["valley"],
                                   fe=False, weights=base.wgt.values)
        comp = full["init_comp"]
        real = s15.pooled_estimator(base, comp)
        alpha = s15.calibrate_alpha(base, cfg["target_frac"])
        arms = {}
        for arm, beta in (("null", 0.0), ("signal", -0.14)):
            print(f"  [{sample} cleaned {arm}]", flush=True)
            rdf = run_cells(base, comp, alpha, beta, REPS[sample],
                            SEEDS[(sample, arm)], workers)
            rdf.insert(0, "arm", arm)
            rdf.insert(0, "mode", "legacy")
            rdf.insert(0, "sample", sample)
            rdf["beta_inj"] = beta
            rdf.to_csv(OUT / f"ruwe_{sample}_{arm}_raw.csv", index=False)
            frames.append(rdf)
            b = rdf["beta_hat"].values
            b = b[np.isfinite(b)]
            arms[arm] = {"n": int(len(b)), "median": float(np.median(b)),
                         "sd": float(np.std(b, ddof=1)),
                         "q025": float(np.quantile(b, 0.025)),
                         "q975": float(np.quantile(b, 0.975))}
        sep = ((arms["signal"]["median"] - arms["null"]["median"])
               / arms["null"]["sd"])
        p_null_out = not (arms["null"]["q025"] <= real[0]
                          <= arms["null"]["q975"])
        in_sig = (arms["signal"]["q025"] <= real[0]
                  <= arms["signal"]["q975"])
        summary[sample] = {
            "quality": qs,
            "real_beta_cleaned": real[0], "real_se_cleaned": real[1],
            "arms": arms, "arm_separation_nullSD": float(sep),
            "real_pct_null": float(np.mean(
                frames[-2]["beta_hat"].values <= real[0])),
            "PASS": bool(in_sig and p_null_out),
            "verdict": "PASS" if (in_sig and p_null_out)
            else "FAIL(inconclusive)"}
    full = pd.concat(frames, ignore_index=True)
    full.to_csv(OUT / "ruwe_gate_all_raw.csv", index=False)
    summary["meta"] = {"reps": REPS, "rule": "RUWE<=1.4 & parallax_SNR>5",
                       "mode": "legacy", "estimator": "pooled"}
    with open(OUT / "ruwe_clean_summary.json", "w") as fh:
        json.dump(summary, fh, indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
