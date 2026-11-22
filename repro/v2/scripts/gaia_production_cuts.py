#!/usr/bin/env python3
"""gaia_production_cuts.py — production Gaia quality cuts for the radius-valley paper.

WHAT THIS DOES
--------------
1. VERIFY the archived *post-hoc* RUWE / parallax-S/N sensitivity run
   (frozen results/appendix/ruwe_clean_summary.json + ruwe_*_raw.csv):
   recompute every summary statistic from the raw replicate CSVs and from the
   frozen base tables + Gaia quality table.
2. PRODUCTION CUT: apply the pre-specified rule
       keep  RUWE <= 1.4  AND  parallax_SNR = parallax / parallax_error > 5
   to the frozen samples *before* the gate, then recompute the sample counts,
   the completeness-weight effective sample size (ESS) and the real pooled
   beta/SE using the identical weighting, 95%-weight trim, design matrix and
   fixed-component sandwich-SE estimator as the frozen pipeline.
3. PARALLAX ERROR PROPAGATION: for d = 1/parallax, sigma_d/d = sigma_pi/pi;
   report median / 95th percentile over all Gaia hosts and over the SNR<=5
   failures.

READ-ONLY POLICY: the frozen project is never written to. All outputs go to
this manuscript's results/audit/ directory.

USAGE
-----
  python gaia_production_cuts.py                       # verify + recompute (fast)
  python gaia_production_cuts.py --gate-reps 200       # ALSO rerun pooled arms

Multiprocessing (optional arm rerun) is guarded for Windows with the
``if __name__ == '__main__':`` pattern below.
"""
import argparse
import importlib.util
import json
import multiprocessing as mp
import os
import pathlib
import sys
import time

import numpy as np
import pandas as pd

DEFAULT_FROZEN = pathlib.Path(r"REPO\radius-valley-age-evolution")
HERE = pathlib.Path(__file__).resolve().parent
DEFAULT_OUT = HERE.parent / "results" / "audit"

RULE = "RUWE <= 1.4 AND parallax/parallax_error > 5"
RUWE_MAX = 1.4
SNR_MIN = 5.0

# ---------------------------------------------------------------- frozen code
_S15 = None


def load_frozen(frozen):
    """Import the frozen 15_stratified_gate.py (and its m9) read-only."""
    global _S15
    if _S15 is None:
        spec = importlib.util.spec_from_file_location("s15_frozen", frozen / "code" / "15_stratified_gate.py")
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        _S15 = mod
    return _S15


def set_cfg(s15, sample):
    cfg = s15.SAMPLES[sample]
    s15.CFG = cfg
    s15.TARGET_FRAC, s15.VALLEY = cfg["target_frac"], cfg["valley"]
    s15.GAMMA_LOGP = cfg["gamma_logp"]
    return cfg


def ess(w):
    """Kish effective sample size — identical to s15 (line 469) and 14_one_per_host."""
    w = np.asarray(w, float)
    return float(w.sum() ** 2 / (w ** 2).sum())


# ---------------------------------------------------------- quality / cut ---
def quality_table_check(frozen):
    """DR3.csv (used by the archived appendix) vs DR3_v2.csv (frozen v2 table)."""
    a = pd.read_csv(frozen / "data" / "gaia_quality_DR3.csv")
    b = pd.read_csv(frozen / "data" / "gaia_quality_DR3_v2.csv")
    m = a.merge(b, on="source_id", suffixes=("_a", "_b"))
    diffs = {c: float(np.abs(m[f"{c}_a"] - m[f"{c}_b"]).max())
             for c in ("ruwe", "parallax", "parallax_error")}
    return {"rows_DR3": int(len(a)), "rows_DR3_v2": int(len(b)),
            "columns_identical": list(a.columns) == list(b.columns),
            "max_abs_diff": diffs,
            "identical": bool(len(a) == len(b) and all(v == 0.0 for v in diffs.values()))}


def build_cut_frame(s15, frozen, sample, quality_file="gaia_quality_DR3_v2.csv"):
    """Frozen base -> merged frame with RUWE/parallax-SNR, plus pre-cut frame."""
    cfg = set_cfg(s15, sample)
    base = s15.load_control(cfg)                       # applies Kepler-only + 95% weight trim
    q = pd.read_csv(frozen / "data" / quality_file)
    g = pd.read_csv(frozen / "data" / f"gaia_hosts_{cfg['suffix']}.csv",
                    usecols=["hostname", "source_id"])
    b = base.merge(g, on="hostname", how="left")
    n_nomatch = int(b.source_id.isna().sum())
    b = b.merge(q[["source_id", "ruwe", "parallax", "parallax_error"]],
                on="source_id", how="left")
    b["psnr"] = b["parallax"] / b["parallax_error"]
    keep = b[(b.ruwe <= RUWE_MAX) & (b.psnr > SNR_MIN)].copy().reset_index(drop=True)
    stats = {
        "planets_pre": int(len(base)), "hosts_pre": int(base.hostname.nunique()),
        "planets_no_source_id": n_nomatch,
        "planets_no_quality_row": int(b.ruwe.isna().sum()),
        "ruwe_gt1.4": int((b.ruwe > RUWE_MAX).sum()),
        "psnr_le5": int((b.psnr <= SNR_MIN).sum()),
        "dropped_union": int(((b.ruwe > RUWE_MAX) | (b.psnr <= SNR_MIN)).sum()),
        "dropped_intersection": int(((b.ruwe > RUWE_MAX) & (b.psnr <= SNR_MIN)).sum()),
        "ruwe_median": float(b.ruwe.median()), "psnr_median": float(b.psnr.median()),
        "planets_post": int(len(keep)), "hosts_post": int(keep.hostname.nunique()),
        "planets_dropped": int(len(base) - len(keep)),
        "hosts_dropped": int(base.hostname.nunique() - keep.hostname.nunique()),
        "vtan_median_kept": float(keep.vtan.median()),
        "vtan_median_dropped": float(b.loc[~b.index.isin(
            b[(b.ruwe <= RUWE_MAX) & (b.psnr > SNR_MIN)].index), "vtan"].median()),
    }
    return base, keep, stats


def real_beta(s15, frame, valley):
    """Production IDENTICAL: rebuild design on this frame, pin its own components,
    refit real pooled beta and cluster-robust sandwich SE (per IQR)."""
    m9 = s15.m9
    full = m9.build_design(frame, "vtan", valley=valley, fe=False,
                           weights=frame.wgt.values)
    return s15.pooled_estimator(frame, full["init_comp"])


# ------------------------------------------------------- archive verification
def verify_archive(frozen):
    """Recompute arm statistics and separation from the archived raw CSVs."""
    arch = json.loads((frozen / "results" / "appendix" / "ruwe_clean_summary.json").read_text())
    out = {}
    for sample in ("M", "FGK"):
        arms, arrays = {}, {}
        for arm in ("null", "signal"):
            raw = pd.read_csv(frozen / "results" / "appendix" / f"ruwe_{sample}_{arm}_raw.csv")
            b = raw.beta_hat.values[np.isfinite(raw.beta_hat.values)]
            arrays[arm] = b
            arms[arm] = {"n": int(len(b)), "median": float(np.median(b)),
                         "sd": float(np.std(b, ddof=1)),
                         "q025": float(np.quantile(b, 0.025)),
                         "q975": float(np.quantile(b, 0.975))}
        sep = float((np.median(arrays["signal"]) - np.median(arrays["null"]))
                    / np.std(arrays["null"], ddof=1))
        rp = float(np.mean(arrays["null"] <= arch[sample]["real_beta_cleaned"]))
        out[sample] = {
            "arms": arms, "arm_separation_nullSD": sep, "real_pct_null": rp,
            "archive_arms": arch[sample]["arms"],
            "archive_arm_separation_nullSD": arch[sample]["arm_separation_nullSD"],
            "archive_real_pct_null": arch[sample]["real_pct_null"],
            "max_abs_diff_arm_stats": float(max(
                abs(arms[a][k] - arch[sample]["arms"][a][k])
                for a in ("null", "signal")
                for k in ("median", "sd", "q025", "q975"))),
            "sep_diff": float(abs(sep - arch[sample]["arm_separation_nullSD"])),
            "archive_verdict": arch[sample]["verdict"]}
    # baseline (no-cut) production arm separations from confirm_summary.json
    confirm = json.loads((frozen / "results" / "appendix" / "confirm_summary.json").read_text())
    out["baseline_production_arm"] = confirm["meta"]["production_legacy"]
    return out


# -------------------------------------------------- parallax propagation ---
def parallax_quantiles(frozen):
    """sigma_d/d ~ sigma_pi/pi = parallax_error/parallax (small-angle, d=1/pi).
    Uses the DR3 v2 quality table; cross-checks parallax against gaia_hosts_*.csv."""
    q = pd.read_csv(frozen / "data" / "gaia_quality_DR3_v2.csv")
    q["snr"] = q.parallax / q.parallax_error
    q["frac"] = np.abs(q.parallax_error / q.parallax)   # magnitude of fractional d error
    gm = pd.read_csv(frozen / "data" / "gaia_hosts_M.csv",
                     usecols=["source_id", "parallax"]); gm["cls"] = "M"
    gf = pd.read_csv(frozen / "data" / "gaia_hosts_FGK.csv",
                     usecols=["source_id", "parallax"]); gf["cls"] = "FGK"
    hosts = pd.concat([gm, gf], ignore_index=True).drop_duplicates("source_id")
    hq = q.merge(hosts[["source_id", "cls", "parallax"]], on="source_id",
                 suffixes=("_q", "_h"))

    def qtile(sub):
        f = sub.frac.dropna().values
        return {"n": int(len(sub)), "median": float(np.median(f)),
                "p95": float(np.quantile(f, 0.95)), "p05": float(np.quantile(f, 0.05)),
                "max": float(f.max())}

    return {
        "definition": "d=1/pi (mas->kpc); sigma_d/d = sigma_pi/pi = parallax_error/parallax",
        "quality_table": "data/gaia_quality_DR3_v2.csv",
        "quality_rows": int(len(q)), "unique_source_ids": int(q.source_id.nunique()),
        "nonpositive_parallax": int((q.parallax <= 0).sum()),
        "parallax_vs_gaia_hosts_max_abs_diff": float(
            np.abs(hq.parallax_q - hq.parallax_h).max()),
        "all_hosts": qtile(q),
        "M_hosts": qtile(hq[hq.cls == "M"]),
        "FGK_hosts": qtile(hq[hq.cls == "FGK"]),
        "snr_le5_failures": qtile(q[q.snr <= SNR_MIN]),
        "snr_gt5_passes": qtile(q[q.snr > SNR_MIN]),
        "failures_by_class": {k: int(v) for k, v in
                              hq[hq.snr <= SNR_MIN].cls.value_counts().items()},
    }


# ----------------------------------------------- optional full arm rerun ----
_W = {}


def _init_worker(base, comp, alpha, beta):
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[v] = "1"
    _W.update(base=base, comp=comp, alpha=alpha, beta=beta)


def _one_rep(task):
    rep, seed = task
    rng = np.random.default_rng(seed)
    d = _S15.synth(_W["base"], _W["beta"], _W["alpha"], rng)
    r = _S15.pooled_estimator(d, _W["comp"])
    return rep, (float(r[0]) if r else float("nan"))


def rerun_arms(frozen, reps, workers):
    """Optional pooled-only null/signal arms with the production cut applied.
    Same synth/estimator/calibration/seeds as frozen appendix_C_ruwe.py."""
    global _S15
    _S15 = load_frozen(frozen)
    out = {}
    SEEDS = {("M", "null"): 7201, ("M", "signal"): 7202,
             ("FGK", "null"): 7211, ("FGK", "signal"): 7212}
    for sample in ("M", "FGK"):
        cfg = set_cfg(_S15, sample)
        base, keep, _ = build_cut_frame(_S15, frozen, sample)
        comp = _S15.m9.build_design(keep, "vtan", valley=cfg["valley"],
                                    fe=False, weights=keep.wgt.values)["init_comp"]
        alpha = _S15.calibrate_alpha(keep, cfg["target_frac"])
        out[sample] = {}
        for arm, beta in (("null", 0.0), ("signal", -0.14)):
            children = np.random.SeedSequence(SEEDS[(sample, arm)]).spawn(reps)
            todo = [(i, children[i]) for i in range(reps)]
            t0 = time.time()
            if workers <= 1:
                _init_worker(keep, comp, alpha, beta)
                rows = [_one_rep(t) for t in todo]
            else:
                method = "fork" if "fork" in mp.get_all_start_methods() else "spawn"
                with mp.get_context(method).Pool(workers, initializer=_init_worker,
                                                 initargs=(keep, comp, alpha, beta)) as pool:
                    rows = pool.map(_one_rep, todo, chunksize=8)
            b = np.array([r[1] for r in rows]); b = b[np.isfinite(b)]
            out[sample][arm] = {"n": int(len(b)), "median": float(np.median(b)),
                                "sd": float(np.std(b, ddof=1)),
                                "q025": float(np.quantile(b, 0.025)),
                                "q975": float(np.quantile(b, 0.975))}
            print(f"  {sample}/{arm}: {len(b)} reps in {time.time()-t0:.1f}s", flush=True)
        nul, sig = out[sample]["null"], out[sample]["signal"]
        out[sample]["arm_separation_nullSD"] = float(
            (sig["median"] - nul["median"]) / nul["sd"])
    return out


# ------------------------------------------------------------------- main ---
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--frozen-root", type=pathlib.Path, default=DEFAULT_FROZEN)
    ap.add_argument("--out", type=pathlib.Path, default=DEFAULT_OUT)
    ap.add_argument("--gate-reps", type=int, default=0,
                    help="rerun pooled arms with production cut (0 = carry archive)")
    ap.add_argument("--workers", type=int, default=0)
    args = ap.parse_args()
    frozen, out = args.frozen_root, args.out
    out.mkdir(parents=True, exist_ok=True)

    print("quality table check ...", flush=True)
    qt = quality_table_check(frozen)
    s15 = load_frozen(frozen)
    arch = json.loads((frozen / "results" / "appendix" / "ruwe_clean_summary.json").read_text())
    verify = verify_archive(frozen)
    para = parallax_quantiles(frozen)

    prod = {}
    for sample in ("M", "FGK"):
        cfg = set_cfg(s15, sample)
        base, keep, qs = build_cut_frame(s15, frozen, sample)
        rb_pre = real_beta(s15, base, cfg["valley"])
        rb_post = real_beta(s15, keep, cfg["valley"])
        prod[sample] = {
            "quality_posthoc_recomputed": qs,
            "quality_matches_archive": all(
                qs[k] == arch[sample]["quality"].get(k)
                for k in ("planets_pre", "hosts_pre", "ruwe_gt1.4", "psnr_le5",
                          "planets_post", "hosts_post", "planets_dropped", "hosts_dropped")),
            "baselineARCHIVE": {
                "N": arch[sample]["quality"]["planets_pre"],
                "hosts": arch[sample]["quality"]["hosts_pre"],
                "ESS": 139.50744116694585 if sample == "M" else 1644.987404492584,
                "real_beta": rb_pre[0], "real_se": rb_pre[1]},
            "production_cut_RECOMPUTED": {
                "N": qs["planets_post"], "hosts": qs["hosts_post"],
                "ESS_recomputed": ess(keep.wgt),
                "ESS_retention": ess(keep.wgt) / len(keep),
                "real_beta_recomputed": rb_post[0],
                "real_se_recomputed": rb_post[1]},
            "archive_real_beta_cleaned": arch[sample]["real_beta_cleaned"],
            "archive_real_se_cleaned": arch[sample]["real_se_cleaned"],
            "real_beta_matches_archive": bool(
                abs(rb_post[0] - arch[sample]["real_beta_cleaned"]) < 1e-12),
        }
        print(f"{sample}: N {qs['planets_pre']}->{qs['planets_post']}  "
              f"ESS {prod[sample]['baselineARCHIVE']['ESS']:.2f}->"
              f"{prod[sample]['production_cut_RECOMPUTED']['ESS_recomputed']:.2f}  "
              f"real beta {rb_pre[0]:+.5f}->{rb_post[0]:+.5f}", flush=True)

    arms = None
    if args.gate_reps > 0:
        workers = args.workers or max(1, (os.cpu_count() or 1) - 1)
        print(f"rerunning arms: reps={args.gate_reps} workers={workers}", flush=True)
        arms = rerun_arms(frozen, args.gate_reps, workers)

    summary = {
        "meta": {"rule": RULE, "estimator": "pooled", "mode": "legacy",
                 "frozen_root": str(frozen),
                 "quality_file_used": "data/gaia_quality_DR3_v2.csv (identical to DR3.csv used by archive)",
                 "ESS_definition": "sum(w)^2 / sum(w^2) on 95%-trimmed completeness weights",
                 "arms": ("RECOMPUTED from archived raw ruwe_*_raw.csv" if arms is None
                          else f"RERUN with production cut, reps={args.gate_reps}"),
                 "recomputed_vs_carried": {
                     "RECOMPUTED": ["post-hoc quality counts", "arm statistics/separation from raw CSVs",
                                    "production-cut N/hosts/ESS", "production-cut real beta/SE",
                                    "parallax error quantiles"],
                     "CARRIED_FROM_ARCHIVE": ["baseline (no-cut) ESS values",
                                              "baseline production arm separations (-0.19 M / -1.01 FGK)",
                                              "null/signal arm draws when --gate-reps=0",
                                              "verdict labels (FAIL/inconclusive)"]}},
        "quality_table_check": qt,
        "posthoc_verification": verify,
        "production_cut": prod,
        "parallax": para,
        "optional_arm_rerun": arms,
    }
    path = out / "gaia_cut_summary.json"
    path.write_text(json.dumps(summary, indent=2))
    print(f"wrote {path}", flush=True)


if __name__ == "__main__":
    main()
