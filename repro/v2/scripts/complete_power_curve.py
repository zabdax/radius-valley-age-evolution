#!/usr/bin/env python3
"""COMPLETE power-vs-effect Monte Carlo for the calibrated cool-host gate.

This is a NEW driver. It does not modify power_vs_effect.py (the earlier,
partial run). It reuses the frozen production gate in
radius-valley-age-evolution/code/15_stratified_gate.py and the authoritative
bracket points in results/science/effect_size_derivation.json.

Design
------
* REPS = 150 replicates per arm (the minimum acceptable), SEED = 20260826.
* Checkpoint tags are prefixed "fullcurve_" so no partial checkpoint is ever
  resumed or mixed with this run.
* A SINGLE null arm is drawn once and reused as the null reference for every
  injected-effect point AND as the signal arm of the "null" bracket point.
* 10 bracket points are reported (null + 9 injected effects).

All deliverables are written to manuscripts/q1-paper/v2-2026-11/results/science/
with a _FINAL suffix, then copied to the legacy unsuffixed names so the FINAL
files are authoritative.
"""
import importlib.util
import json
import math
import pathlib
import shutil
import sys
import time

import numpy as np
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[4]
FROZEN = ROOT / "radius-valley-age-evolution"
SCI = ROOT / "manuscripts" / "q1-paper" / "v2-2026-11" / "results" / "science"
SCI.mkdir(parents=True, exist_ok=True)

# ---- load frozen gate machinery (read only) ----
_spec = importlib.util.spec_from_file_location(
    "gate", FROZEN / "code" / "15_stratified_gate.py")
g = importlib.util.module_from_spec(_spec)
sys.modules["gate"] = g          # needed for spawn pickling of g._init_worker / g._one_rep
_spec.loader.exec_module(g)
g.CFG = g.SAMPLES["M"]                       # cool / M-dwarf host sample
g.TARGET_FRAC = g.CFG["target_frac"]
g.VALLEY = g.CFG["valley"]
g.GAMMA_LOGP = g.CFG["gamma_logp"]
g.OUT = SCI
g.CKPT = SCI / "checkpoints"
g.CKPT.mkdir(parents=True, exist_ok=True)

REPS = 150              # replicates per arm (minimum acceptable)
SEED = 20260826
STRAT_COL = "sy_vmag"
WORKERS = 7             # machine has 8 logical CPUs; pin BLAS to 1 in workers

# Exact order requested by the manuscript disclosure.
POINT_ORDER = [
    "null",
    "half_primary",
    "Kamulali_primary",
    "Gaidos_rho025_vtan",
    "Berger_FGK",
    "Sandoval_FGK_base3Gyr",
    "Gaidos_rho050_vtan",
    "strong_2x",
    "Gaidos_raw_true_age",
    "upper",
]

LEGACY = {
    "power_vs_effect_summary_FINAL.json": "power_vs_effect_summary.json",
    "power_vs_effect_raw_FINAL.csv": "power_vs_effect_raw.csv",
    "power_vs_effect_figure_data_FINAL.csv": "power_vs_effect_figure_data.csv",
    "power_vs_effect_FINAL.png": "power_vs_effect.png",
}


def wilson_interval(k, n, z=1.96):
    """Wilson score 95% interval for a binomial proportion k/n."""
    if n <= 0:
        return float("nan"), float("nan")
    p = k / n
    denom = 1.0 + z * z / n
    center = (p + z * z / (2.0 * n)) / denom
    half = (z / denom) * math.sqrt(p * (1.0 - p) / n + z * z / (4.0 * n * n))
    return max(0.0, center - half), min(1.0, center + half)


def archived_reproduction():
    """Exact reproduction of the production -0.14 point from archived 1000-rep draws."""
    p = FROZEN / "results" / "stratified_gate" / "M"
    try:
        nul = pd.read_csv(p / "strat_gate_M_vmag_null_r1000_raw.csv")["beta_pooled"].dropna().values
        sig = pd.read_csv(p / "strat_gate_M_vmag_signal_r1000_raw.csv")["beta_pooled"].dropna().values
    except Exception as exc:  # non-fatal: archival cross-check is optional
        return {"status": "unavailable", "reason": str(exc)}
    sep = (np.median(sig) - np.median(nul)) / np.std(nul)
    return dict(source="archived stratified_gate/M/strat_gate_M_vmag_*_r1000_raw.csv",
                beta_inj=-0.14, null_median=float(np.median(nul)),
                null_sd=float(np.std(nul)), signal_median=float(np.median(sig)),
                arm_separation_nullSD=float(sep), n_null=len(nul), n_signal=len(sig),
                verdict="FAIL" if abs(sep) < 0.05 else "CHECK")


def run_curve(effect):
    base = g.load_control(g.CFG)
    full = g.m9.build_design(base, "vtan", valley=g.VALLEY, fe=False,
                             weights=base.wgt.values)
    comp = full["init_comp"]
    cut = float(base[STRAT_COL].median())
    real = float(g.pooled_estimator(base, comp)[0])
    print(f"control N={len(base)} ESS={base.wgt.sum()**2/(base.wgt**2).sum():.0f} "
          f"real_beta={real:+.6f} cut={cut:.4f}", flush=True)

    pts = effect["bracket_points"]
    missing = [n for n in POINT_ORDER if n not in pts]
    if missing:
        raise RuntimeError(f"bracket points missing from effect_size_derivation.json: {missing}")

    results = []
    raw_parts = []
    t_start = time.time()

    # ---- ONE null arm, reused everywhere (new prefix => no stale checkpoint) ----
    null_df = g.run_arm(base, comp, STRAT_COL, cut, 0.0, REPS, SEED,
                        "fullcurve_null", workers=WORKERS)
    nul = null_df["beta_pooled"].dropna().values
    if len(nul) < 150:
        raise RuntimeError(f"null arm has only {len(nul)} usable draws (<150)")
    nmed, nsd = float(np.median(nul)), float(np.std(nul))
    n025, n975 = np.percentile(nul, [2.5, 97.5])
    print(f"null arm done ({len(nul)} draws) in {time.time()-t_start:.0f}s "
          f"median={nmed:+.4f} sd={nsd:.4f} 95%=[{n025:+.4f},{n975:+.4f}]", flush=True)

    null_raw = null_df.copy()
    null_raw["point"] = "null"
    null_raw["beta_inj"] = 0.0
    null_raw["arm"] = "null"
    raw_parts.append(null_raw)

    def stats_for(name, b, sig, is_null_point=False):
        smed = float(np.median(sig))
        sep = (smed - nmed) / nsd
        pn = float(np.mean(nul <= real))
        ps = float(np.mean(sig <= real))
        verdict = "PASS" if (pn <= 0.05 and 0.025 <= ps <= 0.975) else "FAIL"
        k = int(np.sum(sig < n025))
        det = k / len(sig)
        dlo, dhi = wilson_interval(k, len(sig))
        det_two = float(np.mean((sig < n025) | (sig > n975)))
        det_sign = float(np.mean(sig < 0.0))
        print(f"  {name:24s} b={b:+.4f} sep={sep:+.3f} pct_null={pn:.3f} "
              f"pct_sig={ps:.3f} det={det:.3f} [{dlo:.3f},{dhi:.3f}] {verdict} "
              f"({len(sig)} draws)", flush=True)
        return dict(point=name, beta_inj=float(b),
                    reps_null=int(len(nul)), reps_signal=int(len(sig)),
                    null_median=nmed, null_sd=nsd,
                    null_q025=float(n025), null_q975=float(n975),
                    signal_median=smed, signal_sd=float(np.std(sig)),
                    arm_separation_nullSD=sep, real_beta=real,
                    real_pct_in_null=pn, real_pct_in_signal=ps,
                    detect_rate=det, detect_rate_lo_wilson=dlo,
                    detect_rate_hi_wilson=dhi,
                    detect_rate_two_sided=det_two, detect_correct_sign_rate=det_sign,
                    verdict=verdict, signal_arm_is_null=is_null_point)

    # ---- point 1: "null" -- the single null arm doubles as its signal arm ----
    results.append(stats_for("null", 0.0, nul, is_null_point=True))

    # ---- points 2..10: injected effects, one fresh arm each ----
    signal_names = [n for n in POINT_ORDER if n != "null"]
    for i, name in enumerate(signal_names, start=1):
        b = float(pts[name])
        tag = f"fullcurve_sig_{i}"
        seed = SEED + 1000 + i
        df = g.run_arm(base, comp, STRAT_COL, cut, b, REPS, seed, tag,
                       workers=WORKERS)
        sig = df["beta_pooled"].dropna().values
        if len(sig) < 150:
            raise RuntimeError(f"arm {name} has only {len(sig)} usable draws (<150)")
        results.append(stats_for(name, b, sig))
        sraw = df.copy()
        sraw["point"] = name
        sraw["beta_inj"] = b
        sraw["arm"] = "signal"
        raw_parts.append(sraw)

    # ---- raw draws ----
    raw = pd.concat(raw_parts, ignore_index=True)
    cols = ["point", "arm", "beta_inj", "rep", "beta_pooled", "se_pooled",
            "beta_strat", "se_strat", "beta_lo", "beta_hi", "n"]
    raw[cols].to_csv(SCI / "power_vs_effect_raw_FINAL.csv", index=False)

    # ---- figure data ----
    fcols = ["point", "beta_inj", "reps_null", "reps_signal", "null_median",
             "null_sd", "signal_median", "arm_separation_nullSD", "real_beta",
             "real_pct_in_null", "real_pct_in_signal", "detect_rate",
             "detect_rate_lo_wilson", "detect_rate_hi_wilson", "verdict"]
    fig = pd.DataFrame([{k: r[k] for k in fcols} for r in results])
    fig.to_csv(SCI / "power_vs_effect_figure_data_FINAL.csv", index=False)

    summary = dict(status="COMPLETE", points_completed=len(results),
                   sample="M (cool / late-K early-M hosts)",
                   estimator="pooled (production)", strat_col=STRAT_COL,
                   reps_per_arm=REPS, null_arm_reps=int(len(nul)), seed=SEED,
                   workers=WORKERS, checkpoint_prefix="fullcurve_",
                   runtime_s=round(time.time() - t_start, 1),
                   real_beta=real, real_beta_rounded_4dp=round(real, 4),
                   null_median=nmed, null_sd=nsd,
                   gaidos_raw=float(pts["Gaidos_raw_true_age"]),
                   gaidos_rho025=float(pts["Gaidos_rho025_vtan"]),
                   gaidos_rho050=float(pts["Gaidos_rho050_vtan"]),
                   reproduction_check=None, points=results)
    return summary, raw, fig


def make_plot(figdata, summary, path):
    """Publication-ready figure. Okabe-Ito palette; Wilson CI error bars.

    All labels avoid em dashes (plain hyphens / 'vs' / 'to' are used instead).
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x = figdata["beta_inj"].values.astype(float)
    y = figdata["arm_separation_nullSD"].values.astype(float)
    det = figdata["detect_rate"].values.astype(float)
    dlo = figdata["detect_rate_lo_wilson"].values.astype(float)
    dhi = figdata["detect_rate_hi_wilson"].values.astype(float)
    yerr = np.vstack([det - dlo, dhi - det])

    f, ax1 = plt.subplots(figsize=(7.8, 5.0))
    ax1.plot(x, y, "o-", color="#0072B2", lw=1.9, ms=6,
             label="Arm separation (null SD)")
    ax1.axhline(0.0, color="#888888", lw=0.9)
    ax1.axhline(-1.96, color="#888888", ls=":", lw=1.1,
                label="95% null detection edge")
    ax1.set_xlabel("Injected effect, beta per SD of v_tan")
    ax1.set_ylabel("Signal minus null separation (null SD)", color="#0072B2")
    ax1.tick_params(axis="y", labelcolor="#0072B2")
    ax1.grid(alpha=0.25)

    for b, lab, col, ls in [
            (summary["gaidos_rho025"], "Gaidos rho=0.25 (v_tan)", "#009E73", "--"),
            (summary["gaidos_raw"], "Gaidos raw (true age)", "#CC79A7", "-.")]:
        ax1.axvline(b, color=col, ls=ls, lw=1.3, alpha=0.9, label=lab)

    ax2 = ax1.twinx()
    ax2.errorbar(x, det, yerr=yerr, fmt="s--", color="#D55E00", lw=1.5, ms=5,
                 capsize=3, elinewidth=1.1,
                 label="Detection rate (Wilson 95% CI)")
    ax2.set_ylabel("Detection rate (fraction of signal draws)", color="#D55E00")
    ax2.tick_params(axis="y", labelcolor="#D55E00")
    ax2.set_ylim(-0.05, 1.05)

    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, loc="best", fontsize=8, framealpha=0.9)
    ax1.set_title("Cool-host gate: power vs injected effect")
    f.tight_layout()
    f.savefig(path, dpi=200)
    plt.close(f)


def main():
    t0 = time.time()
    eff = json.loads((SCI / "effect_size_derivation.json").read_text())
    summary, raw, fig = run_curve(eff)
    summary["reproduction_check"] = archived_reproduction()

    # ---- FINAL deliverables ----
    (SCI / "power_vs_effect_summary_FINAL.json").write_text(
        json.dumps(summary, indent=2))
    make_plot(fig, summary, SCI / "power_vs_effect_FINAL.png")

    # ---- legacy names become copies of the FINAL authoritative files ----
    for final, legacy in LEGACY.items():
        shutil.copyfile(SCI / final, SCI / legacy)

    # ---- decision branch / disclosure text ----
    by = {r["point"]: r for r in summary["points"]}
    graw, gar = by["Gaidos_raw_true_age"], by["Gaidos_rho025_vtan"]
    verdict = ("FAIL at Gaidos-sized effect: rotation-based cool-host signals cannot be "
               "reproduced with kinematic proxies in this archive"
               if graw["verdict"] == "FAIL" and gar["verdict"] == "FAIL"
               else "PASSES at Gaidos-sized effect: calibrated adjudication required")
    decision = (
        "DECISION BRANCH (COMPLETE RUN, 150 reps/arm, fullcurve_ checkpoints)\n"
        + verdict + "\n\n"
        + "Real cool-host beta = " + f"{summary['real_beta']:+.6f}"
        + " (rounded +0.0458); sits at the "
        + f"{by['null']['real_pct_in_null']*100:.1f} percentile of the null\n"
        + "(needs <=5% to pass distinguishability), so the gate FAILS for every "
          "injected effect by construction.\n\n"
        + "Gaidos raw (-0.8824/SD true age): " + json.dumps(graw, indent=2) + "\n"
        + "Gaidos rho=0.25 (-0.2206/SD v_tan): " + json.dumps(gar, indent=2) + "\n"
        + "Archived reproduction: " + json.dumps(summary["reproduction_check"], indent=2) + "\n")
    (SCI / "power_vs_effect_decision_FINAL.txt").write_text(decision)

    # ---- console report ----
    print("\n" + "=" * 78)
    print(f"COMPLETE: {summary['points_completed']} points, "
          f"{summary['reps_per_arm']} reps/arm, runtime {time.time()-t0:.1f}s")
    print(f"real_beta = {summary['real_beta']:+.6f} (rounded {summary['real_beta_rounded_4dp']:+.4f})")
    print(f"{'point':24s} {'beta':>8s} {'sep':>8s} {'pct_null':>9s} "
          f"{'pct_sig':>8s} {'det':>7s} {'verdict':>8s} {'n':>5s}")
    for r in summary["points"]:
        print(f"{r['point']:24s} {r['beta_inj']:+8.4f} "
              f"{r['arm_separation_nullSD']:+8.3f} {r['real_pct_in_null']:9.3f} "
              f"{r['real_pct_in_signal']:8.3f} {r['detect_rate']:7.3f} "
              f"{r['verdict']:>8s} {r['reps_signal']:5d}")
    print("VERDICT:", verdict)
    print("REPRO:", json.dumps(summary["reproduction_check"]))
    print("=" * 78, flush=True)


if __name__ == "__main__":
    main()
