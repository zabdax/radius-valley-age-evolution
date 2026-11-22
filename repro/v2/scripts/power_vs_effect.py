#!/usr/bin/env python3
"""Power-vs-effect curve for the calibrated cool-host decision gate.

Writes all deliverables to manuscripts/q1-paper/v2-2026-11/results/science/.
Reuses the frozen production gate in radius-valley-age-evolution/code/
15_stratified_gate.py (fast mirrored synth, same RNG stream as 13_calibration_fgk.py).
"""
import functools
import importlib.util
import json
import math
import os
import pathlib
import sys
import time

import numpy as np
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[4]
FROZEN = ROOT / "radius-valley-age-evolution"
SCI = (ROOT / "manuscripts" / "q1-paper" / "v2-2026-11" / "results" / "science")
SCI.mkdir(parents=True, exist_ok=True)

# ---- load frozen gate machinery (read only) ----
_spec = importlib.util.spec_from_file_location("gate", FROZEN / "code" / "15_stratified_gate.py")
g = importlib.util.module_from_spec(_spec)
sys.modules["gate"] = g          # needed for spawn pickling of g._init_worker
_spec.loader.exec_module(g)
g.CFG = g.SAMPLES["M"]                      # cool / M-dwarf host sample
g.TARGET_FRAC = g.CFG["target_frac"]
g.VALLEY = g.CFG["valley"]
g.GAMMA_LOGP = g.CFG["gamma_logp"]
g.OUT = SCI
g.CKPT = SCI / "checkpoints"
g.CKPT.mkdir(parents=True, exist_ok=True)

REPS = 300          # replicates per arm (minimum 150 required)
SEED = 20260826
STRAT_COL = "sy_vmag"
CUT = None
NORM_SEP = 1.60     # group-mean separation for a median split (SD units)


def logit(p):
    return math.log(p / (1.0 - p))


def beta_from_se_sn(r_young, r_old, sep=NORM_SEP):
    fy, fo = 1.0 / (1.0 + r_young), 1.0 / (1.0 + r_old)
    d = logit(fo) - logit(fy)
    return dict(f_young=fy, f_old=fo, logit_young=logit(fy), logit_old=logit(fo),
                delta_logit=d, beta_per_sd=d / sep)


def beta_from_sn_se(r_young, r_old, sep=NORM_SEP):
    fy, fo = r_young / (1.0 + r_young), r_old / (1.0 + r_old)
    d = logit(fo) - logit(fy)
    return dict(f_young=fy, f_old=fo, logit_young=logit(fy), logit_old=logit(fo),
                delta_logit=d, beta_per_sd=d / sep)


def derive_effect_sizes():
    d = {}
    d["formulas"] = {
        "SE_SN_ratio": "f_SN = 1/(1 + SE/SN); logit gap = logit(f_old) - logit(f_young)",
        "SN_SE_ratio": "f_SN = (SN/SE)/(1 + SN/SE); logit gap as above",
        "beta_per_SD": "beta = delta_logit / group_mean_separation; separation = 1.60 SD for a median split",
        "proxy_attenuation": "beta_in_vtan = rho * beta_true; cool-host rho ~= 0.25 (not significant, p~0.43)",
    }
    d["median_split_separation_SD"] = NORM_SEP
    d["sources"] = {}

    # Kamulali+26 -- primary production value
    k = beta_from_se_sn(0.51, 0.64)
    k.update(published_ratio="SE/SN 0.51 (<3 Gyr) -> 0.64 (>3 Gyr)",
             source_note="research/notes/astronomy-astrophysics-manuscript-no-aanda.md",
             split="3 Gyr (not a median split; 1.60 divisor is an approximation)",
             role="PRIMARY injected effect used in production (rounded to -0.14/SD)")
    d["sources"]["Kamulali+26"] = k

    # Berger+20
    b = beta_from_se_sn(0.61, 1.00)
    b.update(published_ratio="SE/SN 0.61 (<1 Gyr) -> 1.00 (>1 Gyr)",
             source_note="research/notes/arxiv190505900v1-astro-phep-15-may-2019.md; ratio quoted in research/notes/draft-version-february-23-2021.md line 25-26",
             split="1 Gyr (not a median split; 1.60 divisor approximate)",
             role="FGK literature anchor, larger than primary")
    d["sources"]["Berger+20"] = b

    # Sandoval+21 -- exact ratios ARE available in the frozen note
    s3 = beta_from_se_sn(0.82, 0.97)   # base sample, 3 Gyr split (lines 404-405)
    s1 = beta_from_se_sn(0.76, 0.94)   # base sample, 1 Gyr split (lines 398-399)
    sh = beta_from_se_sn(0.78, 1.00)   # Teff>5500 K, 3 Gyr split (lines 412-414)
    d["sources"]["Sandoval+21"] = {
        "base_3Gyr": dict(s3, published_ratio="SE/SN 0.82 (<3 Gyr) +/- 0.10 -> 0.97 (>3 Gyr) +/- 0.04",
                          role="FGK anchor (exact ratio found; no fallback needed)"),
        "base_1Gyr": dict(s1, published_ratio="SE/SN 0.76 (<1 Gyr) +/- 0.08 -> 0.94 (>1 Gyr) +/- 0.07"),
        "hot_3Gyr": dict(sh, published_ratio="SE/SN 0.78 (<3 Gyr) +/- 0.17 -> 1.00 (>3 Gyr) +/- 0.04, Teff>5500 K"),
        "source_note": "research/notes/draft-version-february-23-2021.md lines 398-414",
        "note": "trend ~2 sigma over 1-10 Gyr; ratio available, so Kamulali fallback not required",
    }

    # Gaidos+24
    gd = beta_from_sn_se(2.38, 0.58)
    gd.update(published_ratio="SN/SE 2.38 (mean age 2.6 Gyr) -> 0.58 (mean age 5.4 Gyr)",
              source_note="research/notes/mnras-000-114.md lines 805-816 (median gyro age split 3.86 Gyr, KS p=0.08)",
              split="3.86 Gyr median split -- exact median split, 1.60 divisor exact",
              role="Gaidos-sized cool-host signal",
              beta_true_age_per_SD=gd["beta_per_sd"],
              caveat="rho=0.25 cool-host proxy correlation is not significant (p~0.43); raw -0.88 is an optimistic upper bracket")
    gd["beta_vtan_rho0_25"] = 0.25 * gd["beta_per_sd"]
    gd["beta_vtan_rho0_50"] = 0.50 * gd["beta_per_sd"]
    d["sources"]["Gaidos+24"] = gd

    d["bracket_points"] = {
        "null": 0.0,
        "half_primary": -0.07,
        "Kamulali_primary": round(k["beta_per_sd"], 4),
        "Gaidos_rho025_vtan": round(0.25 * gd["beta_per_sd"], 4),
        "Berger_FGK": round(b["beta_per_sd"], 4),
        "Sandoval_FGK_base3Gyr": round(s3["beta_per_sd"], 4),
        "Gaidos_rho050_vtan": round(0.50 * gd["beta_per_sd"], 4),
        "strong_2x": -0.60,
        "Gaidos_raw_true_age": round(gd["beta_per_sd"], 4),
        "upper": -1.20,
    }
    d["rho_cool_host"] = 0.25
    d["cross_check"] = "manuscripts/q1-paper/v2-2026-11/results/science/effect_size_derivation_orchestrator_check.md"
    (SCI / "effect_size_derivation.json").write_text(json.dumps(d, indent=2))
    return d


def run_curve(effect):
    base = g.load_control(g.CFG)
    full = g.m9.build_design(base, "vtan", valley=g.VALLEY, fe=False, weights=base.wgt.values)
    comp = full["init_comp"]
    cut = float(base[STRAT_COL].median())
    real = g.pooled_estimator(base, comp)[0]
    print(f"control N={len(base)} ESS={base.wgt.sum()**2/(base.wgt**2).sum():.0f} "
          f"real_beta={real:+.4f} cut={cut}", flush=True)

    pts = effect["bracket_points"]
    # order: null first, then signal points sorted by magnitude
    names = ["null"] + [n for n in pts if n != "null"]
    results = []
    t_start = time.time()

    null_df = g.run_arm(base, comp, STRAT_COL, cut, 0.0, REPS, SEED,
                        "curve_null", workers=7)
    nul = null_df["beta_pooled"].dropna().values
    print(f"null arm done ({len(nul)} draws) in {time.time()-t_start:.0f}s", flush=True)
    nmed, nsd = float(np.median(nul)), float(np.std(nul))
    n025, n975 = np.percentile(nul, [2.5, 97.5])

    for i, name in enumerate(names[1:]):
        b = float(pts[name])
        tag = f"curve_sig_{i}"
        df = g.run_arm(base, comp, STRAT_COL, cut, b, REPS,
                       SEED + 1000 + i, tag, workers=7)
        sig = df["beta_pooled"].dropna().values
        smed = float(np.median(sig))
        sep = (smed - nmed) / nsd
        pn = float(np.mean(nul <= real))
        ps = float(np.mean(sig <= real))
        verdict = "PASS" if (pn <= 0.05 and 0.025 <= ps <= 0.975) else "FAIL"
        # detection rate: signal draws distinguishable from null (outside null 95%)
        det_two = float(np.mean((sig < n025) | (sig > n975)))
        det_one = float(np.mean(sig < n025))  # negative-effect tail
        det_sign = float(np.mean(sig < 0))
        print(f"  {name:24s} b={b:+.4f} sep={sep:+.3f} pct_null={pn:.3f} "
              f"pct_sig={ps:.3f} det={det_one:.3f} {verdict} "
              f"({len(sig)} draws)", flush=True)
        results.append(dict(point=name, beta_inj=b, reps_null=len(nul), reps_signal=len(sig),
                            null_median=nmed, null_sd=nsd, null_q025=float(n025), null_q975=float(n975),
                            signal_median=smed, signal_sd=float(np.std(sig)),
                            arm_separation_nullSD=sep, real_beta=real,
                            real_pct_in_null=pn, real_pct_in_signal=ps,
                            detect_rate_one_sided=det_one, detect_rate_two_sided=det_two,
                            detect_correct_sign_rate=det_sign, verdict=verdict))
        # raw rows for this point
        raw = df.copy()
        raw["point"] = name
        raw["beta_inj"] = b
        raw["arm"] = "signal"
        raw.to_csv(SCI / f"_tmp_raw_{i}.csv", index=False)

    # raw file: null draws + all signal draws
    null_raw = null_df.copy()
    null_raw["point"] = "null"
    null_raw["beta_inj"] = 0.0
    null_raw["arm"] = "null"
    parts = [null_raw]
    for i in range(len(names) - 1):
        parts.append(pd.read_csv(SCI / f"_tmp_raw_{i}.csv"))
        (SCI / f"_tmp_raw_{i}.csv").unlink()
    raw = pd.concat(parts, ignore_index=True)
    cols = ["point", "arm", "beta_inj", "rep", "beta_pooled", "se_pooled",
            "beta_strat", "se_strat", "beta_lo", "beta_hi", "n"]
    raw[cols].to_csv(SCI / "power_vs_effect_raw.csv", index=False)

    # figure data
    fig = pd.DataFrame([{k: r[k] for k in ("point", "beta_inj", "arm_separation_nullSD",
                        "real_pct_in_null", "real_pct_in_signal", "detect_rate_one_sided",
                        "detect_rate_two_sided", "verdict", "reps_signal")} for r in results])
    fig.to_csv(SCI / "power_vs_effect_figure_data.csv", index=False)

    summary = dict(sample="M (cool / late-K early-M hosts)", estimator="pooled (production)",
                   strat_col=STRAT_COL, reps_per_arm=REPS, seed=SEED,
                   runtime_s=round(time.time() - t_start, 1),
                   real_beta=real, null_median=nmed, null_sd=nsd,
                   gaidos_raw=pts["Gaidos_raw_true_age"],
                   gaidos_rho025=pts["Gaidos_rho025_vtan"],
                   gaidos_rho050=pts["Gaidos_rho050_vtan"],
                   reproduction_check=None, points=results)
    return summary, raw, fig


def archived_reproduction():
    """Exact reproduction of the production -0.19 point from archived 1000-rep draws."""
    p = FROZEN / "results" / "stratified_gate" / "M"
    nul = pd.read_csv(p / "strat_gate_M_vmag_null_r1000_raw.csv")["beta_pooled"].dropna().values
    sig = pd.read_csv(p / "strat_gate_M_vmag_signal_r1000_raw.csv")["beta_pooled"].dropna().values
    sep = (np.median(sig) - np.median(nul)) / np.std(nul)
    return dict(source="archived stratified_gate/M/strat_gate_M_vmag_*_r1000_raw.csv",
                beta_inj=-0.14, null_median=float(np.median(nul)), null_sd=float(np.std(nul)),
                signal_median=float(np.median(sig)), arm_separation_nullSD=float(sep),
                n_null=len(nul), n_signal=len(sig), verdict="FAIL" if abs(sep) < 0.05 else "CHECK")


def make_plot(fig, summary):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    x = fig["beta_inj"].values
    y = fig["arm_separation_nullSD"].values
    det = fig["detect_rate_one_sided"].values
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    ax.plot(x, y, "o-", color="#1f77b4", lw=1.8, ms=6, label="Arm separation (null SD)")
    ax.axhline(0, color="0.6", lw=0.8)
    ax.axhline(-1.96, color="0.4", ls=":", lw=1, label="95% null detection edge")
    for b, lab in [(summary["gaidos_rho025"], "Gaidos rho=0.25"),
                   (summary["gaidos_raw"], "Gaidos raw")]:
        ax.axvline(b, color="#d62728", ls="--", lw=1, alpha=0.8)
        ax.annotate(lab, (b, ax.get_ylim()[0]), xytext=(4, 4), textcoords="offset points",
                    fontsize=7, color="#d62728", rotation=90, va="bottom")
    ax.set_xlabel("Injected effect, beta per SD of v_tan (logit SN fraction)")
    ax.set_ylabel("Signal minus null arm separation (null SD)")
    ax.set_title("Cool-host gate power vs injected effect")
    ax.grid(alpha=0.25)
    ax1 = ax.twinx()
    ax1.plot(x, det, "s--", color="#ff7f0e", lw=1.4, ms=5, alpha=0.8,
             label="Detection rate (one-sided)")
    ax1.set_ylabel("Detection rate (fraction of signal draws)")
    ax1.set_ylim(-0.05, 1.05)
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax1.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="lower left", fontsize=7)
    fig.tight_layout()
    fig.savefig(SCI / "power_vs_effect.png", dpi=200)
    plt.close(fig)


def main():
    t0 = time.time()
    effect = derive_effect_sizes()
    summary, raw, fig = run_curve(effect)
    rep = archived_reproduction()
    summary["reproduction_check"] = rep
    (SCI / "power_vs_effect_summary.json").write_text(json.dumps(summary, indent=2))
    make_plot(fig, summary)
    # decision branch text
    gsep = [r for r in summary["points"] if r["point"] == "Gaidos_raw_true_age"][0]
    gar = [r for r in summary["points"] if r["point"] == "Gaidos_rho025_vtan"][0]
    verdict = ("FAIL at Gaidos-sized effect: rotation-based cool-host signals cannot be "
               "reproduced with kinematic proxies in this archive"
               if gsep["verdict"] == "FAIL" and gar["verdict"] == "FAIL"
               else "PASSES at Gaidos-sized effect: calibrated adjudication required")
    (SCI / "power_vs_effect_decision.txt").write_text(
        "DECISION BRANCH\n" + verdict + "\n\n"
        "Gaidos raw (-0.88/SD true age): " + json.dumps(gsep, indent=2) + "\n"
        "Gaidos rho=0.25 (-0.22/SD v_tan): " + json.dumps(gar, indent=2) + "\n"
        "Archived reproduction: " + json.dumps(rep, indent=2) + "\n"
        "Alternatives engaged: (1) injected effect too small relative to true rotation-age signal; "
        "(2) composition effect near the water-condensation boundary (Luque & Palle 2022), which would "
        "produce a real radius-population shift that a kinematic proxy cannot recover.\n")
    print("TOTAL RUNTIME", round(time.time() - t0, 1), "s")
    print("VERDICT", verdict)
    print("REPRO", rep)


if __name__ == "__main__":
    main()
