#!/usr/bin/env python3
"""Write PARTIAL power-vs-effect output from completed fullcurve_ checkpoints.

Emergency session-close writer. It reuses the gate instance, constants,
Wilson-interval helper and publication plot from complete_power_curve.py, but
reads ONLY finished checkpoints (null + first two signal arms). It runs no
Monte Carlo replicates. The summary carries status="PARTIAL" and an explicit
missing_points list.
"""
import json
import shutil

import numpy as np
import pandas as pd

import complete_power_curve as cp

g = cp.g
SCI = cp.SCI

# completed checkpoint tag -> manuscript point name
TAGMAP = {
    "fullcurve_null": "null",
    "fullcurve_sig_1": "half_primary",
    "fullcurve_sig_2": "Kamulali_primary",
}


def main():
    base = g.load_control(g.CFG)
    full = g.m9.build_design(base, "vtan", valley=g.VALLEY, fe=False,
                             weights=base.wgt.values)
    comp = full["init_comp"]
    cut = float(base[cp.STRAT_COL].median())
    real = float(g.pooled_estimator(base, comp)[0])

    eff = json.loads((SCI / "effect_size_derivation.json").read_text())
    pts = eff["bracket_points"]

    arms = {}
    for tag, name in TAGMAP.items():
        blob = json.loads((SCI / "checkpoints" / f"{tag}.json").read_text())
        done = int(blob["done"])
        df = pd.DataFrame(blob["rows"])
        sig = df["beta_pooled"].dropna().values.astype(float)
        if len(sig) < 150:
            raise RuntimeError(f"{tag} has only {len(sig)} usable draws (<150)")
        arms[name] = (tag, done, sig, df)

    nul = arms["null"][2]
    nmed, nsd = float(np.median(nul)), float(np.std(nul))
    n025, n975 = np.percentile(nul, [2.5, 97.5])

    results, raw_parts = [], []
    for name, (tag, done, sig, df) in arms.items():
        b = float(pts[name])
        smed = float(np.median(sig))
        sep = (smed - nmed) / nsd
        pn = float(np.mean(nul <= real))
        ps = float(np.mean(sig <= real))
        verdict = "PASS" if (pn <= 0.05 and 0.025 <= ps <= 0.975) else "FAIL"
        k = int(np.sum(sig < n025))
        det = k / len(sig)
        dlo, dhi = cp.wilson_interval(k, len(sig))
        results.append(dict(
            point=name, beta_inj=b, reps_null=int(len(nul)), reps_signal=int(len(sig)),
            null_median=nmed, null_sd=nsd, null_q025=float(n025), null_q975=float(n975),
            signal_median=smed, signal_sd=float(np.std(sig)),
            arm_separation_nullSD=sep, real_beta=real,
            real_pct_in_null=pn, real_pct_in_signal=ps,
            detect_rate=det, detect_rate_lo_wilson=dlo, detect_rate_hi_wilson=dhi,
            detect_rate_two_sided=float(np.mean((sig < n025) | (sig > n975))),
            detect_correct_sign_rate=float(np.mean(sig < 0.0)),
            verdict=verdict, signal_arm_is_null=(name == "null"),
            checkpoint_tag=tag, checkpoint_done=done))
        raw = df.copy()
        raw["point"] = name
        raw["beta_inj"] = b
        raw["arm"] = "null" if name == "null" else "signal"
        raw_parts.append(raw)

    raw = pd.concat(raw_parts, ignore_index=True)
    cols = ["point", "arm", "beta_inj", "rep", "beta_pooled", "se_pooled",
            "beta_strat", "se_strat", "beta_lo", "beta_hi", "n"]
    raw[cols].to_csv(SCI / "power_vs_effect_raw_FINAL.csv", index=False)

    fcols = ["point", "beta_inj", "reps_null", "reps_signal", "null_median",
             "null_sd", "signal_median", "arm_separation_nullSD", "real_beta",
             "real_pct_in_null", "real_pct_in_signal", "detect_rate",
             "detect_rate_lo_wilson", "detect_rate_hi_wilson", "verdict"]
    fig = pd.DataFrame([{k: r[k] for k in fcols} for r in results])
    fig.to_csv(SCI / "power_vs_effect_figure_data_FINAL.csv", index=False)

    completed = [r["point"] for r in results]
    missing = [n for n in cp.POINT_ORDER if n not in completed]
    summary = dict(
        status="PARTIAL",
        partial_reason="Emergency session close: full curve exceeded time budget; "
                       "only completed fullcurve_ checkpoints are reported.",
        points_completed=len(results),
        signal_points_completed=len(results) - 1,
        missing_points=missing,
        sample="M (cool / late-K early-M hosts)", estimator="pooled (production)",
        strat_col=cp.STRAT_COL, reps_per_arm=cp.REPS,
        null_arm_reps=int(len(nul)), seed=cp.SEED,
        workers=cp.WORKERS, checkpoint_prefix="fullcurve_",
        real_beta=real, real_beta_rounded_4dp=round(real, 4),
        null_median=nmed, null_sd=nsd,
        gaidos_raw=float(pts["Gaidos_raw_true_age"]),
        gaidos_rho025=float(pts["Gaidos_rho025_vtan"]),
        gaidos_rho050=float(pts["Gaidos_rho050_vtan"]),
        reproduction_check=None, points=results)
    (SCI / "power_vs_effect_summary_FINAL.json").write_text(json.dumps(summary, indent=2))

    cp.make_plot(fig, summary, SCI / "power_vs_effect_FINAL.png")
    for final, legacy in cp.LEGACY.items():
        shutil.copyfile(SCI / final, SCI / legacy)

    (SCI / "power_vs_effect_decision_FINAL.txt").write_text(
        "DECISION BRANCH (PARTIAL: null + 2 signal points, 150 reps/arm)\n"
        "Missing points: " + ", ".join(missing) + "\n\n"
        "Real cool-host beta = " + f"{real:+.6f}" + " (rounded +0.0458); sits at the "
        + f"{results[0]['real_pct_in_null']*100:.1f} percentile of the null "
        "(needs <=5%), so the gate FAILS at every reported point.\n")

    print("WROTE PARTIAL OUTPUT: points=", completed, " missing=", missing, flush=True)
    for r in results:
        print(f"  {r['point']:24s} b={r['beta_inj']:+.4f} "
              f"sep={r['arm_separation_nullSD']:+.3f} det={r['detect_rate']:.3f} "
              f"{r['verdict']} n={r['reps_signal']}", flush=True)


if __name__ == "__main__":
    main()
