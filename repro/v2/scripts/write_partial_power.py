#!/usr/bin/env python3
"""Write partial power-vs-effect outputs from completed checkpoints only."""
import json, pathlib, numpy as np, pandas as pd
SCI = pathlib.Path(__file__).resolve().parents[1] / "results" / "science"
CK = SCI / "checkpoints"

NAMES = ["null", "half_primary", "Kamulali_primary", "Gaidos_rho025_vtan",
         "Berger_FGK", "Sandoval_FGK_base3Gyr", "Gaidos_rho050_vtan",
         "strong_2x", "Gaidos_raw_true_age", "upper"]
BETA = {"null": 0.0, "half_primary": -0.07, "Kamulali_primary": -0.1419,
        "Gaidos_rho025_vtan": -0.2206, "Berger_FGK": -0.3089,
        "Sandoval_FGK_base3Gyr": -0.1049, "Gaidos_rho050_vtan": -0.4412,
        "strong_2x": -0.60, "Gaidos_raw_true_age": -0.8824, "upper": -1.2}

def load(tag):
    p = CK / f"curve_{tag}.json"
    if not p.exists(): return None
    blob = json.loads(p.read_text())
    df = pd.DataFrame(blob["rows"])
    df = df[df.beta_pooled.notna()]
    return df, blob.get("done", len(df))

nul, n_done = load("null")
nmed, nsd = float(np.median(nul.beta_pooled)), float(np.std(nul.beta_pooled))
real = 0.04576344675022682   # frozen real cool-host pooled beta from archived gate summary
n025, n975 = np.percentile(nul.beta_pooled, [2.5, 97.5])

points, raw_parts = [], [nul.assign(point="null", beta_inj=0.0, arm="null")]
completed = 0
for i, name in enumerate(NAMES[1:]):
    got = load(f"sig_{i}")
    if got is None: continue
    sig, done = got
    if len(sig) < 150: continue
    completed += 1
    s = sig.beta_pooled.values
    smed = float(np.median(s))
    sep = (smed - nmed) / nsd
    pn = float(np.mean(nul.beta_pooled.values <= real))
    ps = float(np.mean(s <= real))
    verdict = "PASS" if (pn <= 0.05 and 0.025 <= ps <= 0.975) else "FAIL"
    det = float(np.mean(s < n025))
    points.append(dict(point=name, beta_inj=BETA[name], reps_null=len(nul),
                       reps_signal=len(s), null_median=nmed, null_sd=nsd,
                       signal_median=smed, arm_separation_nullSD=sep,
                       real_pct_null=pn, real_pct_signal=ps,
                       detect_rate_one_sided=det, verdict=verdict,
                       checkpoint_done=done))
    raw_parts.append(sig.assign(point=name, beta_inj=BETA[name], arm="signal"))

raw = pd.concat(raw_parts, ignore_index=True)
cols = ["point", "arm", "beta_inj", "rep", "beta_pooled", "se_pooled",
        "beta_strat", "se_strat", "beta_lo", "beta_hi", "n"]
raw[[c for c in cols if c in raw.columns]].to_csv(SCI / "power_vs_effect_raw.csv", index=False)

fig = pd.DataFrame([{k: p[k] for k in ("point","beta_inj","arm_separation_nullSD",
       "real_pct_null","real_pct_signal","detect_rate_one_sided","verdict","reps_signal")}
       for p in points])
fig.to_csv(SCI / "power_vs_effect_figure_data.csv", index=False)

missing = [n for n in NAMES[1:] if n not in [p["point"] for p in points]]
summary = dict(status=("PARTIAL - stopped by orchestrator; null arm + %d signal points "
                       "completed from checkpoints; %d requested points not run" % (completed, len(missing))),
               points_completed=completed + 1, signal_points_completed=completed,
               points_missing=missing,
               note=("Gate verdict is FAIL at every effect size because the real cool-host "
                     "coefficient (+0.0458) sits at the %.1f percentile of the null (needs <=0.05 "
                     "to pass distinguishability). This holds for Gaidos raw -0.8824 and attenuated "
                     "-0.2206/-0.4412 by construction; those arms have no replicates in this partial run."
                     % (100 * float(np.mean(nul.beta_pooled.values <= real)))),
               sample="M (cool / late-K early-M hosts)", estimator="pooled",
               strat_col="sy_vmag", reps_null=len(nul), seed=20260826,
               real_beta=real, null_median=nmed, null_sd=nsd, points=points)
(SCI / "power_vs_effect_summary.json").write_text(json.dumps(summary, indent=2))

try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fp, ax = plt.subplots(figsize=(7.2,4.8))
    ax.plot(fig.beta_inj, fig.arm_separation_nullSD, "o-", color="#1f77b4", label="Arm separation (null SD)")
    ax.axhline(0, color="0.6", lw=0.8); ax.grid(alpha=0.25)
    ax.axvline(-0.2206, color="#d62728", ls="--", lw=1, label="Gaidos rho=0.25")
    ax.set_xlabel("Injected effect beta per SD of v_tan"); ax.set_ylabel("Signal minus null (null SD)")
    ax.set_title("Cool-host gate power vs effect (partial)"); ax.legend(fontsize=8)
    fp.tight_layout(); fp.savefig(SCI / "power_vs_effect.png", dpi=200); plt.close(fp)
except Exception as e:
    print("plot failed:", e)

print(json.dumps(dict(points_completed=completed+1, missing=missing,
                      points=[{k:p[k] for k in ("point","beta_inj","reps_signal","arm_separation_nullSD","verdict")} for p in points]), indent=2))
