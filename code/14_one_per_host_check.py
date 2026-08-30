"""
14_one_per_host_check.py — Bounded robustness check: within-host clustering.

Review Issue C: planets sharing a host star are not independent draws (shared
host age, formation history, disk composition). Does the FGK gate's
inconclusive verdict depend on how within-host structure is handled?

Code-inspection context (documented, verified in 09/11/13):
  - Production fits use fit(fixed_su=True): the mixture likelihood's host
    random intercept is fixed at sigma_u = 0 (working independence), chosen
    deliberately to avoid the slope-vs-random-effect identifiability
    competition documented in 09's docstrings.
  - All reported SEs are HOST-LEVEL CLUSTER-ROBUST: 09.sandwich_se()
    aggregates per-planet scores to per-host scores (bincount on hidx)
    before forming the sandwich meat, V = H^-1 (sum_g s_g s_g^T) H^-1.
  - The injection machinery (11_power_analysis.py, 13_calibration_fgk.py)
    bootstraps WHOLE HOSTS and draws ONE shared host effect per host
    (sigma_u = 1.5 logits), so null/signal distributions inherit
    within-host correlation.

This check empirically tests materiality: re-run the primary FGK gate using
ONLY ONE PLANET PER HOST and compare beta/SE/percentile placement against
the full-sample result (beta = -0.030 +- 0.106 per IQR; 43.3rd percentile of
the null-injected distribution; verdict FAIL(inconclusive)).

Selection rule (documented): PRIMARY = LARGEST planet (max pl_rade) per host,
maximizing per-planet S/N; SENSITIVITY = uniformly random planet per host
(fixed seed). Bounded scope: one real fit per sample + REPS null-injected
replicates for the primary sample only. Verdict criterion identical to 13:
distinguishable from null iff P(null >= beta_real) <= 0.05.

Run:  py -3 14_one_per_host_check.py [reps] [smoke]
Outputs: results/one_per_host_check.json + console summary.
"""
import importlib.util
import json
import pathlib
import sys
import numpy as np
import pandas as pd

_spec = importlib.util.spec_from_file_location(
    "h9", pathlib.Path(__file__).parent / "09_hierarchical_model.py")
m9 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m9)

D = pathlib.Path(__file__).resolve().parents[1] / "data"
R = pathlib.Path(__file__).resolve().parents[1] / "results"
SMOKE = "smoke" in [a.lower() for a in sys.argv[1:]]
REPS_NULL = 3 if SMOKE else int(next((a for a in sys.argv[1:]
                                      if a.isdigit()), 120))
GAMMA_LOGP = 2.3        # measured gamma_logP, unweighted FGK fit (as in 13)
SIGMA_U = 1.5           # shared host-effect SD, as in 11/13
VALLEY = 1.88
RNG_SEED = 20260829


def load_weighted_fgk_kepler():
    """Identical to 13_calibration_fgk.py: FGK, Kepler-only, trim 0.95."""
    pl = pd.read_csv(D / "planet_sample_FGK.csv")
    k = pd.read_csv(D / "hosts_kinematics_FGK.csv")
    w = pd.read_csv(D / "completeness_weights_FGK.csv")
    df = pl.merge(k[["hostname", "vtan", "W"]], on="hostname", how="inner")
    df = df.merge(w[["pl_name", "w", "p_det"]], on="pl_name", how="left")
    df["w"] = df.w.fillna(1.0)
    df = df[df.hostname.str.startswith("Kepler")].copy()
    cap = np.quantile(df.w, 0.95)
    df["w"] = np.minimum(df.w, cap)
    return df.rename(columns={"w": "wgt"}).reset_index(drop=True)


def one_per_host(df, mode="largest", seed=RNG_SEED):
    """Reduce to one planet per host. mode='largest' (primary) keeps the
    largest-radius planet per host; mode='random' keeps a uniformly random
    one (documented sensitivity, fixed seed)."""
    if mode == "largest":
        sub = (df.sort_values("pl_rade", ascending=False)
                 .drop_duplicates("hostname")
                 .sort_values("hostname"))
    else:
        rng = np.random.default_rng(seed)
        rows = [g.iloc[rng.integers(0, len(g))]
                for _, g in df.groupby("hostname", sort=True)]
        sub = pd.DataFrame(rows)
    return sub.reset_index(drop=True)


def fit_real(df):
    """Production-identical real fit; returns per-IQR beta and sandwich SE."""
    mdl = m9.build_design(df, "vtan", valley=VALLEY, fe=False,
                          weights=df.wgt.values)
    theta, _, _, _ = m9.fit(mdl, fixed_su=True, fix_comp=True)
    se = m9.sandwich_se(mdl, theta, fixed_su=True)
    sd = df.vtan.std()
    iqr = df.vtan.quantile(0.75) - df.vtan.quantile(0.25)
    ess = float(df.wgt.sum() ** 2 / (df.wgt ** 2).sum())
    return (theta[1] * iqr / sd, se[1] * iqr / sd,
            len(df), int(mdl["n_hosts"]), ess)


def calibrate_alpha(df, target_frac):
    """Deterministic bisection on alpha (identical to 13)."""
    lo, hi = -6.0, 6.0
    for _ in range(28):
        mid = 0.5 * (lo + hi)
        frac = np.mean(1 / (1 + np.exp(-(mid + GAMMA_LOGP *
                    (np.log10(df.pl_orbper) - np.log10(df.pl_orbper).mean())
                    / np.log10(df.pl_orbper).std()))))
        if frac < target_frac:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def synth_null(base, alpha, rng):
    """One null synthetic realization; host-bootstrap + one shared host
    effect per host (identical structure to 13_calibration_fgk.synth with
    beta_inj = 0)."""
    picked = rng.choice(base.hostname.unique(),
                        size=base.hostname.nunique(), replace=True)
    pos = {h: i for i, h in enumerate(picked)}
    parts = [base[base.hostname == h] for h in picked]
    d = pd.concat(parts, ignore_index=True)
    zp = ((np.log10(d.pl_orbper.values)
           - np.log10(base.pl_orbper).mean())
          / np.log10(base.pl_orbper).std())
    eta = alpha + GAMMA_LOGP * zp
    uh = rng.normal(0, SIGMA_U, len(picked))
    eta += uh[d.hostname.map(pos).values]
    pi = 1 / (1 + np.exp(-np.clip(eta, -30, 30)))
    sn = rng.uniform(size=len(d)) < pi
    r_true = np.where(sn, rng.normal(2.10, 0.35, len(d)),
                      rng.normal(1.20, 0.18, len(d)))
    e = d.pl_radeerr1.abs()
    e = np.where(np.isfinite(e) & (e > 0.01), e, 0.05 * d.pl_rade)
    e = np.clip(e, 0.01 * d.pl_rade, 0.30 * d.pl_rade)
    d["pl_rade"] = r_true + rng.normal(0, 1, len(d)) * e
    return d


def null_gate(df, b_real, iqr, sd, reps):
    """Null-injected distribution; percentile placement of beta_real,
    exactly the 13 protocol's null arm."""
    frac_obs = float((df.pl_rade > VALLEY).mean())
    alpha = calibrate_alpha(df, min(0.495, frac_obs / 0.90))
    rng = np.random.default_rng(RNG_SEED)
    bets = []
    for rep in range(reps):
        ds = synth_null(df, alpha, rng)
        md = m9.build_design(ds, "vtan", valley=VALLEY, fe=False,
                             weights=ds.wgt.values)
        th, _, _, _ = m9.fit(md, fixed_su=True, fix_comp=True)
        bets.append(th[1] * iqr / sd)
        if rep % 10 == 0:
            print(f"    rep {rep}/{reps}", flush=True)
    bets = np.array(bets)
    pct_below = float(np.mean(bets <= b_real))
    return {"reps": reps, "median": float(np.median(bets)),
            "q025": float(np.quantile(bets, .025)),
            "q975": float(np.quantile(bets, .975)),
            "pct_rank_of_real": pct_below,
            "p_upper": float(1.0 - pct_below),
            "distinguishable_from_null": bool(1.0 - pct_below <= 0.05),
            "null_betas": [float(b) for b in bets]}


def main():
    out = {"protocol": "13_calibration_fgk: FGK, Kepler-only, trim 0.95, "
                       "fixed_su sandwich SE, per-IQR beta",
           "selection_primary": "largest planet per host (max pl_rade)",
           "selection_sensitivity": "random planet per host, seed 20260829",
           "reps_null": REPS_NULL, "smoke": SMOKE}
    full = load_weighted_fgk_kepler()
    b0, se0, n0, h0, ess0 = fit_real(full)
    print(f"FULL sample: N={n0}, ESS={ess0:.0f} (expect 2135 / 1645)  "
          f"beta={b0:+.3f}+-{se0:.3f}", flush=True)
    out["full"] = {"n_planets": n0, "n_hosts": h0, "ess": ess0,
                   "beta_per_iqr": b0, "se_per_iqr": se0, "z": b0 / se0}

    for mode in (["largest"] if SMOKE else ["largest", "random"]):
        sub = one_per_host(full, mode)
        b, se, n, h, ess = fit_real(sub)
        print(f"\nONE-PER-HOST [{mode}]: N={n} ({h} hosts), ESS={ess:.0f}  "
              f"beta={b:+.3f}+-{se:.3f}  z={b / se:+.2f}  "
              f"SE ratio vs full = {se / se0:.2f}", flush=True)
        entry = {"n_planets": n, "n_hosts": h, "ess": ess,
                 "beta_per_iqr": b, "se_per_iqr": se, "z": b / se,
                 "se_ratio_vs_full": se / se0}
        if mode == "largest":
            print(f"  null gate ({REPS_NULL} reps) starting...", flush=True)
            iqr = sub.vtan.quantile(0.75) - sub.vtan.quantile(0.25)
            sd = sub.vtan.std()
            g = null_gate(sub, b, iqr, sd, REPS_NULL)
            entry["null_gate"] = {k: v for k, v in g.items()
                                  if k != "null_betas"}
            entry["verdict"] = ("PASS" if g["distinguishable_from_null"]
                                else "FAIL(inconclusive)")
            print(f"  null median {g['median']:+.3f} "
                  f"[{g['q025']:+.3f},{g['q975']:+.3f}] | "
                  f"real percentile {g['pct_rank_of_real']:.3f} | "
                  f"p_upper {g['p_upper']:.3f} -> verdict "
                  f"{entry['verdict']}", flush=True)
        out[f"one_per_host_{mode}"] = entry
        R.mkdir(exist_ok=True)
        with open(R / "one_per_host_check.json", "w") as fh:
            json.dump(out, fh, indent=1)

    print("\nsaved results/one_per_host_check.json")


if __name__ == "__main__":
    main()