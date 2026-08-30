"""
B2_isochrone_gate.py — Phase B calibrated gate using published isochrone ages
(specifically controlling for the period-mixing confound).

Re-runs the exact mechanism from Task 2 (13_calibration_fgk.py) with TWO changes:
  1) The age sequence is now st_age (published isochrone/gyro ages) instead of v_tan.
  2) log_10(P) is added as a strict linear covariate in the hierarchical model,
     both for the real fits and synthetic calibration.
"""
import pathlib
import numpy as np
import pandas as pd
import json

import importlib.util
spec = importlib.util.spec_from_file_location("h9", pathlib.Path(__file__).parent / "09_hierarchical_model.py")
m9 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m9)

D = pathlib.Path(__file__).resolve().parents[1] / "data"
R = pathlib.Path(__file__).resolve().parents[1] / "results"
B_INJ = -0.14
REPS = 120

def load_data():
    pl = pd.read_csv(D / "planet_sample_FGK.csv")
    w = pd.read_csv(D / "completeness_weights_FGK.csv")
    ag = pd.read_csv(D / "published_ages_FGK.csv")

    # The archive-wide ages table is per-planet, not per-host -> collapse to
    # one row per hostname (taking the FIRST age and SMALLEST positive error
    # across duplicate host rows, which are identical for nominal FGK Kepler
    # hosts anyway). Without this dedup, the inner join multiplies planets.
    ag_k = (ag.sort_values(["hostname", "st_ageerr1"])
              .drop_duplicates(subset="hostname", keep="first")
              [["hostname", "st_age", "st_ageerr1", "st_ageerr2"]])

    df = pl.merge(w[["pl_name", "w", "p_det"]], on="pl_name", how="left")
    df = df.merge(ag_k[["hostname", "st_age", "st_ageerr1", "st_ageerr2"]],
                  on="hostname", how="inner")

    df["wgt"] = df.w.fillna(1.0)
    df = df[df.hostname.str.startswith("Kepler")]
    df = df.dropna(subset=["st_age"]).reset_index(drop=True)

    # Per-host age error (prefer upper-sigma column); fallback to ~25% of age
    df["st_age_err"] = df.st_ageerr1.fillna(df.st_ageerr2)
    df["st_age_err"] = df.st_age_err.fillna(0.25 * df.st_age).clip(lower=1e3)

    # Sample definition (final v3): FGK Kepler hosts with both non-null
    # st_age and st_age_err. Earlier B1 had 1065/1250 because the original
    # `ps` join returned fewer rows on st_age only; the present pull adds
    # ~600 more hosts. Documented deviation from the protocol's "preserve
    # sample definition" -- flagged for the report.
    df = df.dropna(subset=["st_age_err"]).reset_index(drop=True)

    # Trim weights (per protocol, untouched)
    cap = np.quantile(df.wgt, 0.95)
    df["wgt"] = np.minimum(df.wgt, cap)

    return df

def synth(base, beta_inj, alpha, rng):
    picked = rng.choice(base.hostname.unique(), size=base.hostname.nunique(), replace=True)
    pos = {h: i for i, h in enumerate(picked)}
    parts = []
    for h in picked:
        p_df = base[base.hostname == h]
        parts.append(p_df)
    d = pd.concat(parts, ignore_index=True)
    d = d.reset_index(drop=True)
    
    # Draw a host-level age shift consistent with each host's measured
    # Berger+20 age uncertainty (per-protocol age uncertainty model)
    age_shift = np.array([
        rng.normal(0.0, base.loc[base.hostname==h, "st_age_err"].values[0] / 1000.0)
        for h in d.hostname.values
    ])
    d["st_age"] = d.st_age.values + age_shift
    
    za = (d.st_age.values - d.st_age.mean()) / d.st_age.std()
    zp = (np.log10(d.pl_orbper.values) - np.log10(d.pl_orbper).mean()) / np.log10(d.pl_orbper).std()
    lum = np.clip(d.st_rad.fillna(0.6), 0.1, 4)**2 * (np.clip(d.st_teff.fillna(4000), 2500, 7500)/5772.0)**4
    a_au = np.cbrt(np.clip(d.st_mass.fillna(0.5), 0.08, 1.6)) * (d.pl_orbper.values / 365.25)**(2/3)
    zs = (np.log10(lum/a_au**2) - np.log10(lum/a_au**2).mean()) / np.log10(lum/a_au**2).std()
    eta = alpha + beta_inj * za + 2.30 * zp - 0.49 * zs - 4.08 * (d.hostname.str.startswith("Kepler").astype(float))
    
    uh = rng.normal(0, 1.5, len(picked))
    eta += uh[d.hostname.map(pos).values]
    
    pi = 1 / (1 + np.exp(-np.clip(eta, -30, 30)))
    sn = rng.uniform(size=len(d)) < pi
    r_true = np.where(sn, rng.normal(2.10, 0.35, len(d)), rng.normal(1.20, 0.18, len(d)))
    
    e = d.pl_radeerr1.abs()
    e = np.where(np.isfinite(e) & (e > 0.01), e, 0.05 * d.pl_rade)
    e = np.clip(e, 0.01 * d.pl_rade, 0.30 * d.pl_rade)  # outlier protection
    d["pl_rade"] = r_true + rng.normal(0, 1, len(d)) * e
    
    return d

def calibrate_alpha(df, target_frac, rng):
    lo, hi = -6.0, 6.0
    for _ in range(28):
        mid = 0.5 * (lo + hi)
        frac = np.mean(1 / (1 + np.exp(-(mid + 2.38 * 
                     (np.log10(df.pl_orbper) - np.log10(df.pl_orbper).mean()) / np.log10(df.pl_orbper).std()))))
        if frac < target_frac: lo = mid
        else: hi = mid
    return 0.5 * (lo + hi)

def main():
    df = load_data()
    print(f"Phase B Control Set (Isochrone Ages): N={len(df)} planets / {df.hostname.nunique()} hosts")
    
    # 1. REAL FIT
    mdl = m9.build_design(df, "st_age", valley=1.88, fe=False, weights=df.wgt.values)
    theta, se0, _, _ = m9.fit(mdl, fixed_su=True, fix_comp=True)
    se = m9.sandwich_se(mdl, theta, fixed_su=True)
    sd = df.st_age.std(); iqr = df.st_age.quantile(0.75) - df.st_age.quantile(0.25)
    
    b_real = theta[1] * (iqr / sd)
    se_real = se[1] * (iqr / sd)
    b_per = theta[2]  # Gamma for logP
    print(f"\nREAL FIT: beta_age/IQR = {b_real:+.3f} +- {se_real:.3f} | beta_logP = {b_per:+.3f}")

    # 2. CALIBRATION
    rng = np.random.default_rng(20261109)
    res = {"beta_real": float(b_real), "se_real": float(se_real), "reps": REPS,
           "per_rep_betas": {"null": [], "signal": []},
           "per_rep_ses":   {"null": [], "signal": []}}

    for name, binj in [("null", 0.0), ("signal", B_INJ)]:
        a = calibrate_alpha(df, 0.408, rng)
        bets = []
        ses  = []
        for rep in range(REPS):
            ds = synth(df, binj, a, rng)
            md = m9.build_design(ds, "st_age", valley=1.88, fe=False, weights=ds.wgt.values)
            th, _, _, _ = m9.fit(md, fixed_su=True, fix_comp=True)
            se_rep = m9.sandwich_se(md, th, fixed_su=True)
            bets.append(th[1] * (iqr / sd))
            ses.append(se_rep[1] * (iqr / sd))
            if rep % 40 == 0: print(f"  [{name}] rep {rep}/{REPS}", flush=True)

        bets = np.array(bets)
        ses  = np.array(ses)
        res["per_rep_betas"][name] = bets.tolist()
        res["per_rep_ses"][name]   = ses.tolist()
        pct_below = float(np.mean(bets <= b_real))
        res[name] = {"median": float(np.median(bets)),
                     "q025": float(np.quantile(bets, .025)),
                     "q975": float(np.quantile(bets, .975)),
                     "pct_rank_of_real": pct_below,
                     "per_rep_se_mean": float(np.mean(ses)),
                     "per_rep_se_median": float(np.median(ses)),
                     "empirical_sd": float(np.std(bets, ddof=1))}
        print(f"[{name}-injected] median {np.median(bets):+.3f} "
              f"[{np.quantile(bets,.025):+.3f},{np.quantile(bets,.975):+.3f}]"
              f" | real percentile {pct_below:.3f}")

    # 3. VERDICT
    p_null_upper = 1.0 - res["null"]["pct_rank_of_real"]
    inside_sig = (res["signal"]["q025"] <= b_real <= res["signal"]["q975"])
    res["verdict"] = {"distinguishable_from_null": bool(p_null_upper <= 0.05),
                      "consistent_with_signal": bool(inside_sig)}
    res["verdict"]["PASS"] = bool(res["verdict"]["distinguishable_from_null"] and inside_sig)
    
    print("\n=== PHASE B VERDICT ===")
    print(f"P(null-injected >= real): {p_null_upper:.4f}")
    print(f"real inside 95% of signal-injected distribution: {inside_sig}")
    print("VERDICT:", "PASS" if res["verdict"]["PASS"] else
          ("FAIL(inconclusive)" if not res["verdict"]["distinguishable_from_null"] else "FAIL(null-like)"))

if __name__ == "__main__":
    main()
