"""
13_calibration_fgk.py — TASK 2 GATE: calibrated FGK control re-validation.

Protocol (per directive):
  REAL   : re-run the headline age-evolution fit on the FGK control WITH
           completeness weighting applied -> beta_real.
  (a)    : inject a literature-magnitude synthetic effect into synthetic
           samples mirroring the weighted FGK noise structure; fit each ->
           signal-injected distribution of beta_hat.
  (b)    : inject zero under identical conditions -> null-injected
           distribution.
  Compare beta_real against BOTH distributions (empirical upper-tail p vs
  null; containment within signal distribution).
  PASS = consistent with (a) AND distinguishable from (b).
  FAIL = indistinguishable from (b), or from BOTH -> "inconclusive".

Injection magnitude translation (documented assumption): SWEET-Cat SE/SN
ratio 0.51 -> 0.64 between age halves == -0.224 logits spread over ~1.6 SD
of v_tan == **-0.14 per SD**, the same value used in the pre-registered
power analysis ("SWEET" scenario).

Also runs Task 2b part 1: weighted vs unweighted BINNED tests for
reconciliation with the hierarchical result.
"""
import importlib.util
import json
import pathlib
import numpy as np
import pandas as pd
from scipy.stats import norm

_spec = importlib.util.spec_from_file_location(
    "h9", pathlib.Path(__file__).parent / "09_hierarchical_model.py")
m9 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m9)

D = pathlib.Path(__file__).resolve().parents[1] / "data"
R = pathlib.Path(__file__).resolve().parents[1] / "results"
B_INJ = -0.14          # SWEET-Cat-equivalent, per SD of v_tan
REPS = 120


def load_weighted(tag="FGK"):
    pl = pd.read_csv(D / f"planet_sample_{tag}.csv")
    k = pd.read_csv(D / f"hosts_kinematics_{tag}.csv")
    w = pd.read_csv(D / f"completeness_weights_{tag}.csv")
    df = pl.merge(k[["hostname", "vtan", "W"]], on="hostname", how="inner")
    df = df.merge(w[["pl_name", "w", "p_det"]], on="pl_name", how="left")
    df["w"] = df.w.fillna(1.0)
    return df


def weighted_binned(df, label):
    """Horvitz-Thompson weighted two-proportion test across vtan median."""
    med = df.vtan.median()
    out = {}
    for nm, g in [("young", df[df.vtan <= med]), ("old", df[df.vtan > med])]:
        y = g.is_sn.values * g.wgt.values
        sw = g.wgt.sum()
        p = y.sum() / sw
        var = (g.wgt.values ** 2 * (g.is_sn.values - p) ** 2).sum() / sw ** 2
        out[nm] = (p, np.sqrt(var))
    z = (out["old"][0] - out["young"][0]) / np.sqrt(out["young"][1] ** 2
                                                    + out["old"][1] ** 2)
    pval = 2 * norm.sf(abs(z))
    print(f"[{label}] SN young={out['young'][0]:.3f}+-{out['young'][1]:.3f} "
          f"old={out['old'][0]:.3f}+-{out['old'][1]:.3f} z={z:+.2f} p={pval:.4f}")
    return z, pval


def synth(base, beta_inj, alpha, rng):
    picked = rng.choice(base.hostname.unique(),
                        size=base.hostname.nunique(), replace=True)
    pos = {h: i for i, h in enumerate(picked)}
    parts = [base[base.hostname == h] for h in picked]
    d = pd.concat(parts, ignore_index=True)
    zv = (d.vtan.values - base.vtan.mean()) / base.vtan.std()
    zp = (np.log10(d.pl_orbper.values) - np.log10(base.pl_orbper).mean()) \
        / np.log10(base.pl_orbper).std()
    # reuse measured gamma_logP from the unweighted FGK fit (~+2.3)
    eta = alpha + beta_inj * zv + 2.3 * zp
    uh = rng.normal(0, 1.5, len(picked))
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


def calibrate_alpha(df, target_frac, rng):
    lo, hi = -6.0, 6.0
    for _ in range(28):
        mid = 0.5 * (lo + hi)
        frac = np.mean(1 / (1 + np.exp(-(mid + 2.3 *
                     (np.log10(df.pl_orbper) - np.log10(df.pl_orbper).mean())
                     / np.log10(df.pl_orbper).std()))))
        if frac < target_frac:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def main():
    import sys
    subset = sys.argv[1] if len(sys.argv) > 1 else "all"       # all | kepler_only
    trim = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0     # e.g. 0.95
    df = load_weighted("FGK").rename(columns={"w": "wgt"})
    if subset == "kepler_only":
        df = df[df.hostname.str.startswith("Kepler")]
    if trim > 0:
        cap = np.quantile(df.wgt, trim)
        df["wgt"] = np.minimum(df.wgt, cap)
    ess = df.wgt.sum() ** 2 / (df.wgt ** 2).sum()
    print(f"weighted FGK [{subset}, trim={trim}]: N={len(df)}, ESS={ess:.0f}")
    globals()["TAG_LABEL"] = f"{subset}_trim{trim}"

    # ---- REAL weighted fit ----
    mdl = m9.build_design(df, "vtan", valley=1.88, fe=False,
                          weights=df.wgt.values)
    theta, se0, _, _ = m9.fit(mdl, fixed_su=True, fix_comp=True)
    se = m9.sandwich_se(mdl, theta, fixed_su=True)
    sd = df.vtan.std(); iqr = df.vtan.quantile(0.75) - df.vtan.quantile(0.25)
    b_real = theta[1] * (iqr / sd)
    se_real = se[1] * (iqr / sd)
    print(f"\nREAL weighted FGK: beta_age/IQR = {b_real:+.3f} +- {se_real:.3f}")

    # ---- calibration distributions ----
    rng = np.random.default_rng(20260826)
    res = {"beta_real": float(b_real), "se_real": float(se_real),
           "reps": REPS}
    for name, binj in [("null", 0.0), ("signal", B_INJ)]:
        a = calibrate_alpha(df, 0.408, rng)
        bets = []
        for rep in range(REPS):
            ds = synth(df, binj, a, rng)
            md = m9.build_design(ds, "vtan", valley=1.88, fe=False,
                                 weights=ds.wgt.values)
            th, s0, _, _ = m9.fit(md, fixed_su=True, fix_comp=True)
            bets.append(th[1] * (iqr / sd))
            if rep % 40 == 0:
                print(f"  [{name}] rep {rep}/{REPS}", flush=True)
        bets = np.array(bets)
        # empirical position of beta_real in this distribution
        pct_below = float(np.mean(bets <= b_real))
        res[name] = {"median": float(np.median(bets)),
                     "q025": float(np.quantile(bets, .025)),
                     "q975": float(np.quantile(bets, .975)),
                     "pct_rank_of_real": pct_below}
        print(f"[{name}-injected] median {np.median(bets):+.3f} "
              f"[{np.quantile(bets,.025):+.3f},{np.quantile(bets,.975):+.3f}]"
              f" | real percentile {pct_below:.3f}")

    # ---- PASS/FAIL per protocol ----
    p_null_upper = 1.0 - res["null"]["pct_rank_of_real"]      # P(null >= real)
    inside_sig = (res["signal"]["q025"] <= b_real <= res["signal"]["q975"])
    res["verdict"] = {
        "distinguishable_from_null": bool(p_null_upper <= 0.05),
        "consistent_with_signal": bool(inside_sig),
    }
    res["verdict"]["PASS"] = bool(res["verdict"]["distinguishable_from_null"]
                                  and inside_sig)
    print("\n=== TASK 2 VERDICT ===")
    print(f"P(null-injected >= real): {p_null_upper:.4f}")
    print(f"real inside 95% of signal-injected distribution: {inside_sig}")
    print("VERDICT:", "PASS" if res["verdict"]["PASS"] else
          ("FAIL(inconclusive)" if not
           res["verdict"]["distinguishable_from_null"] else
           "FAIL(null-like)"))

    # ---- Task 2b part 1: weighted vs unweighted binned tests ----
    print("\n=== TASK 2b: binned-test reconciliation ===")
    df["is_sn"] = df.pl_rade > 1.88
    weighted_binned(df, "weighted")
    dfu = df.copy(); dfu["wgt"] = 1.0
    weighted_binned(dfu, "unweighted")

    lab = globals().get("TAG_LABEL", "all_trim0")
    with open(R / f"task2_calibration_{lab}.json", "w") as fh:
        json.dump(res, fh, indent=1)
    print(f"saved {R/f'task2_calibration_{lab}.json'}")


if __name__ == "__main__":
    main()
