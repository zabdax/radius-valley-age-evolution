"""
11_power_analysis.py — Could ANY feasible M-dwarf sample detect the FGK-size
radius-valley age effect?

Method
------
Generative engine calibrated on the REAL M sample:
  - bootstrap hosts (with their planets, covariates, v_tan, radius errors)
    at scale factors k in {0.5, 1, 2, 4, 8} x the observed 417 planets,
    preserving cluster structure;
  - inject an age effect beta_inj (per SD of v_tan) into the linear
    predictor alongside the covariate coefficients measured on real data;
  - draw class labels, then observed radii from overlapping class-
    conditionals (N(1.30,0.22) / N(2.35,0.45)) plus per-planet errors ->
    soft-classification weights behave exactly as in the real pipeline;
  - estimate with the full hierarchical pipeline (09) and record whether
    the 95% CI excludes zero with the injected sign ("detection").

Effect sizes (per SD of v_tan):
  null        :  0.00
  SWEET-Cat   : -0.14   (= SN fraction falling ~5 pp across age halves,
                          logit shift 0.22 spread over 1.6 sigma)
  strong      : -0.28   (2x SWEET-Cat, Berger/David-flavoured)

Outputs: results/power_analysis.json + console table.
Caveat (documented): upsampling resamples the empirical covariate/host
distribution rather than creating independent new hosts -- optimistic by
some unknown amount; treat required-N as lower bounds.
"""
import importlib.util
import json
import pathlib
import numpy as np
import pandas as pd

_spec = importlib.util.spec_from_file_location(
    "h9", pathlib.Path(__file__).parent / "09_hierarchical_model.py")
m9 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m9)

D = pathlib.Path(__file__).resolve().parents[1] / "data"
R = pathlib.Path(__file__).resolve().parents[1] / "results"

# generative coefficients (from the real M x vtan fit, 09)
GAMMA_LOGP, GAMMA_LOGS = 1.613, -0.032
GAMMA_KEP, GAMMA_K2 = -1.594, 0.893
SIGMA_U = np.exp(-1.5)          # placeholder; overwritten below if available
# reality-matched class-conditionals (audit finding: previous N(1.30,0.22)/
# N(2.35,0.45) gave ~2x too little cross-valley ambiguity vs the real
# M-sample modes 1.23/2.08 and boundary-window mass ~28%)
MU_SE, S_SE = 1.20, 0.18
MU_SN, S_SN = 2.10, 0.35
VALLEY = 1.85


def load_base():
    pl = pd.read_csv(D / "planet_sample_M.csv")
    k = pd.read_csv(D / "hosts_kinematics_M.csv")
    df = pl.merge(k[["hostname", "vtan", "W"]], on="hostname", how="inner").reset_index(drop=True)
    err = df.pl_radeerr1.abs()
    e = np.where(np.isfinite(err) & (err > 0.01), err, 0.05 * df.pl_rade)
    df["e_r"] = e
    df["mission"] = df.hostname.str.extract(
        r"^(TOI|Kepler|KIC|K2|HIP|HD|GJ|Wolf|LHS|TRAPPIST|EPIC)",
        expand=False).fillna("other")
    df["fam"] = df.mission.map(lambda s: s if s in ("Kepler", "K2", "TOI") else "other")
    return df


def zscore_within(df, col):
    v = df[col].values.astype(float)
    return (v - v.mean()) / (v.std() + 1e-12)


def synth(base, scale, beta_inj, alpha, rng):
    """Draw one synthetic dataset."""
    hosts = base.hostname.unique()
    n_hosts = max(1, int(round(scale * len(hosts))))
    picked = rng.choice(hosts, size=n_hosts, replace=True)
    idx = np.concatenate([np.flatnonzero(base.hostname.values == h) for h in picked])
    d = base.iloc[idx].reset_index(drop=True)
    zv = zscore_within(d, "vtan")
    zp = zscore_within(d.assign(lp=np.log10(d.pl_orbper)), "lp")
    lum = np.clip(d.st_rad.fillna(0.6), 0.1, 4) ** 2 * \
        (np.clip(d.st_teff.fillna(4000), 2500, 7500) / 5772.0) ** 4
    a_au = np.cbrt(np.clip(d.st_mass.fillna(0.5), 0.08, 1.6)) * \
        (d.pl_orbper.values / 365.25) ** (2 / 3)
    d["ls"] = np.log10(lum / a_au ** 2)
    zs = zscore_within(d, "ls")
    eta = (alpha + beta_inj * zv + GAMMA_LOGP * zp + GAMMA_LOGS * zs
           + GAMMA_KEP * (d.fam == "Kepler").values.astype(float)
           + GAMMA_K2 * (d.fam == "K2").values.astype(float))
    # map rows to their DRAWN host index (not categorical-sorted codes!)
    pos = {h: i for i, h in enumerate(picked)}
    codes = d.hostname.map(pos).values
    uh = rng.normal(0, SIGMA_U, len(picked))
    eta += uh[codes]
    pi = 1 / (1 + np.exp(-np.clip(eta, -30, 30)))
    sn = rng.uniform(size=len(d)) < pi
    r_true = np.where(sn, rng.normal(MU_SN, S_SN, len(d)),
                      rng.normal(MU_SE, S_SE, len(d)))
    d["pl_rade"] = r_true + rng.normal(0, 1, len(d)) * d.e_r.values
    return d


def calibrate_alpha(base, beta_inj, target_frac, rng, mission_offset=0.0):
    """Bisect alpha so the pre-radius SN probability averages target_frac.
    Uses one planet per host (ignores cluster size) -- adequate for a
    calibration constant. mission_offset adds the mean contribution of the
    mission dummies (audit finding: omitting it biased the realized rate)."""
    lo, hi = -6.0, 6.0
    for _ in range(30):
        mid = 0.5 * (lo + hi)
        hosts = base.hostname.unique()
        picked = rng.choice(hosts, size=len(hosts), replace=True)
        first_idx = [np.flatnonzero(base.hostname.values == h)[0] for h in picked]
        d1 = base.iloc[first_idx]
        zv = (d1.vtan.values - base.vtan.mean()) / base.vtan.std()
        zp = (np.log10(d1.pl_orbper.values) - np.log10(base.pl_orbper).mean()) \
            / np.log10(base.pl_orbper).std()
        eta = (mid + mission_offset + beta_inj * zv + GAMMA_LOGP * zp)
        frac = np.mean(1 / (1 + np.exp(-eta)))
        if frac < target_frac:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def main(reps_null=120, reps_eff=150):
    base = load_base()
    real_frac = float((base.pl_rade > VALLEY).mean())
    # calibrate pre-radius prob slightly high: class-overlap noise in the
    # radius draw pulls the OBSERVED valley-crossing fraction below it
    calib_frac = min(0.495, real_frac / 0.90)
    print(f"real M SN fraction: {real_frac:.3f}; "
          f"{base.hostname.nunique()} hosts, {len(base)} planets")

    # Generative host heterogeneity scenario. The v2 model fixes sigma_u=0
    # (sandwich SEs), so there is no fitted value to read; we adopt the
    # host-overdispersion scale measured by the v1 hierarchical fit
    # (~1.5 logits, itself a sqrt2-deflated estimate of ~2.2) as the base
    # scenario. Larger true heterogeneity would LOWER power further.
    global SIGMA_U
    SIGMA_U = 1.5
    print(f"generative sigma_u scenario: {SIGMA_U:.2f} logits")

    cells = []
    # stable seeds (audit finding: salted hash() was non-reproducible)
    eff_seed = {"null": 0, "SWEET": 1, "strong": 2}
    for scale, label in [(0.5, "0.5x"), (1.0, "1x"), (2.0, "2x"),
                         (4.0, "4x"), (8.0, "8x")]:
        for eff_name, beta_inj in [("null", 0.0),
                                   ("SWEET", -0.14),
                                   ("strong", -0.28)]:
            reps = reps_null if eff_name == "null" else reps_eff
            seed = eff_seed[eff_name] * 1000 + int(round(scale * 100))
            cells.append((scale, label, eff_name, beta_inj, reps, seed))

    results = []
    f_kep = float((base.fam == "Kepler").mean())
    f_k2 = float((base.fam == "K2").mean())
    mission_offset = GAMMA_KEP * f_kep + GAMMA_K2 * f_k2
    print(f"mission offset added to calibration: {mission_offset:+.3f} logits")
    for scale, label, eff_name, beta_inj, reps, seed in cells:
        rng = np.random.default_rng(seed)
        alpha = calibrate_alpha(base, beta_inj, calib_frac, rng,
                                mission_offset=mission_offset)
        det = cov_cnt = 0
        bets = []
        for r in range(reps):
            d = synth(base, scale, beta_inj, alpha, rng)
            mdl = m9.build_design(d, "vtan", VALLEY, fe=False)
            theta, _, _, _ = m9.fit(mdl, fixed_su=True)
            se = m9.sandwich_se(mdl, theta, fixed_su=True)
            b, s = theta[1], se[1]
            bets.append(b)
            # CI excludes zero iff both bounds share a sign
            excl0 = (b - 1.96 * s) * (b + 1.96 * s) > 0
            covers = (b - 1.96 * s <= beta_inj) and (beta_inj <= b + 1.96 * s)
            if beta_inj == 0:
                det += excl0
            else:
                det += excl0 and (np.sign(b) == np.sign(beta_inj))
            cov_cnt += covers
        results.append({
            "scale": label, "effect": eff_name, "beta_inj": beta_inj,
            "reps": reps, "detect_rate": det / reps,
            "coverage": cov_cnt / reps,
            "median_beta_hat": float(np.median(bets))})
        print(f"[{label:>4} | {eff_name:>6}] detect={det/reps:5.1%} "
              f"coverage={cov_cnt/reps:5.1%} med_bhat={np.median(bets):+.3f}",
              flush=True)
        with open(R / "power_analysis.json", "w") as fh:
            json.dump(results, fh, indent=1)


if __name__ == "__main__":
    main()
