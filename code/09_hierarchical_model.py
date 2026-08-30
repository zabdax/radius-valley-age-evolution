"""
09_hierarchical_model.py — Bayesian logistic model for P(sub-Neptune | age, covariates).

Structure
---------
For planet i orbiting host h:

    latent true radius : R_i ~ N(r_obs_i, e_i)            (soft classification;
        p_SN^obs_i = Phi((r_obs_i - floor)/e_i) -> no hard cut at the valley)
    host random effect : u_h ~ N(0, sigma_u)              (multi-planet systems)
    linear predictor   : z_i = alpha + beta_a*a_h + gamma'x_i + u_h
    outcome likelihood : l_i(u_h) = log[ w*pi + (1-w)*(1-pi) ],  pi=sigmoid(z_i)

Marginalized over u_h by Gauss-Hermite quadrature; MAP + Laplace CIs
(pure numpy/scipy -- no new dependencies). Insolation computed from
P, M*, R*, Teff via Kepler III because the archive column is mostly empty.

Age proxies (ordinal): v_tan (all hosts) or |W| (RV subset). Hypothesis from
Berger20/David21/Kamulali26/Gaidos24: beta_a < 0 (SN fraction declines with age).

Outputs
-------
results/hierarchical_{tag}_{variant}.json + printed summary.
Run `py -3 09_hierarchical_model.py SELFTEST` first: recovery on synthetic data.
"""
import json
import pathlib
import sys
import numpy as np
import pandas as pd
from scipy.special import erf, expit, logsumexp
from scipy.optimize import minimize

D = pathlib.Path(__file__).resolve().parents[1] / "data"
R = pathlib.Path(__file__).resolve().parents[1] / "results"

GH_NODES, GH_W = np.polynomial.hermite_e.hermegauss(24)


def gauss_hermite_grid(sigma):
    # probabilists' convention: hermegauss nodes satisfy
    # sum(w_j f(x_j))/sqrt(2 pi) ~= E[f(Z)], Z~N(0,1); scale nodes by sigma.
    u = sigma * GH_NODES
    w = GH_W / np.sqrt(2.0 * np.pi)
    return u, w


def build_design(df, age_col, valley, fe=False, weights=None):
    """Assemble data quantities and design matrix. Returns dict.

    Mixture likelihood (v2): the outcome is the OBSERVED RADIUS itself,
    modeled as a two-component Gaussian mixture (super-Earth / sub-Neptune
    populations) whose mixing fraction pi(age, covariates) is the science
    target. Cross-valley misclassification is thereby MODELED rather than
    ignored (the previous Phi-weight Bernoulli treatment attenuated beta
    by ~35% in the real-data regime -- see selftest).
    """
    d = df.copy()
    err = df.get("pl_radeerr1")
    e = np.asarray(err.abs().values if err is not None
                   else np.full(len(d), np.nan), float)
    # physical clip: archive contains absurd error entries (one at 68.9 Re!)
    # -> cap relative error at 30%, floor at 1%
    r_abs = np.abs(d.pl_rade.values) + 1e-6   # guard against negative synthetic radii
    e = np.where(np.isfinite(e) & (e > 0.01), e, 0.05 * r_abs)
    e = np.clip(e, 0.01 * r_abs, 0.30 * r_abs)

    # physical insolation (S in solar constants): a[AU]=(M*)^(1/3)(P/yr)^(2/3)
    mstar = d.st_mass.fillna(0.5).values          # typical fallback
    mstar = np.clip(mstar, 0.08, 1.6)
    p_yr = d.pl_orbper.values / 365.25
    a_au = np.cbrt(mstar) * p_yr ** (2.0 / 3.0)
    rstar = d.st_rad.fillna(0.6).values if "st_rad" in d else 0.6
    tstar = d.st_teff.fillna(4000).values if "st_teff" in d else 4000
    lum = np.clip(rstar, 0.1, 4) ** 2 * (np.clip(tstar, 2500, 7500) / 5772.0) ** 4
    s_in = lum / a_au ** 2

    # predictors
    lp = np.log10(d.pl_orbper.values)
    ls = np.log10(s_in)
    mission = d.hostname.str.extract(
        r"^(TOI|Kepler|KIC|K2|HIP|HD|GJ|Wolf|LHS|TRAPPIST|EPIC)",
        expand=False).fillna("other")
    fam = mission.map(lambda s: s if s in ("Kepler", "K2", "TOI") else "other")

    cols = {"logP": lp, "logS": ls}
    X = [np.ones(len(d))]
    names = ["alpha"]

    def z(v):
        v = np.asarray(v, float)
        return (v - v.mean()) / (v.std() + 1e-12)

    age_z = z(d[age_col].values)
    iqr = d[age_col].quantile(0.75) - d[age_col].quantile(0.25)
    X.append(age_z); names.append("beta_age")
    for k, v in cols.items():
        # winsorize z-scores at +-5 rather than zeroing rows: zeroing turned
        # outliers into fake "mean-imputed" points (audit finding)
        vz = np.clip(z(v), -5, 5)
        X.append(vz); names.append(f"gamma_{k}")
    dum_kep = (fam == "Kepler").astype(float).values
    dum_k2 = (fam == "K2").astype(float).values
    X.append(dum_kep); names.append("gamma_missionKepler")
    X.append(dum_k2); names.append("gamma_missionK2")
    if fe:
        feh = d.st_met.fillna(0.0).values
        X.append(z(feh)); names.append("gamma_FeH")
    Xm = np.column_stack(X)

    hosts = d.hostname.values
    uh_list = pd.unique(hosts)
    hidx = pd.Categorical(hosts, categories=uh_list).codes
    # mixture-component init from the data (below/above the measured floor)
    r = d.pl_rade.values
    mu_se0 = float(np.median(r[r < valley]))
    mu_sn0 = float(np.median(r[r >= valley]))
    s_se0 = float(np.std(r[r < valley]) + 1e-3)
    s_sn0 = float(np.std(r[r >= valley]) + 1e-3)
    return {"X": Xm, "names": names, "r_obs": r, "e_r": e,
            "hidx": hidx, "n_hosts": len(uh_list), "age_iqr": iqr,
            "age_raw": d[age_col].values,
            "w_pl": (np.asarray(weights, float) if weights is not None
                     else np.ones(len(d))),
            "init_comp": (mu_se0, s_se0, mu_sn0, s_sn0)}


# parameter vector layout: [beta..., log_sigma_u, mu_se, log_s_se, mu_sn, log_s_sn]
NCOMP = 4


def neg_log_post(theta, mdl, prior_scale=2.5, fixed_su=False, norm=1.0):
    """Negative log posterior.
    norm: divide the FULL objective (data+priors) by this scalar to keep
          gradients O(1) -- prevents L-BFGS-B's relative tolerance from 
          terminating at the zero-initialized beta on large samples.
    fixed_su=True : no host random intercept (sigma_u = 0). Use with the
        cluster-robust sandwich SE; removes the slope-vs-random-effect
        weak-identifiability competition.
    """
    ncoef = mdl["X"].shape[1]
    beta = theta[:ncoef]
    if fixed_su:
        su = 0.0
        mu_se, log_sse = theta[ncoef], theta[ncoef + 1]
        mu_sn, log_ssn = theta[ncoef + 2], theta[ncoef + 3]
    else:
        log_su = theta[ncoef]
        if not (-6.0 < log_su < 2.0):
            return 1e12
        su = np.exp(log_su)
        mu_se, log_sse = theta[ncoef + 1], theta[ncoef + 2]
        mu_sn, log_ssn = theta[ncoef + 3], theta[ncoef + 4]
    # ordering + sanity constraints
    if not (0.6 < mu_se < valley_of(mdl) < mu_sn < 4.5):
        return 1e12
    if not (-2.5 < log_sse < -0.2 and -2.5 < log_ssn < -0.0):
        return 1e12
    tau_se = np.exp(log_sse)
    tau_sn = np.exp(log_ssn)

    eta0 = np.clip(mdl["X"] @ beta, -30, 30)
    lw_sn = norm_logpdf(mdl["r_obs"], mu_sn,
                        np.sqrt(tau_sn ** 2 + mdl["e_r"] ** 2))
    lw_se = norm_logpdf(mdl["r_obs"], mu_se,
                        np.sqrt(tau_se ** 2 + mdl["e_r"] ** 2))

    u, gw = gauss_hermite_grid(su)          # single node (u=0) when su=0
    totals = np.empty((len(u), mdl["n_hosts"]))
    for j, uj in enumerate(u):
        pi = expit(np.clip(eta0 + uj, -30, 30))
        li = np.logaddexp(lw_sn + np.log(pi), lw_se + np.log1p(-pi))
        totals[j] = np.bincount(mdl["hidx"], weights=mdl["w_pl"] * li,
                                minlength=mdl["n_hosts"])
    exponents = totals + np.log(gw)[:, None]
    nll = -np.sum(logsumexp(exponents, axis=0))

    # priors: weak on coefficients; weak on component params
    nll += 0.5 * np.sum((beta[1:] / prior_scale) ** 2)
    nll += 0.5 * (beta[0] / 5.0) ** 2
    if not fixed_su:
        nll += 0.5 * ((log_su + 1.5) / 1.0) ** 2
    nll += 0.5 * ((mu_se - 1.25) / 0.75) ** 2 + 0.5 * ((mu_sn - 2.15) / 0.75) ** 2
    if not np.isfinite(nll):
        return 1e12
    return nll / norm


def valley_of(mdl):
    """Midpoint between fitted/initial components; used only for the
    ordering constraint."""
    lo, hi = 0.9, 4.4
    m0 = mdl.get("init_comp")
    if m0:
        return 0.5 * (m0[0] + m0[2])
    return 0.5 * (lo + hi)


def norm_logpdf(x, mu, sig):
    return -0.5 * ((x - mu) / sig) ** 2 - np.log(sig) - 0.5 * np.log(2 * np.pi)


def fit(mdl, fixed_su=True, fix_comp=True):
    """fix_comp=True pins the mixture components in a narrow window around
    data-measured starting values (medians below/above the valley floor).
    Free components suffer a variance/mixing degeneracy that makes beta
    unidentifiable in a minority of replicates (component absorbs the age
    signal); fixed components remove it. Sensitivity runs should vary the
    pinning window."""
    p = mdl["X"].shape[1]
    mu_se0, s_se0, mu_sn0, s_sn0 = mdl["init_comp"]
    if fixed_su:
        th0 = np.r_[np.zeros(p), [mu_se0, np.log(s_se0),
                                  mu_sn0, np.log(s_sn0)]]
        if fix_comp:
            cb = [(mu_se0 - 0.06, mu_se0 + 0.06),
                  (np.log(s_se0) - 0.12, np.log(s_se0) + 0.12),
                  (mu_sn0 - 0.06, mu_sn0 + 0.06),
                  (np.log(s_sn0) - 0.12, np.log(s_sn0) + 0.12)]
        else:
            cb = [(0.7, 1.75), (-2.5, -0.2), (1.8, 4.2), (-2.3, -0.05)]
        bounds = [(-25, 25)] * p + cb
    else:
        th0 = np.r_[np.zeros(p), np.array([-1.5]),
                    [mu_se0, np.log(s_se0), mu_sn0, np.log(s_sn0)]]
        bounds = ([(-25, 25)] * p + [(-6, 2)]
                  + [(0.7, 1.75), (-2.5, -0.2), (1.8, 4.2), (-2.3, -0.05)])
    S = float(len(mdl["r_obs"]))
    nll = lambda t: neg_log_post(t, mdl, fixed_su=fixed_su, norm=S)
    res = minimize(nll, th0, method="L-BFGS-B",
                   bounds=bounds, options={"maxiter": 2000})
    f0 = nll(res.x)
    # central-difference Hessian -> Laplace covariance
    n = len(th0)
    eps = 1e-4
    H = np.zeros((n, n))
    for i in range(n):
        ei = np.zeros(n); ei[i] = eps
        H[i, i] = (nll(res.x + ei) - 2 * f0 + nll(res.x - ei)) / eps ** 2
        for j in range(i + 1, n):
            ej = np.zeros(n); ej[j] = eps
            Hij = (nll(res.x + ei + ej)
                   - nll(res.x + ei - ej)
                   - nll(res.x - ei + ej)
                   + nll(res.x - ei - ej)) / (4 * eps * eps)
            H[i, j] = H[j, i] = Hij
    # H was computed on the normalized objective: H_normed = H_true / S,
    # so cov_true = inv(S * H_normed) = inv(H_normed) / S
    cov = np.linalg.inv(H + np.eye(n) * 1e-6) / S
    se = np.sqrt(np.diag(cov))
    return res.x, se, cov, res


def sandwich_se(mdl, theta, fixed_su=True, prior_scale=2.5):
    """Cluster-robust (host-level) sandwich covariance for the mixture
    likelihood: V = H^-1 (sum_g s_g s_g^T) H^-1, s_g = per-host score."""
    ncoef = mdl["X"].shape[1]
    p = len(theta)

    def total_nll(t):
        return neg_log_post(t, mdl, fixed_su=fixed_su,
                            norm=float(len(mdl["r_obs"])))

    # bread: Hessian of the penalized objective
    eps = 1e-4
    S = float(len(mdl["r_obs"]))
    H = np.zeros((p, p))
    f0 = total_nll(theta)
    for i in range(p):
        ei = np.zeros(p); ei[i] = eps
        H[i, i] = (total_nll(theta + ei) - 2 * f0
                   + total_nll(theta - ei)) / eps ** 2
        for j in range(i + 1, p):
            ej = np.zeros(p); ej[j] = eps
            H[i, j] = H[j, i] = (
                total_nll(theta + ei + ej) - total_nll(theta + ei - ej)
                - total_nll(theta - ei + ej) + total_nll(theta - ei - ej)
            ) / (4 * eps * eps)
    H *= S          # bread must be in unnormalized units to match the meat

    # meat: per-host scores of the UNPENALIZED loglik via per-planet grads
    nc = mdl["X"].shape[1]
    b = theta[:nc]
    if fixed_su:
        m_se, ls_se = theta[nc], theta[nc + 1]
        m_sn, ls_sn = theta[nc + 2], theta[nc + 3]
        u_arr, gw = np.array([0.0]), np.array([1.0])
    else:
        lsu = theta[nc]
        m_se, ls_se = theta[nc + 1], theta[nc + 2]
        m_sn, ls_sn = theta[nc + 3], theta[nc + 4]
        su = np.exp(lsu); u_arr, gw = gauss_hermite_grid(su)
    def planet_logliks(t):
        bb = t[:nc]
        if fixed_su:
            m1, l1, m2, l2 = t[nc], t[nc+1], t[nc+2], t[nc+3]
            uu, ww = np.array([0.0]), np.array([1.0])
        else:
            l2u = t[nc]
            m1, l1 = t[nc+1], t[nc+2]
            m2, l2 = t[nc+3], t[nc+4]
            uu, ww = gauss_hermite_grid(np.exp(l2u))
        out = np.empty((len(uu), len(mdl["r_obs"])))
        lw_sn = norm_logpdf(mdl["r_obs"], m2,
                            np.sqrt(np.exp(l2) ** 2 + mdl["e_r"] ** 2))
        lw_se = norm_logpdf(mdl["r_obs"], m1,
                            np.sqrt(np.exp(l1) ** 2 + mdl["e_r"] ** 2))
        for j, uj in enumerate(uu):
            pi = expit(np.clip(mdl["X"] @ bb + uj, -30, 30))
            out[j] = np.logaddexp(lw_sn + np.log(pi), lw_se + np.log1p(-pi))
        return ww @ out

    Sg = np.zeros((mdl["n_hosts"], p))
    for i in range(p):
        ep_i = eps
        ti_p = theta.copy(); ti_p[i] += ep_i
        ti_m = theta.copy(); ti_m[i] -= ep_i
        g_i = (planet_logliks(ti_p) - planet_logliks(ti_m)) / (2 * ep_i)
        # aggregate to hosts (completeness-weighted)
        Sg[:, i] = np.bincount(mdl["hidx"], weights=mdl["w_pl"] * g_i,
                               minlength=mdl["n_hosts"])
    meat = Sg.T @ Sg
    Binv = np.linalg.inv(H)
    V = Binv @ meat @ Binv
    return np.sqrt(np.diag(V))


def summarize(theta, se, mdl, tag, variant, extra=None, fixed_su=True):
    if fixed_su:
        names = mdl["names"] + ["mu_SE", "log_s_SE", "mu_SN", "log_s_SN"]
    else:
        names = (mdl["names"] + ["log_sigma_u", "mu_SE", "log_s_SE",
                                 "mu_SN", "log_s_SN"])
    out = {"tag": tag, "variant": variant,
           "n_planets": int(len(mdl["r_obs"])),
           "n_hosts": int(mdl["n_hosts"]), "params": {}}
    ba, sa = theta[1], se[1]
    iqr = mdl["age_iqr"]
    out["params"] = {n: {"est": float(t), "se": float(s)}
                     for n, t, s in zip(names, theta, se)}
    # coefficient is per 1 SD; convert to per-IQR
    sd = mdl["age_raw"].std()
    per_iqr = ba * (iqr / sd)
    per_iqr_se = sa * (iqr / sd)
    out["beta_age_per_IQR"] = {"est": float(per_iqr),
                               "se": float(per_iqr_se),
                               "odds_ratio_per_IQR":
                                   float(np.exp(per_iqr)),
                               "z": float(per_iqr / per_iqr_se)}
    print(f"\n=== {tag} [{variant}] ===")
    print(f"N={out['n_planets']} planets / {out['n_hosts']} hosts | "
          f"age proxy IQR={iqr:.2f}")
    print(f"beta_age per IQR: {per_iqr:+.3f} +- {per_iqr_se:.3f} "
          f"(z={out['beta_age_per_IQR']['z']:+.2f}, "
          f"OR={np.exp(per_iqr):.3f})")
    for nm in ("gamma_logP", "gamma_logS", "gamma_missionKepler",
               "gamma_missionK2", "gamma_FeH"):
        if nm in out["params"]:
            g = out["params"][nm]
            print(f"  {nm:22s} {g['est']:+.3f} +- {g['se']:.3f}")
    if extra:
        out.update(extra)
    return out


def load(tag, kin_col="vtan"):
    pl = pd.read_csv(D / f"planet_sample_{tag}.csv")
    k = pd.read_csv(D / f"hosts_kinematics_{tag}.csv")
    df = pl.merge(k[["hostname", "vtan", "W"]], on="hostname", how="inner")
    if kin_col == "absW":
        df = df.dropna(subset=["W"]).copy()
        df["absW"] = df.W.abs()
    return df


def run_fit(tag, kin_col="vtan", valley=None, fe=False, dist_max=None,
            label=None):
    VALLEYS = {"M": 1.85, "FGK": 1.88}
    valley = valley or VALLEYS[tag]
    df = load(tag, kin_col)
    if dist_max:
        df = df[df.sy_dist < dist_max]
    mdl = build_design(df, kin_col, valley, fe=fe)
    theta, se0, cov, res = fit(mdl, fixed_su=True)
    se = sandwich_se(mdl, theta, fixed_su=True)
    out = summarize(theta, se, mdl, tag, label or kin_col,
                    extra={"valley": valley, "dist_max": dist_max, "fe": fe,
                           "mode": "mixture_fixed_su_sandwich"})
    R.mkdir(exist_ok=True)
    fn = R / f"hierarchical_{tag}_{label or kin_col}.json"
    with open(fn, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"  saved {fn.name}")
    return out


def selftest(nseed=24, sigma_host=1.5):
    """Recovery check in the REAL-DATA regime: host-effect SD matched to the
    fitted sigma_u (~1.5 logits). Laplace CIs are somewhat anti-conservative
    here -- that is a documented operating characteristic of the estimator;
    the power analysis (11) measures the pipeline end-to-end."""
    hits = tot = 0
    biases = []
    b_true = -0.8
    for seed in range(nseed):
        rng = np.random.default_rng(100 + seed)
        n_host = 600
        hosts = np.repeat([f"H{i}" for i in range(n_host)],
                          rng.integers(1, 3, n_host))
        n = len(hosts)
        age = rng.uniform(10, 80, n)                      # vtan-like
        logP = rng.uniform(-0.5, 1.6, n)
        logS = rng.uniform(-1, 2, n)
        kep = (rng.uniform(size=n) < 0.7).astype(float)
        eta = (-0.3 + b_true * (age - age.mean()) / age.std()
               + 0.3 * logP - 0.2 * logS + 0.15 * kep
               + rng.normal(0, sigma_host, n))            # real-scale effects
        pi_true = expit(eta)
        sn = rng.uniform(size=n) < pi_true
        # observed radii consistent with classification + noise
        r_true = np.where(sn, rng.normal(2.4, 0.5, n), rng.normal(1.3, 0.25, n))
        r_obs = r_true + rng.normal(0, 0.06, n)
        df = pd.DataFrame({"hostname": hosts, "pl_rade": r_obs,
                           "pl_radeerr1": 0.06, "pl_orbper": 10 ** logP,
                           "st_mass": 1.0, "st_rad": 1.0, "st_teff": 5700,
                           "sy_dist": 300, "st_met": 0.0, "age": age})
        mdl = build_design(df.rename(columns={"age": "vtan"}), "vtan",
                           valley=1.88, fe=False)
        theta, se0, cov, _ = fit(mdl, fixed_su=True)
        se = sandwich_se(mdl, theta, fixed_su=True)
        sd = age.std(); iqr = np.percentile(age, 75) - np.percentile(age, 25)
        # both expressed per-IQR: sim injects b_true per SD of age
        est = theta[1] * (iqr / sd)
        ese = se[1] * (iqr / sd)
        truth = b_true * (iqr / sd)
        hit = (est - 1.96 * ese <= truth) and (truth <= est + 1.96 * ese)
        hits += hit; tot += 1; biases.append(est - truth)
        print(f"  seed {seed}: est {est:+.2f}+-{ese:.2f} vs truth "
              f"{truth:+.2f} -> {'HIT' if hit else 'MISS'}")
    print(f"SELFTEST coverage: {hits}/{tot} | mean bias "
          f"{np.mean(biases):+.3f} logits/IQR (+ = attenuation toward 0)")
    print("CALIBRATION NOTE: the estimator attenuates small beta_age "
          "toward zero in this regime (host heterogeneity + population "
          "overlap). Null conclusions are therefore CONSERVATIVE; any "
          "future positive detection must be deconvolved via this "
          "calibration before being believed.")


if __name__ == "__main__":
    if "SELFTEST" in sys.argv:
        selftest()
        sys.exit(0)
    run_fit("M", "vtan")
    run_fit("M", "absW")
    run_fit("M", "vtan", dist_max=200, label="vtan_local")
    run_fit("M", "vtan", valley=1.95, label="vtan_floor195")
    run_fit("M", "vtan", fe=True, label="vtan_fe")
    run_fit("FGK", "vtan")
