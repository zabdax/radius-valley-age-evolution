"""
colab_highrep_gate.py — HIGHER-REPLICATE Monte Carlo replication on Colab.

Standalone script. Its ONLY job: re-run the project's already-decided
calibrated-gate Monte Carlo (and, budget permitting, the 15-cell power
grid) at a higher replicate count than the local runs (120/150), and dump
RAW replicate-level results for import back into the main project.

It adds NO new methodology, NO new diagnostics, NO new model structure,
and it does NOT touch the rest of the pipeline or PAPER_PUBLICATION.pdf.

MODEL SPECIFICATION (unchanged, per the Issue C investigation):
  - working-independence mixture likelihood, host random intercept FIXED
    at sigma_u = 0  (fit(fixed_su=True, fix_comp=True));
  - host-level cluster-robust sandwich SEs (per-planet scores aggregated
    to per-host scores before the meat matrix);
  - host-bootstrap injections: whole hosts resampled with replacement,
    ONE shared host effect per drawn host, sigma_u = 1.5;
  - verbatim copies of code/09_hierarchical_model.py (model),
    code/13_calibration_fgk.py (gate), code/11_power_analysis.py (grid).

FIDELITY (replicate-count increase only):
  - Gate: same seed 20260826, same arm order (null -> signal), same call
    order, same constants (B_INJ=-0.14/SD, target SN frac 0.408,
    gamma_logP=2.3, sigma_u=1.5, class conditionals N(2.10,0.35)/
    N(1.20,0.18)), same sample (FGK Kepler-only, trim 95th pct ->
    N=2,135, ESS~1,645). With the same numpy version the FIRST 120
    replicates of each arm are bit-identical to the local run; this
    script continues the same stream to GATE_REPS.
  - Grid: same 15 cells, same per-cell seeds (eff_seed*1000+scale*100),
    same per-cell rng order (alpha calibration first -- it CONSUMES rng
    draws -- then replicates), same effect sizes (0 / -0.14 / -0.28 per
    SD), same VALLEY=1.85 M-sample constants. Only the replicate count
    changes (GRID_REPS per cell).

INPUTS: upload colab_gate_data.zip (~0.4 MB, the 5 frozen CSVs) to
  /content/, or set COLAB_GATE_DATA to a folder with the 5 CSVs.
  The script validates the frozen sample and aborts on mismatch.

OUTPUTS (in /content/colab_gate_output/ + one zip):
  task2_calibration_kepler_only_trim0.95_rep<GATE_REPS>.json  (13 schema)
  task2_gate_rep<GATE_REPS>_replicates_raw.csv                (raw rows)
  power_analysis_rep<GRID_REPS>.json                          (11 schema)
  power_grid_rep<GRID_REPS>_replicates_raw.csv                (raw rows)
  checkpoints/  (resumable progress + rng states; operational only)

USAGE (Colab): Runtime > CPU; upload colab_gate_data.zip; then
  !python colab_highrep_gate.py        (or paste into one cell and run)
  Optional: --smoke (tiny plumbing test). If the session dies, re-upload
  colab_gate_output.zip (it contains the checkpoints), rerun -- it
  RESUMES from the last checkpoint with the exact rng stream.
"""
import os
import sys
import json
import time
import shutil
import zipfile
from pathlib import Path

# ------------------------------------------------------------------
# CONFIG (edit here; defaults follow the tasking)
# ------------------------------------------------------------------
GATE_REPS = 1000        # per arm (null-injected, signal-injected); >=1000
RUN_GRID = True         # run the 15-cell power grid after the gate
GRID_REPS = 500         # per cell (was 120 null / 150 effect locally)
GRID_BUDGET_HOURS = 5.0 # wall-clock budget for the GRID ONLY; the gate is
                        # ALWAYS completed first and never shortened.
                        # Cells skipped by the budget are recorded in the
                        # manifest; finished cells are fully usable.
CHECKPOINT_EVERY = 25   # reps between checkpoint writes
SEED_GATE = 20260826    # verbatim from 13_calibration_fgk.py
SMOKE = "--smoke" in sys.argv
if SMOKE:               # plumbing test only; writes to _smoke outputs
    GATE_REPS = 6
    GRID_REPS = 3
    GRID_BUDGET_HOURS = 0.2
    CHECKPOINT_EVERY = 2

# ------------------------------------------------------------------
# Dependencies (Colab has numpy/pandas/scipy; install if missing)
# ------------------------------------------------------------------
def _ensure_packages():
    import importlib
    import subprocess
    for pkg in ("numpy", "pandas", "scipy"):
        try:
            importlib.import_module(pkg)
        except Exception:
            print(f"installing {pkg} ...")
            subprocess.run([sys.executable, "-m", "pip", "install", "-q",
                            pkg], check=False)

_ensure_packages()
import numpy as np
import pandas as pd

# ------------------------------------------------------------------
# Paths: locate the 5 frozen data CSVs
# ------------------------------------------------------------------
REQUIRED_FILES = [
    "planet_sample_FGK.csv",
    "hosts_kinematics_FGK.csv",
    "completeness_weights_FGK.csv",
    "planet_sample_M.csv",
    "hosts_kinematics_M.csv",
]

def _script_dir():
    try:
        return Path(__file__).resolve().parent
    except Exception:
        return Path(".").resolve()

def find_data_dir():
    for zc in [Path("/content/colab_gate_data.zip"),
               _script_dir() / "colab_gate_data.zip"]:
        if zc.exists():
            dest = zc.parent / "colab_gate_data"
            dest.mkdir(exist_ok=True)
            with zipfile.ZipFile(zc) as zf:
                zf.extractall(dest)
            print(f"extracted {zc.name} -> {dest}")
            return dest
    env = os.environ.get("COLAB_GATE_DATA", "")
    cands = ([Path(env)] if env else []) + [
        Path("/content/colab_gate_data"),
        Path("/content"),
        Path("/content/drive/MyDrive/colab_gate_data"),
        Path("/content/drive/MyDrive"),
        _script_dir() / "data",          # legacy local layout
        _script_dir(),
        _script_dir().resolve().parents[1] / "data",   # repo layout (script in colab/)
        _script_dir().resolve().parents[1],
    ]
    for c in cands:
        if all((c / f).exists() for f in REQUIRED_FILES):
            return c
    raise SystemExit(
        "ERROR: could not find the 5 required data CSVs. Upload "
        "colab_gate_data.zip to /content/ or set COLAB_GATE_DATA to the "
        f"folder containing: {REQUIRED_FILES}")

if SMOKE:
    OUT_DIR = _script_dir() / "colab_gate_output_smoke"
else:
    OUT_DIR = Path("/content/colab_gate_output")
CKPT_DIR = OUT_DIR / "checkpoints"
OUT_DIR.mkdir(parents=True, exist_ok=True)
CKPT_DIR.mkdir(parents=True, exist_ok=True)
DATA_DIR = find_data_dir()
print(f"data dir: {DATA_DIR}")
print(f"output dir: {OUT_DIR}")

# ------------------------------------------------------------------
# MODEL CORE (verbatim from code/09_hierarchical_model.py)
# ------------------------------------------------------------------
from scipy.special import expit, logsumexp
from scipy.optimize import minimize

GH_NODES, GH_W = np.polynomial.hermite_e.hermegauss(24)


def gauss_hermite_grid(sigma):
    # probabilists' convention: hermegauss nodes satisfy
    # sum(w_j f(x_j))/sqrt(2 pi) ~= E[f(Z)], Z~N(0,1); scale nodes by sigma.
    u = sigma * GH_NODES
    w = GH_W / np.sqrt(2.0 * np.pi)
    return u, w


def build_design(df, age_col, valley, fe=False, weights=None):
    """Assemble data quantities and design matrix (mixture likelihood v2:
    the outcome is the OBSERVED RADIUS itself, modeled as a two-component
    Gaussian mixture whose mixing fraction pi(age, covariates) is the
    science target)."""
    d = df.copy()
    err = df.get("pl_radeerr1")
    e = np.asarray(err.abs().values if err is not None
                   else np.full(len(d), np.nan), float)
    # physical clip: cap relative error at 30%, floor at 1%
    r_abs = np.abs(d.pl_rade.values) + 1e-6
    e = np.where(np.isfinite(e) & (e > 0.01), e, 0.05 * r_abs)
    e = np.clip(e, 0.01 * r_abs, 0.30 * r_abs)

    # physical insolation (S in solar constants): a[AU]=(M*)^(1/3)(P/yr)^(2/3)
    mstar = d.st_mass.fillna(0.5).values
    mstar = np.clip(mstar, 0.08, 1.6)
    p_yr = d.pl_orbper.values / 365.25
    a_au = np.cbrt(mstar) * p_yr ** (2.0 / 3.0)
    rstar = d.st_rad.fillna(0.6).values if "st_rad" in d else 0.6
    tstar = d.st_teff.fillna(4000).values if "st_teff" in d else 4000
    lum = np.clip(rstar, 0.1, 4) ** 2 * (np.clip(tstar, 2500, 7500) / 5772.0) ** 4
    s_in = lum / a_au ** 2

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
    """Negative log posterior. norm divides the FULL objective to keep
    gradients O(1). fixed_su=True: no host random intercept (sigma_u=0);
    use with the cluster-robust sandwich SE."""
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
    """Midpoint between fitted/initial components; ordering constraint."""
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
    unidentifiable in a minority of replicates; fixed components remove
    it. MAP + central-difference Hessian -> Laplace covariance."""
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
    # H was computed on the normalized objective: cov_true = inv(H)/S
    cov = np.linalg.inv(H + np.eye(n) * 1e-6) / S
    se = np.sqrt(np.diag(cov))
    return res.x, se, cov, res


def sandwich_se(mdl, theta, fixed_su=True, prior_scale=2.5):
    """Cluster-robust (host-level) sandwich covariance for the mixture
    likelihood: V = H^-1 (sum_g s_g s_g^T) H^-1, s_g = per-host score.
    Per-planet scores are aggregated to per-host scores (completeness-
    weighted) BEFORE the meat matrix -- this is what makes the reported
    uncertainties robust to planets sharing a host."""
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

# ------------------------------------------------------------------
# GATE GENERATIVE PROCESS (verbatim from code/13_calibration_fgk.py)
# ------------------------------------------------------------------
B_INJ = -0.14          # SWEET-Cat-equivalent, per SD of v_tan
GATE_VALLEY = 1.88
GATE_TARGET_FRAC = 0.408   # hardcoded in 13's arm loop


def load_weighted_fgk():
    """13.load_weighted('FGK') + 13.main's kepler_only + trim 0.95."""
    pl = pd.read_csv(DATA_DIR / "planet_sample_FGK.csv")
    k = pd.read_csv(DATA_DIR / "hosts_kinematics_FGK.csv")
    w = pd.read_csv(DATA_DIR / "completeness_weights_FGK.csv")
    df = pl.merge(k[["hostname", "vtan", "W"]], on="hostname", how="inner")
    df = df.merge(w[["pl_name", "w", "p_det"]], on="pl_name", how="left")
    df["w"] = df.w.fillna(1.0)
    df = df[df.hostname.str.startswith("Kepler")].copy()
    cap = np.quantile(df.w, 0.95)
    df["w"] = np.minimum(df.w, cap)
    df = df.rename(columns={"w": "wgt"}).reset_index(drop=True)
    ess = float(df.wgt.sum() ** 2 / (df.wgt ** 2).sum())
    return df, ess


def calibrate_alpha_gate(df, target_frac, rng):
    """Verbatim 13.calibrate_alpha (bisection; the rng arg is unused
    there -- calibration is deterministic)."""
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


def synth_gate(base, beta_inj, alpha, rng):
    """Verbatim 13.synth: host bootstrap (size = number of hosts), ONE
    shared host effect per drawn host (sigma_u = 1.5), gamma_logP = 2.3
    hardcoded (measured gamma_logP from the unweighted FGK fit)."""
    picked = rng.choice(base.hostname.unique(),
                        size=base.hostname.nunique(), replace=True)
    pos = {h: i for i, h in enumerate(picked)}
    parts = [base[base.hostname == h] for h in picked]
    d = pd.concat(parts, ignore_index=True)
    zv = (d.vtan.values - base.vtan.mean()) / base.vtan.std()
    zp = (np.log10(d.pl_orbper.values) - np.log10(base.pl_orbper).mean()) \
        / np.log10(base.pl_orbper).std()
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


def fit_real_gate(df):
    """13.main REAL weighted fit: mixture fit + host-cluster sandwich SE."""
    mdl = build_design(df, "vtan", GATE_VALLEY, fe=False,
                       weights=df.wgt.values)
    theta, _, _, _ = fit(mdl, fixed_su=True, fix_comp=True)
    se = sandwich_se(mdl, theta, fixed_su=True)
    sd = df.vtan.std()
    iqr = df.vtan.quantile(0.75) - df.vtan.quantile(0.25)
    return theta[1] * (iqr / sd), se[1] * (iqr / sd)

# ------------------------------------------------------------------
# GRID GENERATIVE PROCESS (verbatim from code/11_power_analysis.py)
# ------------------------------------------------------------------
GAMMA_LOGP, GAMMA_LOGS = 1.613, -0.032
GAMMA_KEP, GAMMA_K2 = -1.594, 0.893
SIGMA_U = 1.5
MU_SE, S_SE = 1.20, 0.18
MU_SN, S_SN = 2.10, 0.35
GRID_VALLEY = 1.85


def load_base_m():
    """Verbatim 11.load_base: M sample + kinematics, e_r, mission fam."""
    pl = pd.read_csv(DATA_DIR / "planet_sample_M.csv")
    k = pd.read_csv(DATA_DIR / "hosts_kinematics_M.csv")
    df = pl.merge(k[["hostname", "vtan", "W"]], on="hostname",
                  how="inner").reset_index(drop=True)
    err = df.pl_radeerr1.abs()
    e = np.where(np.isfinite(err) & (err > 0.01), err, 0.05 * df.pl_rade)
    df["e_r"] = e
    df["mission"] = df.hostname.str.extract(
        r"^(TOI|Kepler|KIC|K2|HIP|HD|GJ|Wolf|LHS|TRAPPIST|EPIC)",
        expand=False).fillna("other")
    df["fam"] = df.mission.map(
        lambda s: s if s in ("Kepler", "K2", "TOI") else "other")
    return df


def zscore_within(df, col):
    v = df[col].values.astype(float)
    return (v - v.mean()) / (v.std() + 1e-12)


def synth_grid(base, scale, beta_inj, alpha, rng):
    """Verbatim 11.synth: bootstrap whole hosts at `scale`x, z-score
    covariates WITHIN the synthetic sample, shared host effect per drawn
    host, reality-matched class conditionals + per-planet radius errors."""
    hosts = base.hostname.unique()
    n_hosts = max(1, int(round(scale * len(hosts))))
    picked = rng.choice(hosts, size=n_hosts, replace=True)
    idx = np.concatenate([np.flatnonzero(base.hostname.values == h)
                          for h in picked])
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


def calibrate_alpha_grid(base, beta_inj, target_frac, rng, mission_offset=0.0):
    """Verbatim 11.calibrate_alpha. NOTE: CONSUMES rng draws (rng.choice
    per bisection iteration) -- call order matters for stream fidelity."""
    lo, hi = -6.0, 6.0
    for _ in range(30):
        mid = 0.5 * (lo + hi)
        hosts = base.hostname.unique()
        picked = rng.choice(hosts, size=len(hosts), replace=True)
        first_idx = [np.flatnonzero(base.hostname.values == h)[0]
                     for h in picked]
        d1 = base.iloc[first_idx]
        zv = (d1.vtan.values - base.vtan.mean()) / base.vtan.std()
        zp = (np.log10(d1.pl_orbper.values)
              - np.log10(base.pl_orbper).mean()) \
            / np.log10(base.pl_orbper).std()
        eta = (mid + mission_offset + beta_inj * zv + GAMMA_LOGP * zp)
        frac = np.mean(1 / (1 + np.exp(-eta)))
        if frac < target_frac:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)

# ------------------------------------------------------------------
# OPERATIONAL HELPERS (checkpoints, raw CSV, validation, estimates)
# -- none of this is methodology; it only enables resume + raw output --
# ------------------------------------------------------------------
import csv as _csv


def save_ckpt(name, obj):
    tmp = CKPT_DIR / (name + ".tmp")
    with open(tmp, "w") as fh:
        json.dump(obj, fh)
    tmp.replace(CKPT_DIR / f"{name}.json")


def load_ckpt(name):
    p = CKPT_DIR / f"{name}.json"
    if p.exists():
        with open(p) as fh:
            return json.load(fh)
    return None


def append_rows(csv_path, rows, header):
    new = not csv_path.exists()
    with open(csv_path, "a", newline="") as fh:
        wr = _csv.writer(fh)
        if new:
            wr.writerow(header)
        wr.writerows(rows)


def rewrite_csv_without_reps(csv_path, key_col, key_val, keep_reps):
    """Drop raw rows beyond the checkpoint for one series (crash-between-
    writes guard). Operational only."""
    if not csv_path.exists():
        return
    # keep_default_na=False: our raw CSVs contain the literal string
    # "null" (the gate arm name) -- default NA parsing would turn it
    # into NaN and silently empty the selection.
    df = pd.read_csv(csv_path, keep_default_na=False)
    mask = (df[key_col] == key_val) & (df["rep"] >= keep_reps)
    if mask.any():
        df = df[~mask]
        df.to_csv(csv_path, index=False)


def validate_gate_df(df, ess):
    if len(df) != 2135:
        raise SystemExit(f"ABORT: FGK Kepler-only trimmed N={len(df)}, "
                         "expected 2135 -- the uploaded data files do not "
                         "match the frozen project sample.")
    if abs(ess - 1645) > 3:
        raise SystemExit(f"ABORT: ESS={ess:.0f}, expected ~1645 -- wrong "
                         "completeness_weights file?")
    frac = float(np.mean(df.pl_rade > GATE_VALLEY))
    print(f"GATE sample validated: N={len(df)}, ESS={ess:.0f}, "
          f"hosts={df.hostname.nunique()}, SN frac={frac:.3f}")


def estimate_runtime(df_gate, base_m):
    """Timing warm-up on throwaway replicates (separate rng; does NOT
    touch the gate stream). Prints estimates before committing."""
    print("\n=== RUNTIME ESTIMATE (warm-up) ===")
    rng = np.random.default_rng(987654321)
    a = calibrate_alpha_gate(df_gate, GATE_TARGET_FRAC, rng)
    t0 = time.time()
    for _ in range(2):
        ds = synth_gate(df_gate, 0.0, a, rng)
        mdl = build_design(ds, "vtan", GATE_VALLEY, fe=False,
                           weights=ds.wgt.values)
        th, _, _, _ = fit(mdl, fixed_su=True, fix_comp=True)
    s_gate = (time.time() - t0) / 2.0
    print(f"measured: {s_gate:.1f} s per gate replicate "
          f"(N={len(df_gate)}: synth + mixture fit)")
    gate_h = 2 * GATE_REPS * s_gate / 3600.0
    print(f"gate: 2 arms x {GATE_REPS} reps -> ~{gate_h:.1f} h")
    total_h = gate_h
    if RUN_GRID:
        d1 = synth_grid(base_m, 1.0, 0.0, 0.0, rng)
        mdl = build_design(d1, "vtan", GRID_VALLEY, fe=False)
        t1 = time.time()
        th, _, _, _ = fit(mdl, fixed_su=True)
        se = sandwich_se(mdl, th, fixed_su=True)
        s_grid_unit = (time.time() - t1) / max(len(d1), 1)
        cells_h = 0.0
        for scale in (0.5, 1.0, 2.0, 4.0, 8.0):
            n_est = scale * len(base_m)
            cells_h += 3 * GRID_REPS * n_est * s_grid_unit / 3600.0
        print(f"measured: {s_grid_unit * 1000:.1f} ms per planet for a "
              "grid replicate (synth + fit + sandwich SE)")
        print(f"grid: 15 cells x {GRID_REPS} reps -> ~{cells_h:.1f} h "
              f"(budget {GRID_BUDGET_HOURS} h; unfinished cells are "
              "skipped, finished cells are complete and usable)")
        total_h += min(cells_h, GRID_BUDGET_HOURS)
    print(f"TOTAL estimated: ~{total_h:.1f} h this session "
          f"(gate always completes; everything is resumable)")

# ------------------------------------------------------------------
# GATE RUNNER (13's protocol; replicate count is the only change)
# ------------------------------------------------------------------
def run_gate(df, ess):
    gate_csv = OUT_DIR / f"task2_gate_rep{GATE_REPS}_replicates_raw.csv"
    summary_path = (OUT_DIR /
                    f"task2_calibration_kepler_only_trim0.95_rep{GATE_REPS}.json")
    gate_header = ["arm", "rep", "beta_inj", "alpha", "beta_hat_per_iqr",
                   "n_planets"]

    print("\n=== GATE: REAL weighted fit (FGK Kepler-only, trim 0.95) ===")
    b_real, se_real = fit_real_gate(df)
    sd = df.vtan.std()
    iqr = df.vtan.quantile(0.75) - df.vtan.quantile(0.25)
    print(f"REAL weighted FGK: beta_age/IQR = {b_real:+.3f} +- {se_real:.3f}")

    res = {"beta_real": float(b_real), "se_real": float(se_real),
           "reps": GATE_REPS,
           "n_planets": int(len(df)), "ess": float(ess),
           "model_spec": ("working-independence mixture (fixed_su=True) + "
                          "host-cluster sandwich SEs; unchanged per Issue C"),
           "seed": SEED_GATE, "generated_by": "colab_highrep_gate.py"}

    ck = load_ckpt("gate")
    if ck and ck.get("reps") == GATE_REPS:
        rng = restore_rng(ck["rng_state"])
        done = ck["done"]
        print(f"resuming gate: null={done.get('null', 0)}/{GATE_REPS}, "
              f"signal={done.get('signal', 0)}/{GATE_REPS}")
    elif ck:
        raise SystemExit("ABORT: existing gate checkpoint was made with a "
                         "different GATE_REPS; move/delete the old output "
                         "folder first.")
    else:
        rng = np.random.default_rng(SEED_GATE)   # verbatim 13 seed
        done = {"null": 0, "signal": 0}

    for name, binj in [("null", 0.0), ("signal", B_INJ)]:
        a = calibrate_alpha_gate(df, GATE_TARGET_FRAC, rng)  # no rng draws
        start = done.get(name, 0)
        if start >= GATE_REPS:
            print(f"[{name}] already complete ({GATE_REPS} reps)")
            continue
        print(f"--- arm [{name}] beta_inj={binj}: reps {start}..{GATE_REPS-1}")
        rewrite_csv_without_reps(gate_csv, "arm", name, start)
        for rep in range(start, GATE_REPS):
            ds = synth_gate(df, binj, a, rng)
            mdl = build_design(ds, "vtan", GATE_VALLEY, fe=False,
                               weights=ds.wgt.values)
            th, _, _, _ = fit(mdl, fixed_su=True, fix_comp=True)
            b_iqr = th[1] * (iqr / sd)
            append_rows(gate_csv, [[name, rep, binj, a, b_iqr, len(ds)]],
                        gate_header)
            done[name] = rep + 1
            if (rep + 1) % CHECKPOINT_EVERY == 0 or (rep + 1) == GATE_REPS:
                save_ckpt("gate", {"reps": GATE_REPS, "done": done,
                                   "next_rep": rep + 1, "arm": name,
                                   "rng_state": rng.bit_generator.state})
            if (rep + 1) % 50 == 0:
                print(f"  [{name}] rep {rep + 1}/{GATE_REPS}", flush=True)
        save_ckpt("gate", {"reps": GATE_REPS, "done": done,
                           "next_rep": done[name], "arm": name,
                           "rng_state": rng.bit_generator.state})

    # ---- summary, schema-superset of 13's task2 JSON ----
    # keep_default_na=False: the arm column literally contains "null".
    raw = pd.read_csv(gate_csv, keep_default_na=False)
    res["reps_completed"] = {a: int((raw.arm == a).sum())
                             for a in ("null", "signal")}
    for name in ("null", "signal"):
        bets = raw[raw.arm == name].sort_values("rep").beta_hat_per_iqr.values
        res[name] = {"median": float(np.median(bets)),
                     "q025": float(np.quantile(bets, .025)),
                     "q975": float(np.quantile(bets, .975)),
                     "pct_rank_of_real": float(np.mean(bets <= b_real))}
        print(f"[{name}-injected] median {res[name]['median']:+.3f} "
              f"[{res[name]['q025']:+.3f},{res[name]['q975']:+.3f}] | "
              f"real percentile {res[name]['pct_rank_of_real']:.3f}")

    # ---- PASS/FAIL per protocol (verbatim 13 logic) ----
    p_null_upper = 1.0 - res["null"]["pct_rank_of_real"]
    inside_sig = (res["signal"]["q025"] <= b_real <= res["signal"]["q975"])
    res["verdict"] = {
        "distinguishable_from_null": bool(p_null_upper <= 0.05),
        "consistent_with_signal": bool(inside_sig)}
    res["verdict"]["PASS"] = bool(
        res["verdict"]["distinguishable_from_null"]
        and res["verdict"]["consistent_with_signal"])
    print("\n=== TASK 2 VERDICT (high-replicate) ===")
    print(f"P(null-injected >= real): {p_null_upper:.4f}")
    print(f"real inside 95% of signal-injected distribution: {inside_sig}")
    print("VERDICT:", "PASS" if res["verdict"]["PASS"] else
          ("FAIL(inconclusive)" if not
           res["verdict"]["distinguishable_from_null"] else
           "FAIL(null-like)"))
    with open(summary_path, "w") as fh:
        json.dump(res, fh, indent=1)
    print(f"saved {summary_path.name}")
    return time.time()


def restore_rng(state):
    """Rebuild a Generator from a JSON-round-tripped PCG64 state."""
    r = np.random.default_rng()
    r.bit_generator.state = state
    return r


# ------------------------------------------------------------------
# GRID RUNNER (11's protocol; replicate count is the only change)
# ------------------------------------------------------------------
def run_grid(base, t_grid_start):
    grid_raw = OUT_DIR / f"power_grid_rep{GRID_REPS}_replicates_raw.csv"
    grid_json = OUT_DIR / f"power_analysis_rep{GRID_REPS}.json"
    grid_header = ["cell", "scale", "effect", "beta_inj", "rep",
                   "beta_hat_per_sd", "se_per_sd", "detect", "cover",
                   "n_planets"]

    real_frac = float((base.pl_rade > GRID_VALLEY).mean())
    calib_frac = min(0.495, real_frac / 0.90)
    f_kep = float((base.fam == "Kepler").mean())
    f_k2 = float((base.fam == "K2").mean())
    mission_offset = GAMMA_KEP * f_kep + GAMMA_K2 * f_k2
    print(f"\n=== GRID: real M SN fraction {real_frac:.3f}; "
          f"{base.hostname.nunique()} hosts, {len(base)} planets ===")
    print(f"mission offset added to calibration: {mission_offset:+.3f} logits")

    eff_seed = {"null": 0, "SWEET": 1, "strong": 2}
    cells = []
    for scale, label in [(0.5, "0.5x"), (1.0, "1x"), (2.0, "2x"),
                         (4.0, "4x"), (8.0, "8x")]:
        for eff_name, beta_inj in [("null", 0.0), ("SWEET", -0.14),
                                   ("strong", -0.28)]:
            seed = eff_seed[eff_name] * 1000 + int(round(scale * 100))
            cells.append((scale, label, eff_name, beta_inj, seed))

    deadline = t_grid_start + GRID_BUDGET_HOURS * 3600.0
    skipped = []
    for scale, label, eff_name, beta_inj, seed in cells:
        ckey = f"grid_{label}_{eff_name}"
        cval = f"{label}|{eff_name}"
        if time.time() > deadline:
            skipped.append(cval)
            print(f"[{label}|{eff_name}] SKIPPED (grid budget exhausted)")
            continue
        ck = load_ckpt(ckey)
        if ck and ck.get("reps") == GRID_REPS:
            if ck.get("next_rep", 0) >= GRID_REPS:
                print(f"[{label}|{eff_name}] already complete")
                continue
            rng = restore_rng(ck["rng_state"])
            alpha = ck["alpha"]
            start = ck["next_rep"]
        else:
            rng = np.random.default_rng(seed)   # verbatim 11 per-cell seed
            alpha = calibrate_alpha_grid(base, beta_inj, calib_frac, rng,
                                         mission_offset=mission_offset)
            start = 0
        rewrite_csv_without_reps(grid_raw, "cell", cval, start)
        for rep in range(start, GRID_REPS):
            d = synth_grid(base, scale, beta_inj, alpha, rng)
            mdl = build_design(d, "vtan", GRID_VALLEY, fe=False)
            th, _, _, _ = fit(mdl, fixed_su=True)
            se = sandwich_se(mdl, th, fixed_su=True)
            b, s = th[1], se[1]
            excl0 = (b - 1.96 * s) * (b + 1.96 * s) > 0
            covers = (b - 1.96 * s <= beta_inj) and (beta_inj <= b + 1.96 * s)
            detect = bool(excl0 if beta_inj == 0
                          else (excl0 and (np.sign(b) == np.sign(beta_inj))))
            append_rows(grid_raw, [[cval, label, eff_name, beta_inj, rep,
                                    b, s, int(detect), int(covers), len(d)]],
                        grid_header)
            if (rep + 1) % CHECKPOINT_EVERY == 0 or (rep + 1) == GRID_REPS:
                save_ckpt(ckey, {"reps": GRID_REPS, "next_rep": rep + 1,
                                 "alpha": alpha,
                                 "rng_state": rng.bit_generator.state})
            if (rep + 1) % 50 == 0:
                print(f"  [{label}|{eff_name}] rep {rep + 1}/{GRID_REPS}",
                      flush=True)
        print(f"[{label:>4} | {eff_name:>6}] done ({GRID_REPS} reps, "
              f"alpha={alpha:+.3f})", flush=True)

    # ---- summary, schema-identical to 11's power_analysis.json ----
    if not grid_raw.exists():
        print("grid: no rows produced")
        return skipped
    raw = pd.read_csv(grid_raw, keep_default_na=False)
    results = []
    for _scale, label, eff_name, beta_inj, _seed in cells:
        sub = raw[raw.cell == f"{label}|{eff_name}"].sort_values("rep")
        if len(sub) == 0:
            continue
        results.append({
            "scale": label, "effect": eff_name,
            "beta_inj": float(beta_inj), "reps": int(len(sub)),
            "detect_rate": float(sub.detect.mean()),
            "coverage": float(sub.cover.mean()),
            "median_beta_hat": float(np.median(sub.beta_hat_per_sd.values))})
        with open(grid_json, "w") as fh:
            json.dump(results, fh, indent=1)
    print(f"saved {grid_json.name} ({len(results)} cells)")
    return skipped

# ------------------------------------------------------------------
# MAIN
# ------------------------------------------------------------------
def main():
    t0 = time.time()
    print("=" * 68)
    print("HIGH-REPLICATE CALIBRATED-GATE MC -- spec unchanged (Issue C):")
    print("  working-independence mixture (fixed_su=True); host-cluster")
    print("  sandwich SEs; host-bootstrap injections (shared sigma_u=1.5)")
    print(f"  gate: {GATE_REPS}/arm | grid: {GRID_REPS}/cell "
          f"(RUN_GRID={RUN_GRID})")
    print("=" * 68)
    df_gate, ess = load_weighted_fgk()
    validate_gate_df(df_gate, ess)
    base_m = load_base_m()
    if not (400 <= len(base_m) <= 430):
        raise SystemExit(f"ABORT: M sample N={len(base_m)}, expected ~417 "
                         "-- wrong data files?")
    print(f"GRID base validated: N={len(base_m)}, "
          f"hosts={base_m.hostname.nunique()}")
    estimate_runtime(df_gate, base_m)

    t_grid_start = run_gate(df_gate, ess)
    skipped = []
    if RUN_GRID:
        skipped = run_grid(base_m, t_grid_start)

    manifest = {
        "generated_by": "colab_highrep_gate.py",
        "gate_reps_per_arm": GATE_REPS,
        "grid_reps_per_cell": GRID_REPS,
        "grid_skipped_cells": skipped,
        "elapsed_hours": (time.time() - t0) / 3600.0,
        "outputs": sorted(p.name for p in OUT_DIR.iterdir() if p.is_file()),
        "model_spec": ("unchanged per Issue C: fixed_su mixture + "
                       "host-cluster sandwich + host-bootstrap injections "
                       "(shared sigma_u=1.5)"),
        "seed_gate": SEED_GATE}
    with open(OUT_DIR / "manifest.json", "w") as fh:
        json.dump(manifest, fh, indent=1)
    shutil.make_archive(str(OUT_DIR.parent / "colab_gate_output"), "zip",
                        OUT_DIR)
    print(f"\nDONE in {(time.time() - t0) / 3600.0:.2f} h. Download "
          f"{OUT_DIR.parent / 'colab_gate_output.zip'}")
    print("Import: copy the two *_replicates_raw.csv and the two summary "
          "JSONs into mdwarf_radius_valley/results/ (JSON schemas are "
          "supersets of the existing task2/power_analysis formats).")


if __name__ == "__main__":
    main()