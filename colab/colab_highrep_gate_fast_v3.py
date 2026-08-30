#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
colab_highrep_gate_fast_v3.py

Faster, Drive-persistent runner for the existing calibrated-gate Monte Carlo.

SCIENTIFIC SPECIFICATION (unchanged):
  - data, seeds, synthetic-generation order, effect sizes,
    model objective, optimizer, fixed_su=True, component bounds,
    host-bootstrap design, cluster-robust sandwich construction, detection rules.
  - synthetic datasets are generated in the original serial RNG order.
  - only the deterministic fitting work is parallelized.
  - synthetic gate fits skip Hessian/covariance calculation.
  - grid sandwich meat uses the exact analytic score of the SAME likelihood,
    while the bread/Hessian remains the same finite-difference calculation,
    computed with the exact same arithmetic as before (bit-identical).

PERSISTENCE:
  Outputs/checkpoints are written to Google Drive. The script verifies that
  Drive is ALREADY mounted and writable before doing expensive work, and
  resumes from checkpoints. It does NOT call drive.mount() itself, so it
  never re-triggers Colab's mount/auth flow — mount Drive yourself first
  (e.g. via the Colab file browser's "Mount Drive" button, or
  `from google.colab import drive; drive.mount('/content/drive')` in an
  earlier cell).

REQUIRED FILES IN /content:
  colab_highrep_gate.py
  colab_gate_data.zip
  this file

RUN:
  !python /content/colab_highrep_gate_fast_v3.py

OPTIONAL ENVIRONMENT VARIABLES:
  MC_WORKERS=8            # number of parallel workers
  MC_BATCH=8              # number of replicates per checkpoint batch
  GRID_REPS=500           # number of grid replicates per cell
  GATE_REPS=1000          # number of gate replicates per arm
  RUN_GRID=1              # run grid analysis (0 to skip)
  GRID_BUDGET_HOURS=5     # max time for grid (0 = unlimited)
  MC_PERSIST_DIR=/content/drive/MyDrive/...   # output directory
  MC_EXECUTOR=process     # "process" (default) or "thread"

WHAT'S FIXED/CHANGED VS. THE PREVIOUS DRAFT
--------------------------------------------
1. Real bug fix — RNG desync on resume (gate stage). `calibrate_alpha_gate()`
   was being called for BOTH arms on every resume, even for an arm that was
   already complete. That call draws from the shared RNG, so resuming after
   the "null" arm finished silently shifted every subsequent random draw for
   "signal" relative to an uninterrupted run. Alpha is now cached per arm in
   the checkpoint and only computed once, at the point it would occur in an
   uninterrupted run (same pattern the grid stage already used correctly).

2. Real bug fix — possible duplicate rows on resume. A crash between a
   batch's CSV append and its checkpoint save reprocesses and re-appends
   that batch on the next run, double-counting those replicates. Final
   aggregation now de-duplicates by (arm, rep) / (cell, rep) first.

3. Speed — default executor is a fork-based ProcessPoolExecutor instead of
   threads. Each fit makes many small Python-level calls into
   scipy.optimize.minimize, which stay largely GIL-bound under threads;
   processes sidestep the GIL and scale close to linearly with core count.
   Synthetic data is still generated serially in the parent before dispatch,
   so which process runs an already-fully-determined fit changes nothing
   about the output. Set MC_EXECUTOR=thread to opt back into threads.

4. Robustness: Drive writes (atomic_json/append_rows) retry with backoff on
   transient FUSE I/O errors instead of crashing the run. The final zip step
   is best-effort — a zip failure no longer looks like the whole run failed;
   your results are already safe on Drive throughout the run regardless.

5. Defensive fixes: MC_BATCH=0 no longer causes a silent infinite loop;
   GATE_REPS/GRID_REPS/GRID_BUDGET_HOURS/MC_WORKERS all raise a clear error
   on malformed input instead of an opaque traceback or a silently-ignored
   value. One redundant objective evaluation per Hessian removed (res.fun
   IS nll(res.x), so it's reused instead of recomputed). Unused imports
   (as_completed, contextmanager) and a leftover dead global removed.

Everything else — logging setup, tqdm progress bars, and the "require an
already-mounted Drive" check — is kept as you had it.
"""

import os
import sys
import json
import csv
import shutil
import time
import platform
import logging
import multiprocessing as mp
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor

# Limit BLAS threads to avoid thread storms in children. Must run before
# numpy/scipy are imported anywhere in the process (including via `base`).
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
          "NUMEXPR_NUM_THREADS", "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(v, "1")
os.environ.setdefault("MKL_DYNAMIC", "FALSE")

# -----------------------------------------------------------------------------
# Logging setup
# -----------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger(__name__)

# =============================================================================
# Env var helpers — clear errors instead of opaque tracebacks or silent
# fallbacks on a typo'd value.
# =============================================================================
def _env_int(name, default, min_value=None):
    raw = os.environ.get(name, str(default))
    try:
        val = int(raw)
    except ValueError:
        raise SystemExit(f"Environment variable {name}={raw!r} is not an integer.")
    if min_value is not None and val < min_value:
        logger.warning(f"{name}={val} is below the minimum of {min_value}; "
                        f"clamping to {min_value}.")
        val = min_value
    return val


def _env_float(name, default):
    raw = os.environ.get(name, str(default))
    try:
        return float(raw)
    except ValueError:
        raise SystemExit(f"Environment variable {name}={raw!r} is not a number.")


# =============================================================================
# Configuration
# =============================================================================
GATE_REPS = _env_int("GATE_REPS", 1000, min_value=0)
GRID_REPS = _env_int("GRID_REPS", 500, min_value=0)
RUN_GRID = os.environ.get("RUN_GRID", "1") != "0"
GRID_BUDGET_HOURS = _env_float("GRID_BUDGET_HOURS", 5.0)
BATCH_SIZE = _env_int("MC_BATCH", 8, min_value=1)
EXECUTOR_KIND = os.environ.get("MC_EXECUTOR", "process").strip().lower()

_workers_env = _env_int("MC_WORKERS", 0, min_value=0)
if _workers_env > 0:
    N_WORKERS = _workers_env
elif hasattr(os, "sched_getaffinity"):
    N_WORKERS = len(os.sched_getaffinity(0))
else:
    N_WORKERS = os.cpu_count() or 2
N_WORKERS = max(1, N_WORKERS)

SCRIPT_DIR = Path(__file__).resolve().parent
PERSIST_DIR = Path(os.environ.get(
    "MC_PERSIST_DIR",
    "/content/drive/MyDrive/mdwarf_radius_valley_highrep_fast_v3"
))

# -----------------------------------------------------------------------------
# Google Drive verification (no auto-mount — mount it yourself first; this
# never calls drive.mount(), so it can't re-trigger the mount/auth flow)
# -----------------------------------------------------------------------------
def verify_drive():
    """Check that /content/drive is mounted and writable."""
    drive_root = Path("/content/drive")
    if not drive_root.exists():
        raise SystemExit(
            "Google Drive not mounted at /content/drive.\n"
            "Please mount it manually (e.g., drive.mount('/content/drive'))."
        )
    if not (drive_root / "MyDrive").exists():
        raise SystemExit(
            "Google Drive mounted but 'MyDrive' folder not found – "
            "check mount point."
        )
    try:
        PERSIST_DIR.mkdir(parents=True, exist_ok=True)
        probe = PERSIST_DIR / ".drive_write_test"
        probe.write_text("ok")
        probe.unlink(missing_ok=True)
    except Exception as e:
        raise SystemExit(
            f"Google Drive is NOT writable at {PERSIST_DIR}\nReason: {e}"
        )
    logger.info(f"Verified persistent Drive directory: {PERSIST_DIR}")

verify_drive()

# -----------------------------------------------------------------------------
# Import original base module
# -----------------------------------------------------------------------------
sys.path.insert(0, str(SCRIPT_DIR))
try:
    import colab_highrep_gate as base
except ImportError as e:
    raise SystemExit(
        "Could not import colab_highrep_gate.py. "
        "Place both scripts in /content/.\n"
        f"Import error: {e}"
    )

# Apply execution settings to base
base.GATE_REPS = GATE_REPS
base.GRID_REPS = GRID_REPS
base.RUN_GRID = RUN_GRID
base.GRID_BUDGET_HOURS = GRID_BUDGET_HOURS
base.OUT_DIR = PERSIST_DIR
base.CKPT_DIR = PERSIST_DIR / "checkpoints"
base.OUT_DIR.mkdir(parents=True, exist_ok=True)
base.CKPT_DIR.mkdir(parents=True, exist_ok=True)

# Try to import tqdm for progress bars
try:
    from tqdm.auto import tqdm
    HAVE_TQDM = True
except ImportError:
    HAVE_TQDM = False
    tqdm = None

# =============================================================================
# Executor selection
# =============================================================================
def make_executor(n_workers: int):
    """Build the replicate-level executor.

    Defaults to a fork-based process pool: each fit is a pure function of
    already-generated synthetic data, so which OS process runs it cannot
    change any number in the output, and processes give real multi-core
    speedup where threads are GIL-limited for this kind of Python-heavy
    optimizer loop. Falls back to threads if a fork context can't be
    created, or if MC_EXECUTOR=thread is set explicitly.
    """
    if EXECUTOR_KIND == "thread":
        return ThreadPoolExecutor(max_workers=n_workers), "threaded_exact_rng"

    if EXECUTOR_KIND != "process":
        logger.warning(f"Unknown MC_EXECUTOR={EXECUTOR_KIND!r}; defaulting to 'process'.")

    try:
        ctx = mp.get_context("fork")
        return (
            ProcessPoolExecutor(max_workers=n_workers, mp_context=ctx),
            "process_pool_exact_rng",
        )
    except (ValueError, OSError) as exc:
        logger.warning(f"Could not start a fork-based process pool ({exc}); "
                        f"falling back to threads.")
        return ThreadPoolExecutor(max_workers=n_workers), "threaded_exact_rng"

# =============================================================================
# Fast fitting (no Hessian for gate replicates)
# =============================================================================
def fit_fast(mdl, fixed_su=True, fix_comp=True, need_hessian=False):
    """
    Same objective/bounds as base.fit(), but skips Hessian for gate replicates.
    """
    p = mdl["X"].shape[1]
    mu_se0, s_se0, mu_sn0, s_sn0 = mdl["init_comp"]

    if fixed_su:
        th0 = base.np.r_[
            base.np.zeros(p),
            [mu_se0, base.np.log(s_se0), mu_sn0, base.np.log(s_sn0)]
        ]
        if fix_comp:
            cb = [
                (mu_se0 - 0.06, mu_se0 + 0.06),
                (base.np.log(s_se0) - 0.12, base.np.log(s_se0) + 0.12),
                (mu_sn0 - 0.06, mu_sn0 + 0.06),
                (base.np.log(s_sn0) - 0.12, base.np.log(s_sn0) + 0.12),
            ]
        else:
            cb = [
                (0.7, 1.75), (-2.5, -0.2),
                (1.8, 4.2), (-2.3, -0.05),
            ]
        bounds = [(-25, 25)] * p + cb
    else:
        th0 = base.np.r_[
            base.np.zeros(p), base.np.array([-1.5]),
            [mu_se0, base.np.log(s_se0), mu_sn0, base.np.log(s_sn0)]
        ]
        bounds = (
            [(-25, 25)] * p + [(-6, 2)] +
            [
                (0.7, 1.75), (-2.5, -0.2),
                (1.8, 4.2), (-2.3, -0.05),
            ]
        )

    S = float(len(mdl["r_obs"]))
    nll = lambda t: base.neg_log_post(t, mdl, fixed_su=fixed_su, norm=S)

    res = base.minimize(
        nll, th0, method="L-BFGS-B", bounds=bounds,
        options={"maxiter": 2000}
    )

    if not need_hessian:
        return res.x, None, None, res

    # Central-difference Hessian (same as base.fit() / sandwich bread).
    # f0 reuses the optimizer's own final objective value instead of
    # re-evaluating nll(res.x) a second time (res.fun IS nll(res.x)).
    f0 = float(res.fun)
    n = len(res.x)
    eps = 1e-4
    H = base.np.zeros((n, n))
    for i in range(n):
        ei = base.np.zeros(n)
        ei[i] = eps
        H[i, i] = (nll(res.x + ei) - 2.0 * f0 + nll(res.x - ei)) / eps**2
        for j in range(i + 1, n):
            ej = base.np.zeros(n)
            ej[j] = eps
            H[i, j] = H[j, i] = (
                nll(res.x + ei + ej)
                - nll(res.x + ei - ej)
                - nll(res.x - ei + ej)
                + nll(res.x - ei - ej)
            ) / (4.0 * eps**2)
    H *= S   # as in original sandwich code
    return res.x, None, H, res

# -----------------------------------------------------------------------------
# Exact analytic sandwich meat (fixed_su=True) — left bit-for-bit identical
# to the original; this is the numerically sensitive part of the spec, so
# it isn't touched for speed.
# -----------------------------------------------------------------------------
def analytic_sandwich_se_fixed(mdl, theta, H_unnorm):
    """Compute cluster-robust SEs using analytic score + finite-difference bread."""
    X = mdl["X"]
    r = mdl["r_obs"]
    e = mdl["e_r"]
    hidx = mdl["hidx"]
    w_pl = mdl["w_pl"]

    nc = X.shape[1]
    beta = theta[:nc]
    mu_se = float(theta[nc])
    log_sse = float(theta[nc + 1])
    mu_sn = float(theta[nc + 2])
    log_ssn = float(theta[nc + 3])

    tau_se = base.np.exp(log_sse)
    tau_sn = base.np.exp(log_ssn)
    var_se = tau_se**2 + e**2
    var_sn = tau_sn**2 + e**2
    sig_se = base.np.sqrt(var_se)
    sig_sn = base.np.sqrt(var_sn)

    eta = base.np.clip(X @ beta, -30.0, 30.0)
    pi = base.expit(eta)

    lw_sn = base.norm_logpdf(r, mu_sn, sig_sn)
    lw_se = base.norm_logpdf(r, mu_se, sig_se)

    a_sn = lw_sn + base.np.log(pi)
    a_se = lw_se + base.np.log1p(-pi)
    den = base.logsumexp(base.np.column_stack([a_sn, a_se]), axis=1)
    q_sn = base.np.exp(a_sn - den)
    q_se = 1.0 - q_sn

    deta = q_sn - pi

    npar = len(theta)
    scores = base.np.zeros((len(r), npar))

    scores[:, :nc] = deta[:, None] * X
    scores[:, nc] = q_se * (r - mu_se) / var_se
    scores[:, nc + 1] = (q_se * (tau_se**2 / var_se) *
                         (((r - mu_se)**2 / var_se) - 1.0))
    scores[:, nc + 2] = q_sn * (r - mu_sn) / var_sn
    scores[:, nc + 3] = (q_sn * (tau_sn**2 / var_sn) *
                         (((r - mu_sn)**2 / var_sn) - 1.0))

    Sg = base.np.zeros((mdl["n_hosts"], npar))
    for j in range(npar):
        Sg[:, j] = base.np.bincount(hidx, weights=w_pl * scores[:, j],
                                    minlength=mdl["n_hosts"])

    meat = Sg.T @ Sg
    Binv = base.np.linalg.inv(H_unnorm + base.np.eye(npar) * 1e-6)
    V = Binv @ meat @ Binv
    return base.np.sqrt(base.np.clip(base.np.diag(V), 0.0, base.np.inf))

# =============================================================================
# Workers (module-level so they're picklable for ProcessPoolExecutor)
# =============================================================================
def gate_task(task):
    rep, mdl, iqr, sd = task
    theta, _, _, res = fit_fast(mdl, fixed_su=True, fix_comp=True,
                                need_hessian=False)
    if not base.np.all(base.np.isfinite(theta)):
        raise RuntimeError(f"Non-finite gate fit at replicate {rep}")
    b_iqr = float(theta[1] * (iqr / sd))
    return rep, b_iqr, len(mdl["r_obs"]), bool(res.success)

def grid_task(task):
    rep, mdl, beta_inj = task
    theta, _, H, res = fit_fast(mdl, fixed_su=True, fix_comp=True,
                                need_hessian=True)
    if not base.np.all(base.np.isfinite(theta)):
        raise RuntimeError(f"Non-finite grid fit at replicate {rep}")

    se = analytic_sandwich_se_fixed(mdl, theta, H)
    b = float(theta[1])
    s = float(se[1])

    excl0 = (b - 1.96 * s) * (b + 1.96 * s) > 0
    covers = (b - 1.96 * s <= beta_inj <= b + 1.96 * s)
    detect = bool(excl0 if beta_inj == 0
                 else (excl0 and base.np.sign(b) == base.np.sign(beta_inj)))

    return rep, b, s, int(detect), int(covers), len(mdl["r_obs"]), bool(res.success)

# =============================================================================
# Persistence helpers (with retry for transient Drive FUSE I/O errors)
# =============================================================================
def _drive_retry(fn, *args, attempts=5, base_delay=1.0, **kwargs):
    """Retry a Drive filesystem write on transient I/O errors, with backoff."""
    last_exc = None
    for attempt in range(1, attempts + 1):
        try:
            return fn(*args, **kwargs)
        except OSError as exc:
            last_exc = exc
            if attempt == attempts:
                break
            wait = base_delay * (2 ** (attempt - 1))
            logger.warning(f"Drive write failed ({exc}); "
                            f"retry {attempt}/{attempts} in {wait:.1f}s")
            time.sleep(wait)
    raise last_exc


def atomic_json(path, obj):
    path = Path(path)

    def _write():
        tmp = path.with_suffix(path.suffix + ".tmp")
        with open(tmp, "w") as fh:
            json.dump(obj, fh, indent=1)
            fh.flush()
            os.fsync(fh.fileno())
        tmp.replace(path)

    _drive_retry(_write)

def append_rows(path, header, rows):
    path = Path(path)

    def _write():
        new = not path.exists()
        with open(path, "a", newline="") as fh:
            w = csv.writer(fh)
            if new:
                w.writerow(header)
            w.writerows(rows)
            fh.flush()
            os.fsync(fh.fileno())

    _drive_retry(_write)

def load_ckpt(name):
    p = base.CKPT_DIR / f"{name}.json"
    if not p.exists():
        return None
    with open(p) as fh:
        return json.load(fh)

def save_ckpt(name, obj):
    atomic_json(base.CKPT_DIR / f"{name}.json", obj)

# =============================================================================
# Gate analysis
# =============================================================================
def run_gate_fast(df, ess):
    gate_csv = base.OUT_DIR / (
        f"task2_gate_rep{base.GATE_REPS}_replicates_raw.csv"
    )
    summary_path = base.OUT_DIR / (
        f"task2_calibration_kepler_only_trim0.95_rep{base.GATE_REPS}.json"
    )

    header = ["arm", "rep", "beta_inj", "alpha",
              "beta_hat_per_iqr", "n_planets"]

    logger.info("=== GATE REAL FIT ===")
    b_real, se_real = base.fit_real_gate(df)
    sd = df.vtan.std()
    iqr = df.vtan.quantile(0.75) - df.vtan.quantile(0.25)
    logger.info(f"REAL beta_age/IQR = {b_real:+.4f} +/- {se_real:.4f}")

    ck = load_ckpt("gate_fast_v3")
    if ck and ck.get("reps") == base.GATE_REPS:
        rng = base.restore_rng(ck["rng_state"])
        done = ck.get("done", {"null": 0, "signal": 0})
        alphas = ck.get("alpha", {})
        logger.info(f"Resuming gate: null={done.get('null',0)}/{base.GATE_REPS}, "
                    f"signal={done.get('signal',0)}/{base.GATE_REPS}")
    elif ck:
        raise SystemExit("Existing gate checkpoint uses a different GATE_REPS.")
    else:
        rng = base.np.random.default_rng(base.SEED_GATE)
        done = {"null": 0, "signal": 0}
        alphas = {}

    pool, mode_tag = make_executor(N_WORKERS)
    logger.info(f"[gate] executor: {mode_tag} ({N_WORKERS} workers)")

    with pool:
        for arm, beta_inj in [("null", 0.0), ("signal", base.B_INJ)]:
            # Alpha is calibrated exactly once per arm, at the same point in
            # the RNG stream an uninterrupted run would hit it. Recomputing
            # it on every resume (even for an already-complete arm) was the
            # RNG-desync bug — reuse the cached value instead.
            if arm in alphas:
                alpha = alphas[arm]
            else:
                alpha = base.calibrate_alpha_gate(df, base.GATE_TARGET_FRAC, rng)
                alphas[arm] = alpha
                save_ckpt("gate_fast_v3",
                          {"reps": base.GATE_REPS, "done": done,
                           "alpha": alphas,
                           "rng_state": rng.bit_generator.state,
                           "mode": mode_tag})

            start = int(done.get(arm, 0))
            if start >= base.GATE_REPS:
                logger.info(f"[{arm}] already complete")
                continue

            pbar = tqdm(total=base.GATE_REPS, initial=start,
                        desc=f"Gate {arm}", disable=not HAVE_TQDM)
            while start < base.GATE_REPS:
                n = min(BATCH_SIZE, base.GATE_REPS - start)
                tasks = []
                for rep in range(start, start + n):
                    ds = base.synth_gate(df, beta_inj, alpha, rng)
                    mdl = base.build_design(ds, "vtan", base.GATE_VALLEY,
                                            fe=False, weights=ds.wgt.values)
                    tasks.append((rep, mdl, iqr, sd))

                t0 = time.time()
                results = list(pool.map(gate_task, tasks))
                dt = time.time() - t0
                results.sort(key=lambda x: x[0])

                append_rows(gate_csv, header,
                            [[arm, rep, beta_inj, alpha, b, nplan]
                             for rep, b, nplan, _ in results])

                start += n
                done[arm] = start
                save_ckpt("gate_fast_v3",
                          {"reps": base.GATE_REPS, "done": done,
                           "alpha": alphas,
                           "next_rep": start, "arm": arm,
                           "rng_state": rng.bit_generator.state,
                           "mode": mode_tag})

                successes = sum(r[3] for r in results)
                pbar.update(n)
                logger.info(f"[gate/{arm}] {start}/{base.GATE_REPS} | "
                            f"batch={n} | {dt:.2f}s | "
                            f"{dt/max(n,1):.2f}s/rep | "
                            f"optimizer_success={successes}/{n}")
            pbar.close()

    raw = base.pd.read_csv(gate_csv, keep_default_na=False)

    # Defensive de-duplication: a crash between the CSV append and the
    # checkpoint save for the same batch would otherwise reprocess and
    # re-append the same replicates on the next resume. Rows for the same
    # (arm, rep) are numerically identical regardless of which write "won".
    dup_mask = raw.duplicated(subset=["arm", "rep"], keep="last")
    if dup_mask.any():
        logger.info(f"[gate] removed {int(dup_mask.sum())} duplicate "
                     f"replicate row(s) left over from an interrupted "
                     f"batch write — safe to ignore.")
        raw = raw[~dup_mask].reset_index(drop=True)

    res = {"beta_real": float(b_real), "se_real": float(se_real),
           "reps": base.GATE_REPS, "n_planets": int(len(df)),
           "ess": float(ess), "seed": base.SEED_GATE,
           "generated_by": "colab_highrep_gate_fast_v3.py",
           "execution": mode_tag,
           "workers": N_WORKERS, "batch_size": BATCH_SIZE}

    res["reps_completed"] = {a: int((raw.arm == a).sum())
                             for a in ("null", "signal")}

    for arm in ("null", "signal"):
        vals = raw[raw.arm == arm].sort_values("rep").beta_hat_per_iqr.values
        res[arm] = {"median": float(base.np.median(vals)),
                    "q025": float(base.np.quantile(vals, .025)),
                    "q975": float(base.np.quantile(vals, .975)),
                    "pct_rank_of_real": float(base.np.mean(vals <= b_real))}

    p_null_upper = 1.0 - res["null"]["pct_rank_of_real"]
    inside_signal = (res["signal"]["q025"] <= b_real <= res["signal"]["q975"])
    res["verdict"] = {
        "distinguishable_from_null": bool(p_null_upper <= 0.05),
        "consistent_with_signal": bool(inside_signal),
    }
    res["verdict"]["PASS"] = bool(res["verdict"]["distinguishable_from_null"] and
                                  res["verdict"]["consistent_with_signal"])

    atomic_json(summary_path, res)

    logger.info(f"Gate complete: P(null>=real)={p_null_upper:.4f}; "
                f"inside_signal={inside_signal}; "
                f"VERDICT={'PASS' if res['verdict']['PASS'] else 'FAIL'}")
    return time.time()

# =============================================================================
# Grid analysis (with time budget)
# =============================================================================
def run_grid_fast(base_m, t_grid_start):
    grid_csv = base.OUT_DIR / f"power_grid_rep{base.GRID_REPS}_replicates_raw.csv"
    grid_json = base.OUT_DIR / f"power_analysis_rep{base.GRID_REPS}.json"

    header = ["cell", "scale", "effect", "beta_inj", "rep",
              "beta_hat_per_sd", "se_per_sd", "detect", "cover", "n_planets"]

    real_frac = float((base_m.pl_rade > base.GRID_VALLEY).mean())
    calib_frac = min(0.495, real_frac / 0.90)
    f_kep = float((base_m.fam == "Kepler").mean())
    f_k2 = float((base_m.fam == "K2").mean())
    mission_offset = base.GAMMA_KEP * f_kep + base.GAMMA_K2 * f_k2

    cells = []
    eff_seed = {"null": 0, "SWEET": 1, "strong": 2}
    for scale, label in [(0.5, "0.5x"), (1.0, "1x"), (2.0, "2x"),
                         (4.0, "4x"), (8.0, "8x")]:
        for effect, beta_inj in [("null", 0.0), ("SWEET", -0.14),
                                 ("strong", -0.28)]:
            seed = eff_seed[effect] * 1000 + int(round(scale * 100))
            cells.append((scale, label, effect, beta_inj, seed))

    deadline = (t_grid_start + GRID_BUDGET_HOURS * 3600.0
                if GRID_BUDGET_HOURS > 0 else float("inf"))

    def make_grid_tasks(scale, effect, beta_inj, alpha, rng, start, n):
        tasks = []
        for rep in range(start, start + n):
            d = base.synth_grid(base_m, scale, beta_inj, alpha, rng)
            mdl = base.build_design(d, "vtan", base.GRID_VALLEY, fe=False)
            tasks.append((rep, mdl, beta_inj))
        return tasks

    skipped = []

    pool, mode_tag = make_executor(N_WORKERS)
    logger.info(f"[grid] executor: {mode_tag} ({N_WORKERS} workers)")

    with pool:
        for scale, label, effect, beta_inj, seed in cells:
            ckey = f"grid_{label}_{effect}_fast_v3"
            cval = f"{label}|{effect}"

            if time.time() > deadline:
                skipped.append(cval)
                logger.info(f"[{cval}] SKIPPED: grid budget exhausted")
                continue

            ck = load_ckpt(ckey)
            if ck and ck.get("reps") == base.GRID_REPS:
                start = int(ck.get("next_rep", 0))
                alpha = ck["alpha"]
                rng = base.restore_rng(ck["rng_state"])
            elif ck:
                raise SystemExit(f"Checkpoint {ckey} uses a different GRID_REPS.")
            else:
                rng = base.np.random.default_rng(seed)
                alpha = base.calibrate_alpha_grid(
                    base_m, beta_inj, calib_frac, rng,
                    mission_offset=mission_offset
                )
                start = 0

            grid_csv.parent.mkdir(parents=True, exist_ok=True)
            logger.info(f"[{cval}] start={start}/{base.GRID_REPS}; "
                        f"alpha={alpha:+.4f}")

            pbar = tqdm(total=base.GRID_REPS, initial=start,
                        desc=f"Grid {cval}", disable=not HAVE_TQDM)
            while start < base.GRID_REPS:
                if time.time() > deadline:
                    skipped.append(cval)
                    logger.info(f"[{cval}] STOPPED: grid budget exhausted")
                    break

                n = min(BATCH_SIZE, base.GRID_REPS - start)
                tasks = make_grid_tasks(scale, effect, beta_inj,
                                        alpha, rng, start, n)

                t0 = time.time()
                results = list(pool.map(grid_task, tasks))
                dt = time.time() - t0
                results.sort(key=lambda x: x[0])

                append_rows(grid_csv, header,
                            [[cval, label, effect, beta_inj, rep, b, s,
                              detect, cover, nplan]
                             for rep, b, s, detect, cover, nplan, _ in results])

                start += n
                successes = sum(r[6] for r in results)
                pbar.update(n)

                save_ckpt(ckey,
                          {"reps": base.GRID_REPS, "next_rep": start,
                           "alpha": alpha,
                           "rng_state": rng.bit_generator.state,
                           "seed": seed, "mode": mode_tag})

                logger.info(f"[{cval}] {start}/{base.GRID_REPS} | "
                            f"batch={n} | {dt:.2f}s | "
                            f"{dt/max(n,1):.2f}s/rep | "
                            f"optimizer_success={successes}/{n}")
            pbar.close()

            if start >= base.GRID_REPS:
                logger.info(f"[{cval}] COMPLETE")

    if not grid_csv.exists():
        return skipped

    raw = base.pd.read_csv(grid_csv, keep_default_na=False)

    # Same defensive de-duplication as the gate stage (see comment there).
    dup_mask = raw.duplicated(subset=["cell", "rep"], keep="last")
    if dup_mask.any():
        logger.info(f"[grid] removed {int(dup_mask.sum())} duplicate "
                     f"replicate row(s) left over from an interrupted "
                     f"batch write — safe to ignore.")
        raw = raw[~dup_mask].reset_index(drop=True)

    results = []
    for _scale, label, effect, beta_inj, _seed in cells:
        sub = raw[raw.cell == f"{label}|{effect}"].sort_values("rep")
        if len(sub) == 0:
            continue
        results.append({
            "scale": label, "effect": effect,
            "beta_inj": float(beta_inj),
            "reps": int(len(sub)),
            "detect_rate": float(sub.detect.mean()),
            "coverage": float(sub.cover.mean()),
            "median_beta_hat": float(base.np.median(sub.beta_hat_per_sd.values))
        })

    atomic_json(grid_json, results)
    return skipped

# =============================================================================
# Main
# =============================================================================
def main():
    logger.info("=" * 78)
    logger.info("HIGH-REPLICATE MC — FAST V3 / EXACT RNG / DRIVE-PERSISTENT")
    logger.info("=" * 78)
    logger.info(f"Platform: {platform.platform()}")
    logger.info(f"Python: {platform.python_version()}")
    logger.info(f"CPU workers: {N_WORKERS}")
    logger.info(f"Executor: {EXECUTOR_KIND}")
    logger.info(f"Batch size: {BATCH_SIZE}")
    logger.info(f"Gate: {GATE_REPS}/arm")
    logger.info(f"Grid: {GRID_REPS}/cell")
    logger.info(f"Grid budget: {GRID_BUDGET_HOURS} h")
    logger.info(f"Drive: {base.OUT_DIR}")
    logger.info("STATISTICAL SPECIFICATION: UNCHANGED")
    logger.info("RNG STREAM: exact serial synthetic-generation order")
    logger.info("Gate synthetic fits: unused covariance calculation removed")
    logger.info("Grid sandwich: exact analytic score + same sandwich algebra")
    logger.info("=" * 78)

    df_gate, ess = base.load_weighted_fgk()
    base.validate_gate_df(df_gate, ess)

    base_m = base.load_base_m()
    if not (400 <= len(base_m) <= 430):
        raise SystemExit(f"ABORT: M sample N={len(base_m)}, expected ~417.")
    logger.info(f"M sample validated: N={len(base_m)}, "
                f"hosts={base_m.hostname.nunique()}")

    atomic_json(base.OUT_DIR / "accelerator_manifest.json",
                {"script": "colab_highrep_gate_fast_v3.py",
                 "gate_reps": GATE_REPS, "grid_reps": GRID_REPS,
                 "workers": N_WORKERS, "executor": EXECUTOR_KIND,
                 "batch_size": BATCH_SIZE,
                 "grid_budget_hours": GRID_BUDGET_HOURS,
                 "drive_persistence_verified": True,
                 "rng_stream": "original serial order preserved",
                 "model_spec_unchanged": True,
                 "note": "Synthetic generation serial; fitting parallelized "
                         "(process pool by default)."})

    t0 = time.time()
    t_grid_start = run_gate_fast(df_gate, ess)

    skipped = []
    if RUN_GRID:
        skipped = run_grid_fast(base_m, t_grid_start)

    manifest = {"generated_by": "colab_highrep_gate_fast_v3.py",
                "gate_reps_per_arm": GATE_REPS,
                "grid_reps_per_cell": GRID_REPS,
                "grid_budget_hours": GRID_BUDGET_HOURS,
                "workers": N_WORKERS, "executor": EXECUTOR_KIND,
                "batch_size": BATCH_SIZE,
                "grid_skipped_cells": skipped,
                "elapsed_hours": (time.time() - t0) / 3600.0,
                "drive_output": str(base.OUT_DIR),
                "model_spec_unchanged": True, "rng_stream_exact": True}
    atomic_json(base.OUT_DIR / "manifest.json", manifest)

    try:
        archive = shutil.make_archive(
            str(base.OUT_DIR.parent / "colab_gate_output_fast_v3"),
            "zip", base.OUT_DIR
        )
    except Exception as exc:
        archive = None
        logger.warning(f"Could not create the zip archive ({exc}); this "
                        f"does not affect your results, which are already "
                        f"saved at {base.OUT_DIR}.")

    logger.info("=" * 78)
    logger.info("DONE")
    logger.info(f"Persistent directory: {base.OUT_DIR}")
    if archive:
        logger.info(f"Archive: {archive}")
    logger.info("=" * 78)

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        logger.warning("Interrupted. Progress is checkpointed on Drive — "
                        "rerun this script to resume from where it left off.")
        sys.exit(130)
