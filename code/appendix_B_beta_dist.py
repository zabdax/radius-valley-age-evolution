#!/usr/bin/env python3
"""appendix_B_beta_dist.py — beta-distribution injection sensitivity.

Literature effect-size uncertainty: Kamulali+2026 SE/SN ratios
  young 0.51 (+0.11/-0.08), old 0.64 (+0.11/-0.11)
are approximated as independent Normals sd=0.10 each (documented
approximation; asymmetric bounds averaged). Per replicate:
  r1 ~ N(0.51, 0.10), r2 ~ N(0.64, 0.10)  [redrawn if <=0.05]
  f1 = 1/(1+r1), f2 = 1/(1+r2]
  beta_inj = -(logit(f1) - logit(f2)) / 1.6   (median-split conversion)
which yields beta_inj ~ median -0.14, sd ~0.16 (verified in analysis:
median -0.142, sd 0.165, q025 -0.31, q975 +0.04 approx; exact quantiles
recomputed at runtime and saved).

Runs M pooled estimator, 300 replicates, per-replicate beta_inj draws,
legacy host-effect mode (archived behavior) + independent mode (150 each
would double cost; here: 300 legacy + 150 independent).
Compares empirical detection rate under distributed-beta vs point -0.14.

Outputs results/appendix/: beta_dist_M_raw.csv, beta_dist_summary.json.
"""
import json
import multiprocessing as mp
import os
import pathlib
import sys
import time

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
import importlib.util as _ilu
_spec = _ilu.spec_from_file_location("s15", HERE / "15_stratified_gate.py")
s15 = _ilu.module_from_spec(_spec)
_spec.loader.exec_module(s15)

OUT = ROOT / "results" / "appendix"
OUT.mkdir(parents=True, exist_ok=True)

REPS_LEGACY = 300
REPS_INDEP = 150
SEED = 9001
R_SD = 0.10

_W = {}


def draw_beta_inj(rng, n):
    r1 = rng.normal(0.51, R_SD, n)
    r2 = rng.normal(0.64, R_SD, n)
    bad = (r1 <= 0.05) | (r2 <= 0.05) | (~np.isfinite(r1)) | (~np.isfinite(r2))
    while bad.any():
        r1[bad] = rng.normal(0.51, R_SD, bad.sum())
        r2[bad] = rng.normal(0.64, R_SD, bad.sum())
        bad = ((r1 <= 0.05) | (r2 <= 0.05) | (~np.isfinite(r1))
               | (~np.isfinite(r2)))
    f1 = 1 / (1 + r1)
    f2 = 1 / (1 + r2)
    logit = lambda f: np.log(f / (1 - f))
    return -(logit(f1) - logit(f2)) / 1.6


def _init(base, comp, alpha, betas, independent_hosts):
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "NUMEXPR_NUM_THREADS"):
        os.environ[v] = "1"
    _W.update(base=base, comp=comp, alpha=alpha, betas=np.asarray(betas),
              independent_hosts=independent_hosts)


def _one(task):
    rep, seedseq = task
    rng = np.random.default_rng(seedseq)
    bi = float(_W["betas"][rep])
    d = s15.synth(_W["base"], bi, _W["alpha"], rng,
                  independent_hosts=_W["independent_hosts"])
    r = s15.pooled_estimator(d, _W["comp"])
    if r is None:
        return {"rep": rep, "beta_inj": bi, "beta_hat": float("nan"),
                "se": float("nan"), "n": len(d)}
    b, s = r
    excl0 = (b - 1.96 * s) * (b + 1.96 * s) > 0
    det = bool(excl0 and (np.sign(b) == np.sign(bi)) and bi != 0)
    return {"rep": rep, "beta_inj": bi, "beta_hat": b, "se": s,
            "detect_correct_sign": det, "n": len(d)}


def run_arm(sample, independent_hosts, reps, seed):
    cfg = s15.SAMPLES[sample]
    s15.CFG = cfg
    s15.TARGET_FRAC, s15.VALLEY = cfg["target_frac"], cfg["valley"]
    s15.GAMMA_LOGP = cfg["gamma_logp"]
    base = s15.load_control(cfg)
    full = s15.m9.build_design(base, "vtan", valley=cfg["valley"], fe=False,
                               weights=base.wgt.values)
    comp = full["init_comp"]
    alpha = s15.calibrate_alpha(base, cfg["target_frac"])
    rng0 = np.random.default_rng(seed)
    betas = draw_beta_inj(rng0, reps)
    children = np.random.SeedSequence(seed + 1).spawn(reps)
    todo = [(i, children[i]) for i in range(reps)]
    workers = max(1, (os.cpu_count() or 1) - 1)
    t0 = time.time()
    init = (base, comp, alpha, betas, independent_hosts)
    rows = []
    tag = f"betadist_{sample}_{'indep' if independent_hosts else 'legacy'}"
    print(f"[{tag}] reps={reps} beta_inj median {np.median(betas):+.3f} "
          f"sd {np.std(betas, ddof=1):.3f} "
          f"[{np.quantile(betas, .025):+.3f},{np.quantile(betas, .975):+.3f}]",
          flush=True)
    if workers <= 1:
        _init(*init)
        for k, task in enumerate(todo, 1):
            rows.append(_one(task))
            if k % 25 == 0 or k == reps:
                print(f"  [{tag}] {k}/{reps}", flush=True)
    else:
        chunk = max(1, min(8, len(todo) // (workers * 4) or 1))
        method = ("fork" if "fork" in mp.get_all_start_methods() else "spawn")
        with mp.get_context(method).Pool(
                workers, initializer=_init, initargs=init) as pool:
            for k, row in enumerate(
                    pool.imap_unordered(_one, todo, chunksize=chunk), 1):
                rows.append(row)
                if k % 25 == 0 or k == reps:
                    print(f"  [{tag}] {k}/{reps} "
                          f"{(time.time()-t0)/k:.2f}s/rep", flush=True)
    df = pd.DataFrame(rows).sort_values("rep")
    df.insert(0, "mode", "independent" if independent_hosts else "legacy")
    df.insert(0, "sample", sample)
    return df, betas


def main():
    dfs = []
    for indep, reps, seed in ((False, REPS_LEGACY, SEED),
                              (True, REPS_INDEP, SEED + 100)):
        df, betas = run_arm("M", indep, reps, seed)
        dfs.append(df)
        df.to_csv(OUT / f"beta_dist_M_{'indep' if indep else 'legacy'}_raw.csv",
                  index=False)
    full = pd.concat(dfs, ignore_index=True)
    full.to_csv(OUT / "beta_dist_M_all_raw.csv", index=False)
    summ = {}
    for mode in ("legacy", "independent"):
        d = full[full["mode"] == mode]
        bh = d.beta_hat.values
        bi = d.beta_inj.values
        ok = np.isfinite(bh)
        bh, bi = bh[ok], bi[ok]
        det = d.detect_correct_sign.values[ok]
        resid = bh - bi
        summ[mode] = {
            "n": int(ok.sum()),
            "beta_inj": {"median": float(np.median(bi)),
                         "sd": float(np.std(bi, ddof=1)),
                         "q025": float(np.quantile(bi, 0.025)),
                         "q975": float(np.quantile(bi, 0.975)),
                         "frac_near_zero_lt0.05": float(np.mean(abs(bi) < 0.05)),
                         "frac_positive": float(np.mean(bi > 0))},
            "beta_hat": {"median": float(np.median(bh)),
                         "sd": float(np.std(bh, ddof=1)),
                         "q025": float(np.quantile(bh, 0.025)),
                         "q975": float(np.quantile(bh, 0.975))},
            "detect_correct_sign_rate": float(np.mean(det)),
            "resid_hat_minus_inj": {"median": float(np.median(resid)),
                                    "sd": float(np.std(resid, ddof=1))},
            "corr_inj_hat": float(np.corrcoef(bi, bh)[0, 1]) if len(bh) > 2
            else float("nan")}
    summ["meta"] = {"R_SD": R_SD, "conversion": "median-split 1.6 SD",
                    "assumption": "independent Normal(0.10) on each ratio; "
                                  "redraw <=0.05",
                    "point_signal_reference": "M 1x point-beta detection 9.4% "
                                              "(500/cell production grid)"}
    with open(OUT / "beta_dist_summary.json", "w") as fh:
        json.dump(summ, fh, indent=1)
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
