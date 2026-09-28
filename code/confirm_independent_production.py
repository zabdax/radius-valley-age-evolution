#!/usr/bin/env python3
"""confirm_independent_production.py — production-scale independent-host confirmation.

Pooled-only arms with independent_hosts=True:
  M:   1000 reps/arm (null beta=0, signal beta=-0.14)
  FGK: 400 reps/arm
Compares against legacy production values (M arm -0.19, FGK arm -1.01).
Writes results/appendix/confirm_independent_{M,FGK}_raw.csv + confirm_summary.json.
"""
import json
import pathlib
import os
import sys
import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import appendix_A_arm_delta as appA

OUT = appA.OUT
REPS = {"M": 1000, "FGK": 400}
SEEDS = {("M", "null"): 7101, ("M", "signal"): 7102,
         ("FGK", "null"): 7111, ("FGK", "signal"): 7112}
REAL = {"M": 0.0458, "FGK": -0.0296}


def main():
    workers = max(1, (os.cpu_count() or 1) - 1)
    print(f"workers={workers}", flush=True)
    frames = []
    for sample in ("M", "FGK"):
        for arm, beta in (("null", 0.0), ("signal", -0.14)):
            df = appA.run_one(sample, "independent", arm, beta, REPS[sample],
                              SEEDS[(sample, arm)], workers)
            df.to_csv(OUT / f"confirm_independent_{sample}_{arm}_raw.csv",
                      index=False)
            frames.append(df)
    full = pd.concat(frames, ignore_index=True)
    full.to_csv(OUT / "confirm_independent_all_raw.csv", index=False)
    summ = {}
    for sample in ("M", "FGK"):
        sub = full[full["sample"] == sample]
        arms = {}
        for arm in ("null", "signal"):
            b = sub[sub["arm"] == arm]["beta_hat"].values
            b = b[np.isfinite(b)]
            arms[arm] = {"n": int(len(b)), "median": float(np.median(b)),
                         "sd": float(np.std(b, ddof=1)),
                         "q025": float(np.quantile(b, 0.025)),
                         "q975": float(np.quantile(b, 0.975)),
                         "pct_rank_of_real": float(np.mean(b <= REAL[sample]))}
        sep = ((arms["signal"]["median"] - arms["null"]["median"])
               / arms["null"]["sd"])
        summ[sample] = {"arms": arms, "arm_separation_nullSD": float(sep)}
    summ["meta"] = {"reps": REPS, "mode": "independent",
                    "production_legacy": {"M_arm": -0.19, "FGK_arm": -1.01}}
    with open(OUT / "confirm_summary.json", "w") as fh:
        json.dump(summ, fh, indent=1)
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
