import json
import numpy as np
import pandas as pd

import pathlib
APP = pathlib.Path(__file__).resolve().parents[1] / "results" / "appendix"
df = pd.read_csv(APP / "arm_delta_all_raw.csv", keep_default_na=False)
REAL = {"M": 0.0458, "FGK": -0.0296}
summary = {}
for sample in ("M", "FGK"):
    sub = df[df["sample"] == sample]
    out = {}
    for mode in ("legacy", "independent"):
        arms = {}
        for arm in ("null", "signal"):
            sel = sub[(sub["mode"] == mode) & (sub["arm"] == arm)]
            b = sel["beta_hat"].values
            se = sel["se"].values
            b = b[np.isfinite(b)]
            se = se[np.isfinite(se)]
            n = min(len(b), len(se))
            b, se = b[:n], se[:n]
            det = float(np.mean((b - 1.96 * se) * (b + 1.96 * se) > 0))
            arms[arm] = {"n": int(len(b)),
                         "median": float(np.median(b)),
                         "sd": float(np.std(b, ddof=1)) if len(b) > 1 else 0.0,
                         "q025": float(np.quantile(b, 0.025)),
                         "q975": float(np.quantile(b, 0.975)),
                         "pct_rank_of_real": float(np.mean(b <= REAL[sample])),
                         "analytical_det_rate": det}
        armsep = ((arms["signal"]["median"] - arms["null"]["median"])
                  / arms["null"]["sd"]) if arms["null"]["sd"] > 0 else float("nan")
        out[mode] = {"arms": arms, "arm_separation_nullSD": float(armsep)}
    deltas = {}
    for arm in ("null", "signal"):
        for stat in ("median", "sd", "pct_rank_of_real", "analytical_det_rate"):
            deltas[f"{arm}_{stat}"] = float(
                out["independent"]["arms"][arm][stat]
                - out["legacy"]["arms"][arm][stat])
    deltas["arm_separation"] = float(out["independent"]["arm_separation_nullSD"]
                                     - out["legacy"]["arm_separation_nullSD"])
    out["delta_independent_minus_legacy"] = deltas
    summary[sample] = out
summary["meta"] = {"reps": {"M": 200, "FGK": 100}, "B_SIGNAL": -0.14,
                   "REAL": REAL,
                   "note": "pooled estimator only; components pinned; TRIM=0.95; "
                           "read raw CSVs with keep_default_na=False ('null' string)"}
with open(APP / "arm_delta_summary.json", "w") as fh:
    json.dump(summary, fh, indent=1)
print(json.dumps(summary, indent=1))
