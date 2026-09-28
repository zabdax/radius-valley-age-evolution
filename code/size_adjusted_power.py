#!/usr/bin/env python3
"""size_adjusted_power.py — empirical-quantile detection rates for the rep500 grid.

Nominal `detect` uses the analytical sandwich 95% CI (known anti-conservative:
pooled FPR 6.68%, coverage 91.3%). This script recomputes detection per cell
against the empirical null distribution at the same scale (one-sided,
correct sign: beta_hat <= q025 of the null cell for negative injections),
plus Wilson 95% CIs for rates and coverage.

Reads results/highrep/power_grid_rep500_replicates_raw.csv
  (columns: scale, effect, beta_inj, beta_hat_per_sd, se_per_sd, detect, cover)
Writes results/appendix/size_adjusted_power.csv + prints headline cells.
"""
import pathlib
import numpy as np
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
R = ROOT / "results"


def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    m = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (max(0.0, (c - m) / d), min(1.0, (c + m) / d))


df = pd.read_csv(R / "highrep" / "power_grid_rep500_replicates_raw.csv",
                 keep_default_na=False)
rows = []
for scale in ["0.5x", "1x", "2x", "4x", "8x"]:
    null = df[(df["scale"] == scale)
              & (df["effect"].isna() | (df["effect"] == "null"))]
    if len(null) == 0:  # effect column may encode null differently
        null = df[(df["scale"] == scale) & (df["beta_inj"] == 0.0)]
    q025 = float(np.quantile(null["beta_hat_per_sd"].values, 0.025))
    for eff in ["SWEET", "strong"]:
        cell = df[(df["scale"] == scale) & (df["effect"] == eff)]
        if len(cell) == 0:
            continue
        b = cell["beta_hat_per_sd"].values
        adj = b <= q025  # negative injections: correct-sign tail rule
        k = int(adj.sum())
        n = len(b)
        lo, hi = wilson(k, n)
        nom = float(cell["detect"].mean())
        cov = float(cell["cover"].mean())
        clo, chi = wilson(int(cell["cover"].sum()), n)
        rows.append({"scale": scale, "effect": eff, "n": n,
                     "null_q025": round(q025, 4),
                     "adj_detect": round(k / n, 4), "k": k,
                     "adj_wilson_lo": round(lo, 4), "adj_wilson_hi": round(hi, 4),
                     "nominal_detect": round(nom, 4),
                     "coverage": round(cov, 4),
                     "cov_wilson_lo": round(clo, 4),
                     "cov_wilson_hi": round(chi, 4)})
out = pd.DataFrame(rows)
out.to_csv(R / "appendix" / "size_adjusted_power.csv", index=False)
print(out.to_string(index=False))
print("wrote results/appendix/size_adjusted_power.csv")
