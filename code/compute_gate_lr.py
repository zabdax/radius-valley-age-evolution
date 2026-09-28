#!/usr/bin/env python3
"""compute_gate_lr.py — reproduce Table gate likelihood ratios (proof of method).

Method (now documented in Sec gate): Gaussian KDE (scipy gaussian_kde,
Scott's bandwidth) fitted separately to each arm's beta_hat replicates,
ratio signal/null evaluated at the real coefficient.

Sources:
  FGK pooled     : results/highrep/task2_gate_rep1000_replicates_raw.csv
                   (real -0.0296) -> expect ~0.68
  FGK stratified : results/stratified_gate/strat_gate_mag_{null,signal}_r400_raw.csv
                   beta_strat column (real -0.1482) -> expect ~1.74
  M pooled       : results/appendix/arm_delta_M_legacy_{null,signal}_raw.csv
                   200/arm reduced ensemble (real +0.0458) -> expect ~1.0
                   (production 0.93 from unarchived 1000/arm run)
Read replicate CSVs with keep_default_na=False (arm label 'null').
"""
import pathlib
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

ROOT = pathlib.Path(__file__).resolve().parents[1]
R = ROOT / "results"


def lr(null_vals, sig_vals, real):
    kn = gaussian_kde(np.asarray(null_vals, float))
    ks = gaussian_kde(np.asarray(sig_vals, float))
    return float(ks(real)[0] / kn(real)[0]), float(kn.factor)


out = {}
g = pd.read_csv(R / "highrep" / "task2_gate_rep1000_replicates_raw.csv",
                keep_default_na=False)
n = g[g["arm"] == "null"]["beta_hat_per_iqr"].values
s = g[g["arm"] == "signal"]["beta_hat_per_iqr"].values
v, bw = lr(n, s, -0.0296)
out["FGK_pooled_1000arm"] = {"LR": v, "bandwidth_factor": bw, "expect": 0.68}
print(f"FGK pooled (1000/arm archived): LR={v:.3f} (ms Table 0.68)")

n2 = pd.read_csv(R / "stratified_gate" / "strat_gate_mag_null_r400_raw.csv")
s2 = pd.read_csv(R / "stratified_gate" / "strat_gate_mag_signal_r400_raw.csv")
v2, bw2 = lr(n2["beta_strat"].values, s2["beta_strat"].values, -0.1482)
out["FGK_stratified_400arm"] = {"LR": v2, "bandwidth_factor": bw2,
                                "expect": 1.74}
print(f"FGK stratified (400/arm archived): LR={v2:.3f} (ms Table 1.74)")

a = pd.read_csv(R / "appendix" / "arm_delta_all_raw.csv",
                keep_default_na=False)
m0 = a[(a["sample"] == "M") & (a["mode"] == "legacy")
       & (a["arm"] == "null")]["beta_hat"].values
m1 = a[(a["sample"] == "M") & (a["mode"] == "legacy")
       & (a["arm"] == "signal")]["beta_hat"].values
v3, bw3 = lr(m0, m1, 0.0458)
out["M_pooled_200arm_appendix"] = {"LR": v3, "bandwidth_factor": bw3,
                                   "note": "reduced ensemble; prod 0.93"}
print(f"M pooled (200/arm appendix legacy): LR={v3:.3f} (production 0.93)")

import json
with open(R / "appendix" / "gate_lr_check.json", "w") as fh:
    json.dump(out, fh, indent=1)
print("wrote results/appendix/gate_lr_check.json")
