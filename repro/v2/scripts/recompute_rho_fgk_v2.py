#!/usr/bin/env python3
"""recompute_rho_fgk_v2.py -- independent computational verification of the
manuscript's age-proxy Spearman correlations.

Reads ONLY the frozen CSVs under radius-valley-age-evolution/data/ and writes
ONLY to manuscripts/q1-paper/v2-2026-11/results/audit/.

For each sample (FGK control, M/cool stars):
  - host_ages_<TAG>.csv  : hostname, age_gyr, censored, gyro_ok, ...
  - hosts_kinematics_<TAG>.csv : hostname, vtan, ...
  Selection: gyro_ok == True, finite age_gyr AND finite v_tan (inner join on
  hostname, which is unique in both kinematics files).
  Statistic: scipy.stats.spearmanr(vtan, age_gyr) -> exact p-value (t approx
  with df=n-2, which is what scipy reports as `pvalue`).
"""
import json
import pathlib

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = pathlib.Path(r"REPO")
DATA = ROOT / "radius-valley-age-evolution" / "data"
OUT = ROOT / "manuscripts" / "q1-paper" / "v2-2026-11" / "results" / "audit"


def recompute(tag):
    ages = pd.read_csv(DATA / f"host_ages_{tag}.csv")
    kin = pd.read_csv(DATA / f"hosts_kinematics_{tag}.csv")

    # Selection exactly as specified: gyro_ok True and finite age & v_tan.
    sel = ages[(ages["gyro_ok"] == True)].copy()
    sel["age_gyr"] = pd.to_numeric(sel["age_gyr"], errors="coerce")
    sel = sel[np.isfinite(sel["age_gyr"])]

    m = sel.merge(kin[["hostname", "vtan"]], on="hostname", how="inner")
    m["vtan"] = pd.to_numeric(m["vtan"], errors="coerce")
    m = m[np.isfinite(m["vtan"])]

    r = spearmanr(m["vtan"].values, m["age_gyr"].values)
    n_gyro = int((ages["gyro_ok"] == True).sum())
    cens = int((sel["censored"] == True).sum())
    return {
        "n_hosts_gyro_ok": n_gyro,
        "n": int(len(m)),
        "rho": float(r.statistic),
        "p_exact_scipy": float(r.pvalue),
        "censored_among_gyro_ok": cens,
        "frac_censored_among_gyro_ok": float(cens / n_gyro) if n_gyro else None,
        "n_kinematics_hosts": int(kin["hostname"].nunique()),
        "n_age_table_rows": int(len(ages)),
    }


def main():
    fgk = recompute("FGK")
    m = recompute("M")

    result = {
        "method": ("scipy.stats.spearmanr(vtan, age_gyr) on inner join of "
                   "host_ages_<TAG>.csv and hosts_kinematics_<TAG>.csv; "
                   "selection gyro_ok=True & finite age_gyr & finite vtan"),
        "source_files": {
            "FGK": ["data/host_ages_FGK.csv", "data/hosts_kinematics_FGK.csv"],
            "M": ["data/host_ages_M.csv", "data/hosts_kinematics_M.csv"],
        },
        "FGK": fgk,
        "cool_M": m,
        "manuscript_claims": {
            "FGK": {"rho": 0.50, "n": 79, "p": 3e-6},
            "cool_M": {"rho": 0.25, "n": 12, "p": 0.43},
            "host_count": "12 of 303 M-dwarf hosts (planet_sample_M.csv unique hosts = 303)",
            "censoring": "4 of 12 gyro-ok M hosts censored = 33%",
        },
        "support": {
            "FGK_rho_supported": round(fgk["rho"], 2) == 0.50,
            "FGK_n_supported": fgk["n"] == 79,
            "FGK_p_supported": abs(fgk["p_exact_scipy"] / 3e-6 - 1) < 0.20,
            "M_rho_supported": round(m["rho"], 2) == 0.25,
            "M_n_supported": m["n"] == 12,
            "M_p_supported": round(m["p_exact_scipy"], 2) == 0.43,
            "M_censored_33pct_supported": abs(m["frac_censored_among_gyro_ok"] - 1/3) < 1e-9,
        },
        "python": "C:/Users/MIT/AppData/Local/Programs/Python/Python312/python.exe",
        "note": ("Archived results/proxy_validation_fgk_recomputed.json reports "
                 "rho=0.497794, n=79 and p=4.7e-7 using a NORMAL approximation; "
                 "the exact scipy p-value is 3.05e-6, consistent with the "
                 "manuscript's p=3e-6."),
    }
    out_path = OUT / "recompute_rho_fgk_v2.json"
    with open(out_path, "w") as fh:
        json.dump(result, fh, indent=2)
    print(json.dumps(result, indent=2))
    print(f"\nWROTE {out_path}")


if __name__ == "__main__":
    main()
