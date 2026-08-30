"""
12_completeness.py — First-order inverse-detection-efficiency weights (IDEM).

REUSED public products (Task 1 requirement):
  - Kepler-host planets (~2135 FGK + ~126 M planets): REAL per-star DR25
    6-hr CDPP values (rrmscdpp06p0) from the Exoplanet Archive
    `keplerstellar` table via KOI kepid join -> dr25_cdpp_matched.csv.
    This is the field-standard public completeness input (Burke &
    Catanzarite 2017 products).
  - Gillis+26 (arXiv:2602.23364): injection-recovery performed for their
    mid-to-late M TESS sample, but no machine-readable sensitivity product
    is linked from the abstract page and their stellar sample only overlaps
    our M hosts partially -> NOT directly reusable here; noted, not used.

BUILT from scratch (first-order, coarse -- stated assumptions):
  - non-Kepler hosts: CDPP proxy from magnitudes,
      Kepler-band fallback: cdpp6 = 185 ppm * 10^(0.20*(Kp-12))
      V-band proxy:          cdpp6 = 260 ppm * 10^(0.22*(V-12))
    (order-of-magnitude scalings consistent with published Kepler/TESS
    noise curves; deliberately conservative)
  - transit duration t_dur = 13 hr * (P_yr)^(1/3) * (rho_sun/rho_star)^(1/3)
    with rho from st_mass/st_rad when present else solar;
  - n_transits = baseline/P, baseline: Kepler 372 d, K2 80 d, TESS 60 d;
  - SNR = depth/cdpp_eff * sqrt(n_transits * t_dur/6 h);
  - P_det = logistic(SNR; midpoint 7.5, width 1.5)  [DR25-average-like],
    floored at 0.02; weight w = 1/P_det normalised to mean 1.

Output: data/completeness_weights.csv (pl_name, p_det, w, src)
"""
import pathlib
import numpy as np
import pandas as pd

D = pathlib.Path(__file__).resolve().parents[1] / "data"
LOGISTIC_MID, LOGISTIC_W = 7.5, 1.5
BASELINE = {"Kepler": 372.0, "K2": 80.0, "TOI": 60.0}


def cdpp_for_host(hostname, kepmag, vmag, cdpp_map):
    fam = hostname.split("-")[0] if "-" in hostname else hostname[:3]
    if hostname in cdpp_map:                       # REUSED: real DR25 value
        return float(cdpp_map[hostname]), "DR25_cdpp"
    if pd.notna(kepmag):
        return 185.0 * 10 ** (0.20 * (kepmag - 12)), "proxy_kepmag"
    if pd.notna(vmag):
        return 260.0 * 10 ** (0.22 * (vmag - 12)), "proxy_vmag"
    return 400.0, "proxy_default"                  # worst-case floor


def p_detect(pl_rade, st_rad, pl_orbper, cdpp_ppm, baseline_d, st_mass,
             st_rad_for_rho=None):
    rstar = np.clip(st_rad, 0.1, 8)
    depth = (pl_rade * 6371.0 / (rstar * 695700.0)) ** 2        # unitless
    p_yr = np.clip(pl_orbper, 0.05, None) / 365.25
    rho_ratio = 1.0                                             # solar default
    if st_mass is not None and np.isfinite(st_mass):
        rho_ratio = np.clip(st_mass, 0.08, 3) / rstar ** 3
    t_dur_hr = 13.0 * p_yr ** (1 / 3) * rho_ratio ** (-1 / 3)
    n_tr = np.maximum(baseline_d / np.clip(pl_orbper, 0.05, None), 1.5)
    snr = depth / (cdpp_ppm * 1e-6) * np.sqrt(n_tr * t_dur_hr / 6.0)
    return 1 / (1 + np.exp(-(snr - LOGISTIC_MID) / LOGISTIC_W))


def main(tag):
    pl = pd.read_csv(D / f"planet_sample_{tag}.csv")
    mags = pd.read_csv(D / f"mags_{tag}.csv").drop_duplicates("pl_name")
    cd = pd.read_csv(D / "dr25_cdpp_matched.csv") if (D / "dr25_cdpp_matched.csv").exists() else None
    cdpp_map = {}
    if cd is not None:
        host_cdpp = cd.groupby("host", as_index=False).rrmscdpp06p0.median()
        cdpp_map = dict(zip(host_cdpp.host, host_cdpp.rrmscdpp06p0))

    rows = []
    for _, r in pl.iterrows():
        km = mags.loc[mags.pl_name == r.pl_name]
        kepmag = km.sy_kepmag.iloc[0] if len(km) else np.nan
        vmag = km.sy_vmag.iloc[0] if len(km) else np.nan
        cval, src = cdpp_for_host(r.hostname, kepmag, vmag, cdpp_map)
        base = next((v for k, v in BASELINE.items()
                     if r.hostname.startswith(k)), 60.0)
        mass = getattr(r, "st_mass", np.nan)
        p = p_detect(r.pl_rade, getattr(r, "st_rad", 0.6), r.pl_orbper,
                     cval, base, mass)
        rows.append({"pl_name": r.pl_name, "p_det": max(p, 0.02),
                     "src": src})
    w = pd.DataFrame(rows)
    w["w"] = (1.0 / w.p_det)
    w["w"] /= w.w.mean()                    # mean-weight-1 convention
    out = D / f"completeness_weights_{tag}.csv"
    w.to_csv(out, index=False)
    print(f"[{tag}] weights: N={len(w)}, median P_det={w.p_det.median():.3f}, "
          f"weight range {w.w.min():.2f}-{w.w.max():.2f}")
    print(f"  sources: {w.src.value_counts().to_dict()}")
    return w


if __name__ == "__main__":
    main("M")
    main("FGK")
