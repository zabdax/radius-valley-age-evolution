"""
12b_sharpen_weights.py — Option (ii): sharpen non-Kepler P_det model.

Why: the first-order weights gave median P_det = 0.21 (TOI) / 0.17 (K2)
for *detected* planets around M dwarfs -- implausibly low, traced to
proxying TESS/K2 noise from V-band magnitudes. M dwarfs are much brighter
at redder bands (V-T colour ~ 2 mag), so V-based noise wildly overestimates
the difficulty of their transits.

Sharpening (documented assumptions):
  - T = G - 0.5 : standard Gaia-G -> TESS-band approximation (valid to
    ~0.05 mag for most cool stars; G from our Gaia DR3 crossmatch).
  - TESS 6-hr CDPP from log-linear interpolation over published-noise
    anchors (Sullivan+15 / Instrument Handbook order):
      T =   7 ->  40 ppm     T = 10 -> 130 ppm
      T = 12 -> 300 ppm      T = 14 -> 900 ppm
    (ASSUMPTION - flagged for sensitivity testing.)
  - Kepler mag fallback for K2: Kp = G - 0.2; K2 noise = Kepler DR25-style
    magnitude scaling x 2.0 roll/jitter factor.
  - Nearby-named hosts without Gaia G (none expected): fall back to old
    proxy, flagged.

Output: data/completeness_weights_M_sharp.csv
"""
import pathlib
import numpy as np
import pandas as pd

D = pathlib.Path(__file__).resolve().parents[1] / "data"

# log-linear TESS noise anchors: (Tmag, 6-hr CDPP ppm)
ANCH_T = np.array([7.0, 10.0, 12.0, 14.0])
ANCH_C = np.log10(np.array([40.0, 130.0, 300.0, 900.0]))


def tess_cdpp6(t):
    return 10 ** float(np.interp(t, ANCH_T, ANCH_C))


def kep_cdpp6(kp):
    return 185.0 * 10 ** (0.20 * (kp - 12))          # as in 12_completeness


def p_detect(pl_rade, st_rad, pl_orbper, cdpp_ppm, baseline_d):
    rstar = np.clip(st_rad, 0.1, 8)
    depth = (pl_rade * 6371.0 / (rstar * 695700.0)) ** 2
    p_yr = np.clip(pl_orbper, 0.05, None) / 365.25
    t_dur_hr = 13.0 * p_yr ** (1 / 3)
    n_tr = np.maximum(baseline_d / np.clip(pl_orbper, 0.05, None), 1.5)
    snr = depth / (cdpp_ppm * 1e-6) * np.sqrt(n_tr * t_dur_hr / 6.0)
    return 1 / (1 + np.exp(-(snr - 7.5) / 1.5))


def main():
    pl = pd.read_csv(D / "planet_sample_M.csv")
    g = pd.read_csv(D / "gaia_hosts_M.csv")           # has phot_g_mean_mag
    gm = g.drop_duplicates("hostname").set_index("hostname")
    w_old = pd.read_csv(D / "completeness_weights_M.csv")

    rows = []
    for _, r in pl.iterrows():
        G = gm.phot_g_mean_mag.get(r.hostname, np.nan)
        if r.hostname.startswith("TOI"):
            if not np.isfinite(G):
                src, cval, base = "fallback_v", 400.0, 60.0
            else:
                T = G - 0.5
                src, cval, base = "sharpen_TESS_G", tess_cdpp6(T), 60.0
        elif r.hostname.startswith(("K2", "EPIC")):
            if not np.isfinite(G):
                src, cval, base = "fallback_v_k2", 370.0, 80.0
            else:
                kp = G - 0.2
                src, cval, base = "sharpen_K2_G", 2.0 * kep_cdpp6(kp), 80.0
        else:
            # Kepler/nearby keep previous treatment
            old = w_old.loc[w_old.pl_name == r.pl_name]
            rows.append({"pl_name": r.pl_name,
                         "w": float(old.w.iloc[0]) if len(old) else 1.0,
                         "src": "unchanged"})
            continue
        mass = getattr(r, "st_mass", np.nan)
        rstar = np.clip(getattr(r, "st_rad", 0.6) or 0.6, 0.1, 8)
        p = p_detect(r.pl_rade, rstar, r.pl_orbper, cval, base)
        rows.append({"pl_name": r.pl_name, "w": 1.0 / max(p, 0.02),
                     "src": src})
    out = pd.DataFrame(rows)
    out["w"] /= out.w.mean()
    merged = pl[["pl_name"]].merge(out, on="pl_name", how="left")
    merged["w"] = merged.w.fillna(1.0)
    merged.to_csv(D / "completeness_weights_M_sharp.csv", index=False)

    print("source mix:", merged.src.value_counts().to_dict() if "src" in merged
          else "n/a")
    print(f"weights: median {merged.w.median():.2f}, "
          f"p90 {np.quantile(merged.w,.9):.2f}, max {merged.w.max():.2f}")
    ess = merged.w.sum() ** 2 / (merged.w ** 2).sum()
    print(f"full-sample ESS after sharpening: {ess:.0f} / {len(merged)}")


if __name__ == "__main__":
    main()
