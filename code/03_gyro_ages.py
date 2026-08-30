"""
03_gyro_ages.py — Assign rotation-based ages to planet-hosting stars.

Method: Angus et al. (2019) slow-sequence model via `gyrointerp`.
Model validity is Teff 3800-6200 K ENFORCED HERE -- below 3800 K the model
returns NaN (the package docstring's ~3100 K claim does not apply to this
model path; see models.py L534). Ages are right-censored at the 3 Gyr grid.

Input hygiene:
  - archive st_rotperr1 used where finite/positive, clipped to [5%, 35%]
    of Prot (guards against both overconfident and absurd errors);
    fallback 5% otherwise.
  - stars outside Teff validity or with extreme Prot are flagged unusable
    WITHOUT calling the model.
  - host stellar parameters aggregated deterministically (median Teff)
    so multi-planet hosts with inconsistent archive rows cannot silently
    pick arbitrary values.

Output: data/host_ages_{M,FGK}.csv, written incrementally so partial
progress survives a kill. Column `censored` marks ceiling-piled ages.
"""
import pathlib
import numpy as np
import pandas as pd
from gyrointerp.gyro_posterior import gyro_age_posterior

D = pathlib.Path(__file__).resolve().parents[1] / "data"
GRID = np.linspace(0, 3000, 500) / 1000.0   # Gyr

# HARD REALITY CHECK (verified against models.py source):
# gyrointerp slow sequence is only defined for Teff 3800-6200 K.
# Stars cooler than 3800 K return NaN -> flagged unusable here.
# Consequence: archive-rotp gyro ages cover only ~13 M hosts.
# To recover the cool half we must either adopt literature rotation
# catalogs (McQuillan+14, Newton+16/18) with their own age relations,
# measure Prot ourselves from light curves (Gaidos+24 approach),
# or rely on kinematic ages (05_kinematic_ages.py).
TEFF_LO, TEFF_HI = 3800.0, 6200.0           # model validity window
PROT_MAX = 45.0                             # older than grid anyway -> censored


def age_from_posterior(prot, teff, prot_err=None):
    """Return (median, lo, hi) in Gyr; NaNs if unusable."""
    if not np.isfinite(prot) or not np.isfinite(teff):
        return np.nan, np.nan, np.nan
    if not (TEFF_LO <= teff <= TEFF_HI) or prot > PROT_MAX:
        return np.nan, np.nan, np.nan
    # archive error when sane, clipped to [5%, 35%] of Prot; else 5%
    pe = prot_err
    if not (np.isfinite(pe) and pe > 0):
        pe = 0.05 * prot
    pe = float(np.clip(pe, 0.05 * prot, 0.35 * prot))
    pdf = gyro_age_posterior(prot, teff, Prot_err=pe,
                             Teff_err=100.0, n=192)
    pdf = np.asarray(pdf, dtype=float)
    if not np.all(np.isfinite(pdf)) or pdf.sum() <= 0:
        return np.nan, np.nan, np.nan
    cdf = np.cumsum(pdf) / pdf.sum()
    q = lambda p: GRID[min(np.searchsorted(cdf, p), len(GRID) - 1)]
    return q(0.5), q(0.16), q(0.84)


def main(sample_file, out_name):
    df = pd.read_csv(D / sample_file)
    # deterministic host aggregation: median Teff guards against archive
    # rows disagreeing within a multi-planet system (e.g. HD 63433)
    hosts = (df.dropna(subset=["st_rotp"])
               .groupby("hostname", as_index=False)
               .agg(st_teff=("st_teff", "median"),
                    st_rotp=("st_rotp", "first"),
                    st_rotperr1=("st_rotperr1", "mean")))

    outpath = D / out_name
    rows = []
    for i, h in hosts.iterrows():
        med, lo, hi = age_from_posterior(float(h.st_rotp), float(h.st_teff),
                                         prot_err=h.st_rotperr1)
        rows.append({"hostname": h.hostname, "teff": h.st_teff,
                     "prot": h.st_rotp, "age_gyr": med,
                     "age_lo": lo, "age_hi": hi,
                     "censored": bool(np.isfinite(hi) and hi >= 2.95),
                     "gyro_ok": bool(np.isfinite(med))})
        if i % 10 == 0:
            print(f"[{out_name}] {i+1}/{len(hosts)} done", flush=True)
        pd.DataFrame(rows).to_csv(outpath, index=False)   # incremental save

    ok = pd.DataFrame(rows)
    ok = ok[ok.gyro_ok]
    print(f"\n{sample_file}: {len(hosts)} hosts w/ rotation -> {len(ok)} usable")
    print(f"  median age {ok.age_gyr.median():.2f} Gyr "
          f"(range {ok.age_gyr.min():.2f}-{ok.age_gyr.max():.2f})")
    print(f"  at/near 3 Gyr ceiling: {ok.censored.mean():.0%}")


if __name__ == "__main__":
    main("planet_sample_M.csv", "host_ages_M.csv")
    main("planet_sample_FGK.csv", "host_ages_FGK.csv")
