"""
audit3_gyrointerp.py — reproduce a few host_ages entries with the exact code
path of 03_gyro_ages.py, probe edge cases, quantify expected censoring.
"""
import pathlib
import numpy as np
import pandas as pd
from gyrointerp.gyro_posterior import gyro_age_posterior

D = pathlib.Path(__file__).resolve().parents[3] / "data"
GRID = np.linspace(0, 3000, 500) / 1000.0

def age_from_posterior(prot, teff):
    if not np.isfinite(prot) or not np.isfinite(teff):
        return np.nan, np.nan, np.nan
    if not (3800.0 <= teff <= 6200.0) or prot > 45.0:
        return np.nan, np.nan, np.nan
    pdf = gyro_age_posterior(prot, teff, Prot_err=0.05 * prot, Teff_err=100.0, n=192)
    pdf = np.asarray(pdf, float)
    if not np.all(np.isfinite(pdf)) or pdf.sum() <= 0:
        return np.nan, np.nan, np.nan
    cdf = np.cumsum(pdf) / pdf.sum()
    q = lambda p: GRID[np.searchsorted(cdf, p)]
    return q(0.5), q(0.16), q(0.84)

agM = pd.read_csv(D / "host_ages_M.csv")
print("--- reproducibility spot checks (M sample) ---")
for host in ["K2-284", "TOI-544", "Ross 176", "TOI-1238", "TOI-4342"]:
    row = agM[agM.hostname == host].iloc[0]
    med, lo, hi = age_from_posterior(row.prot, row.teff)
    print(f"{host:9s} P={row.prot:6.3f} Teff={row.teff:6.0f}: "
          f"csv=({row.age_gyr:.4f},{row.age_lo:.4f},{row.age_hi:.4f}) "
          f"repro=({med:.4f},{lo:.4f},{hi:.4f}) match={np.isclose(med,row.age_gyr)}")

print("\n--- suspicious-case probes (would-be inputs, ignoring validity gate) ---")
for prot, teff in [(44.0, 3800.0), (40.0, 4089.0), (33.0, 4041.0),
                   (14.691, 3866.0), (8.88, 4140.0)]:
    try:
        pdf = np.asarray(gyro_age_posterior(prot, teff, Prot_err=0.05 * prot,
                                            Teff_err=100.0, n=192), float)
        cdf = np.cumsum(pdf) / pdf.sum()
        q = lambda p: GRID[min(np.searchsorted(cdf, p), len(GRID) - 1)]
        frac_at_top = 1.0 - cdf[-2] if cdf[-1] == 1 else np.nan
        print(f"P={prot:7.3f} Teff={teff:6.0f}: median={q(.5):.3f} "
              f"[{q(.16):.3f},{q(.84):.3f}] Gyr; posterior mass at grid top "
              f"(>=2.994): {max(0.0, 1-cdf[-2] if False else (pdf[-1]/pdf.sum())):.3f}")
    except Exception as e:
        print(f"P={prot} Teff={teff}: model error: {type(e).__name__}: {e}")

# below-3800 behaviour (documented as NaN?)
try:
    out = gyro_age_posterior(30.0, 3600.0, Prot_err=1.5, Teff_err=100.0, n=64)
    print("\nTeff=3600 direct call ->", np.nanmin(out), np.nanmax(out),
          "all-nan:", bool(np.all(np.isnan(out))))
except Exception as e:
    print("\nTeff=3600 direct call raised:", type(e).__name__, e)

# expected censoring: which usable M hosts does the slow sequence place above
# the 3 Gyr grid even before truncation? Use the fraction of posterior mass
# beyond 2.94 Gyr (grid top region) as the censoring indicator.
agok = agM[agM.gyro_ok]
print("\n--- usable M hosts: posterior mass piled at grid top ---")
for _, r in agok.sort_values("prot").iterrows():
    pdf = np.asarray(gyro_age_posterior(r.prot, r.teff, Prot_err=0.05 * r.prot,
                                        Teff_err=100.0, n=192), float)
    mass_top = pdf[GRID >= 2.94].sum() / pdf.sum()
    flag = "CEILING" if mass_top > 0.3 else ""
    print(f"{r.hostname:9s} P={r.prot:6.2f} T={r.teff:6.0f}: mass>2.94Gyr = {mass_top:5.2f} {flag}")
