"""
07_kinematic_ages.py — Velocity-based age indicators for every host.

Inputs : data/gaia_hosts_M.csv   (from 06_gaia_crossmatch.py)
         data/planet_sample_M.csv

Outputs: data/hosts_kinematics_M.csv  (one row per host)
         results/kinematic_validation.txt

Proxies:
  vtan  = 4.74 * mu_total[mas/yr] * d[pc]          (needs no RV -> all hosts)
  U,V,W Galactocentric space velocities             (needs Gaia RV)
Age usage:
  v1 (this file): ORDINAL ranking only (median splits), validated against
     known-young systems and against gyro ages where available.
  v2 (later): quantitative inversion of a published sigma_W(t) law inside
     the hierarchical model, with honest per-star uncertainties.

Validation checks printed:
  1. Known-young hosts must rank among the lowest velocities.
  2. Spearman rank correlation between velocity indicator and
     log(literature/gyro age) over all anchors available.
"""
import pathlib
import sys
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.stats import spearmanr

D = pathlib.Path(__file__).resolve().parents[1] / "data"
R = pathlib.Path(__file__).resolve().parents[1] / "results"
TAG = sys.argv[1] if len(sys.argv) > 1 else "M"   # "M" or "FGK"

# literature ages for anchor systems (Myr) -- well-studied young hosts
ANCHORS = {
    "AU Mic": 22,        # beta Pic mg
    "HIP 67522": 17,
    "DS Tuc": 40,
    "TOI 1227": 8,       # archive name style may vary; checked below
    "K2-25": 650,        # Hyades
    "K2-136": 550,       # Praesepe-ish / Hyades-aged group
}


def load_gaia():
    """Read crossmatch output; tolerate the integer-header written by the
    first (buggy) version of 06 by forcing names and skipping row 1."""
    names = ["hostname", "source_id", "ra", "dec", "parallax", "pmra",
             "pmdec", "radial_velocity", "phot_bp_mean_mag",
             "phot_rp_mean_mag", "phot_g_mean_mag"]
    first = open(D / f"gaia_hosts_{TAG}.csv").readline().strip()
    skip = 1 if first.split(",")[0] == "0" else 0
    g = pd.read_csv(D / f"gaia_hosts_{TAG}.csv", names=names, skiprows=skip)
    for c in names[1:]:
        g[c] = pd.to_numeric(g[c], errors="coerce")
    return g


def compute_kinematics(g):
    ok = g.parallax.notna() & (g.parallax > 0) & g.pmra.notna() & g.pmdec.notna()
    g = g[ok].copy()
    c = SkyCoord(ra=g["ra"].values * u.deg, dec=g["dec"].values * u.deg,
                 distance=(1000.0 / g["parallax"].values) * u.pc,
                 pm_ra_cosdec=g["pmra"].values * u.mas / u.yr,
                 pm_dec=g["pmdec"].values * u.mas / u.yr,
                 radial_velocity=np.nan_to_num(g["radial_velocity"].values,
                                               nan=0.0) * u.km / u.s)
    mu_tot = np.hypot(g["pmra"].values, g["pmdec"].values)
    d_pc = 1000.0 / g["parallax"].values
    out = pd.DataFrame({"hostname": g.hostname,
                        "vtan": 4.74 * mu_tot * d_pc / 1000.0})  # km/s
    has_rv = g["radial_velocity"].notna().values
    gc = c.galactic
    vgal = np.array([gc.velocity.d_x.to_value(u.km / u.s),
                     gc.velocity.d_y.to_value(u.km / u.s),
                     gc.velocity.d_z.to_value(u.km / u.s)])
    U, V, W = vgal   # HELIOCENTRIC Galactic velocities (U->GC, V->rotation,
                     # W->NGP); includes solar motion. NOT Galactocentric --
                     # any future sigma_W(t) inversion must use heliocentric
                     # dispersions to match published calibrations.
    out["has_rv"] = has_rv
    out.loc[has_rv, "W"] = W[has_rv]
    out.loc[has_rv, "vtot"] = np.sqrt(U[has_rv] ** 2 + V[has_rv] ** 2 + W[has_rv] ** 2)
    # fallback total-speed proxy without RV uses vtan only
    return out


def validate(kin):
    lines = []
    norm = lambda s: str(s).strip().replace("-", " ").replace("_", " ").lower()
    kin["_norm"] = kin.hostname.map(norm)
    known = []
    for host, age_myr in ANCHORS.items():
        row = kin[kin._norm == norm(host)]
        if len(row):
            r = row.iloc[0]
            known.append((host, age_myr, r.vtan))
            lines.append(f"{host:<12} lit_age={age_myr:>5} Myr  vtan={r.vtan:6.1f} km/s")
        else:
            lines.append(f"{host:<12} lit_age={age_myr:>5} Myr  NOT IN SAMPLE")
    if len(known) >= 3:
        rho, p = spearmanr([np.log10(a) for _, a, _ in known],
                           [v for _, _, v in known])
        lines.append(f"\nSpearman(log age_literature, vtan) n={len(known)}: "
                     f"rho={rho:.2f} p={p:.3f}")
    else:
        lines.append(f"\nanchor coverage too thin for Spearman (n={len(known)})")
    # vs gyro ages
    try:
        gy = pd.read_csv(D / "host_ages_M.csv")
        gy = gy[gy.gyro_ok]
        mg = kin.merge(gy, on="hostname")
        if len(mg) >= 5:
            ok = mg.age_gyr > 0.05      # exclude floor-piled zeros
            rho2, p2 = spearmanr(np.log10(mg.age_gyr[ok]), mg.vtan[ok])
            lines.append(f"Spearman(log gyro age, vtan)  n={ok.sum()}: "
                         f"rho={rho2:.2f} p={p2:.3f}")
    except FileNotFoundError:
        pass
    return "\n".join(lines)


if __name__ == "__main__":
    g = load_gaia()
    print(f"[{TAG}] gaia rows: {len(g)}, with RV: {g.radial_velocity.notna().sum()}")
    kin = compute_kinematics(g)
    kin.to_csv(D / f"hosts_kinematics_{TAG}.csv", index=False)

    rep = validate(kin)
    print("\n=== VALIDATION ===")
    print(rep)
    open(R / f"kinematic_validation_{TAG}.txt", "w").write(rep)

    vt = kin.vtan.dropna()
    print(f"\nvtan distribution: median {vt.median():.1f} km/s, "
          f"p25 {vt.quantile(.25):.1f}, p75 {vt.quantile(.75):.1f}")
