"""Audit part 2: UVW convention check (heliocentric vs Galactocentric),
AU Mic sanity check, independent Johnson&Soderblom-style computation."""
import numpy as np
import pandas as pd
from pathlib import Path
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import Galactic

D = str(Path(__file__).resolve().parents[3] / "data") + "/"

g = pd.read_csv(D + "gaia_hosts_M.csv", skiprows=1,
                names=["hostname", "source_id", "ra", "dec", "parallax",
                       "pmra", "pmdec", "radial_velocity", "phot_bp_mean_mag",
                       "phot_rp_mean_mag", "phot_g_mean_mag"])
for c in g.columns[1:]:
    g[c] = pd.to_numeric(g[c], errors="coerce")

row = g[g.hostname == "AU Mic"].iloc[0]
print("AU Mic Gaia row:", dict(row))

c = SkyCoord(ra=row.ra * u.deg, dec=row.dec * u.deg,
             distance=(1000.0 / row.parallax) * u.pc,
             pm_ra_cosdec=row.pmra * u.mas / u.yr,
             pm_dec=row.pmdec * u.mas / u.yr,
             radial_velocity=np.nan_to_num(row.radial_velocity) * u.km / u.s)
gc = c.galactic
U = gc.velocity.d_x.to_value(u.km / u.s)
V = gc.velocity.d_y.to_value(u.km / u.s)
W = gc.velocity.d_z.to_value(u.km / u.s)
mu = np.hypot(row.pmra, row.pmdec); d = 1000.0 / row.parallax
vtan = 4.74 * mu * d / 1000.0
print(f"\nastropy c.galactic.velocity: U={U:.2f} V={V:.2f} W={W:.2f}")
print(f"d_pc={d:.3f} mu_tot={mu:.1f} mas/yr -> vtan={vtan:.2f} km/s")
print(f"vtot=sqrt(U^2+V^2+V^2) = {np.sqrt(U**2+V**2+W**2):.2f}")

# ---- independent manual computation (global Galactic basis, ICRS->Galactic) ----
# Galactic north pole (ICRS): ra_p=192.85948 deg, dec_p=27.12825 deg
# Galactic center direction (ICRS): ra_0=266.40500 deg, dec_0=-28.93617 deg
rap, decp = np.radians(192.85948), np.radians(27.12825)
ra0, dec0 = np.radians(266.40500), np.radians(-28.93617)
r = np.radians(row.ra); dd = np.radians(row.dec)

# unit position vector in ICRS
p = np.array([np.cos(dd)*np.cos(r), np.cos(dd)*np.sin(r), np.sin(dd)])
# GLOBAL Galactic axes expressed in ICRS:
k0 = np.array([np.cos(decp)*np.cos(rap), np.cos(decp)*np.sin(rap), np.sin(decp)])  # W axis -> NGP
i0 = np.array([np.cos(dec0)*np.cos(ra0), np.cos(dec0)*np.sin(ra0), np.sin(dec0)])  # U axis -> GC
j0 = np.cross(k0, i0)                                                              # V axis -> rotation

# velocity components in equatorial cartesian (km/s)
vrad = float(row.radial_velocity)
# tangential speed: v[km/s] = 4.74 * mu[arcsec/yr] * d[pc]
v_ra  = 4.74 * row.pmra*1e-3 * d    # km/s along +RA
v_dec = 4.74 * row.pmdec*1e-3 * d   # km/s along +Dec
e_ra  = np.array([-np.sin(r), np.cos(r), 0.0])
e_dec = np.array([-np.sin(dd)*np.cos(r), -np.sin(dd)*np.sin(r), np.cos(dd)])
vec = vrad*p + v_ra*e_ra + v_dec*e_dec
Ui, Vi, Wi = vec @ i0, vec @ j0, vec @ k0
print(f"independent rotation-matrix calc : U={Ui:.2f} V={Vi:.2f} W={Wi:.2f}")
print(f"max diff astropy vs manual: {max(abs(U-Ui),abs(V-Vi),abs(W-Wi)):.4f} km/s")

# ---- heliocentric vs galactocentric discrimination ----
# If these were GALACTOCENTRIC they'd need solar motion removed
# (U_sun~+11, V_sun~+12(+rot), W_sun~+7.5). A thin-disk old star has
# heliocentric |W|~5-15; Galactocentric W would cluster near 0 +- few km/s.
print("\nl,b of AU Mic:", gc.l.deg, gc.b.deg)
print("Interpretation: astropy Galactic frame origin = Sun (heliocentric).")

# quick population stats: W distribution of ALL RV stars (heliocentric)
ok = g.parallax.notna() & (g.parallax > 0) & g.radial_velocity.notna()
gg = g[ok]
cc = SkyCoord(ra=gg.ra.values*u.deg, dec=gg.dec.values*u.deg,
              distance=(1000/gg.parallax.values)*u.pc,
              pm_ra_cosdec=gg.pmra.values*u.mas/u.yr,
              pm_dec=gg.pmdec.values*u.mas/u.yr,
              radial_velocity=np.nan_to_num(gg.radial_velocity.values)*u.km/u.s)
galc = cc.galactic
Ws = galc.velocity.d_z.to_value(u.km/u.s)
Vs = galc.velocity.d_y.to_value(u.km/u.s)
Us = galc.velocity.d_x.to_value(u.km/u.s)
print(f"\nAll-RV-star heliocentric velocities n={len(gg)}:")
print(f"  U: med {np.median(Us):6.2f} std {Us.std():5.2f}")
print(f"  V: med {np.median(Vs):6.2f} std {Vs.std():5.2f}")
print(f"  W: med {np.median(Ws):6.2f} std {Ws.std():5.2f}")
print("(heliocentric W median should be ~ -solar W ~ -7; "
      "Galactocentric W median would be ~0)")
