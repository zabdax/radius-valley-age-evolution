import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u

D = "C:/Users/MIT/Downloads/Astro_research_ASGSR/mdwarf_radius_valley/data/"
names = ["hostname", "source_id", "ra", "dec", "parallax", "pmra", "pmdec",
         "radial_velocity", "phot_bp_mean_mag", "phot_rp_mean_mag",
         "phot_g_mean_mag"]
g = pd.read_csv(D + "gaia_hosts_M.csv", names=names, skiprows=1)
for c in names[1:]:
    g[c] = pd.to_numeric(g[c], errors="coerce")
m = pd.read_csv(D + "planet_sample_M.csv")
hosts = m.sort_values("pl_rade").groupby("hostname", as_index=False).first()[
    ["hostname", "st_teff", "st_rad", "ra", "dec"]]
j = g.merge(hosts, on="hostname", suffixes=("_gaia", "_arch"))
j["sep"] = SkyCoord(j.ra_gaia.values*u.deg, j.dec_gaia.values*u.deg).separation(
    SkyCoord(j.ra_arch.values*u.deg, j.dec_arch.values*u.deg)).arcsec

bad = j[j.hostname.isin(["Kepler-1410", "Kepler-1652"])]
print(bad[["hostname", "source_id", "sep", "parallax", "pmra", "pmdec",
           "radial_velocity", "phot_g_mean_mag", "st_teff"]].to_string(index=False))

# brightest/faintest matched sources overall (G mags)
print("\nG mag range of all matches: %.2f..%.2f" %
      (j.phot_g_mean_mag.min(), j.phot_g_mean_mag.max()))
faintest = j.nlargest(5, "phot_g_mean_mag")[
    ["hostname", "sep", "phot_g_mean_mag", "parallax"]]
print(faintest.to_string(index=False))
