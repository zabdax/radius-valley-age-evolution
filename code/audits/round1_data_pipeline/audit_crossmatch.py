"""Audit part 3: (a) has_rv masking / nan_to_num leakage into stored columns;
(b) 06 nearest-neighbour separation stats, epoch propagation, companion checks;
(c) host-loss accounting sample(303) -> gaia_hosts(?) -> kinematics(301);
(d) 06 resume-read crash reproduction."""
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord

D = "C:/Users/MIT/Downloads/Astro_research_ASGSR/mdwarf_radius_valley/data/"
names = ["hostname", "source_id", "ra", "dec", "parallax", "pmra", "pmdec",
         "radial_velocity", "phot_bp_mean_mag", "phot_rp_mean_mag",
         "phot_g_mean_mag"]

# raw read replicating 07.load_gaia semantics
first = open(D + "gaia_hosts_M.csv").readline().strip()
skip = 1 if first.split(",")[0] == "0" else 0
print(f"gaia_hosts_M first line starts with '{first[:12]}...' -> skiprows={skip}")
g = pd.read_csv(D + "gaia_hosts_M.csv", names=names, skiprows=skip)
for c in names[1:]:
    g[c] = pd.to_numeric(g[c], errors="coerce")
print("gaia rows:", len(g), "| unique hostnames:", g.hostname.nunique(),
      "| dup hostnames:", g.hostname.duplicated().sum(),
      "| dup source_id:", g.source_id.duplicated().sum())

k = pd.read_csv(D + "hosts_kinematics_M.csv")
m = pd.read_csv(D + "planet_sample_M.csv")

# ---------- (a) has_rv mask verification ----------
print("\n=== has_rv MASK CHECK ===")
print("stored kinematics rows:", len(k),
      "| W notna:", k.W.notna().sum(), "| has_rv True:", k.has_rv.sum(),
      "| vtot notna:", k.vtot.notna().sum())
bad_w = k[k.W.notna() & ~k.has_rv]
bad_nw = k[k.W.isna() & k.has_rv]
print("rows with stored W but has_rv False:", len(bad_w))
print("rows with has_rv True but W NaN   :", len(bad_nw))
chk = k.merge(g[["hostname", "radial_velocity"]], on="hostname", how="left")
print("rows with stored W but gaia RV NaN:",
      int((chk.W.notna() & chk.radial_velocity.isna()).sum()))
print("rows with W NaN but gaia RV present:",
      int((chk.W.isna() & chk.radial_velocity.notna()).sum()))

# re-run compute_kinematics fresh and compare every number
ok = g.parallax.notna() & (g.parallax > 0) & g.pmra.notna() & g.pmdec.notna()
gg = g[ok].copy()
c = SkyCoord(ra=gg.ra.values*u.deg, dec=gg.dec.values*u.deg,
             distance=(1000/gg.parallax.values)*u.pc,
             pm_ra_cosdec=gg.pmra.values*u.mas/u.yr,
             pm_dec=gg.pmdec.values*u.mas/u.yr,
             radial_velocity=np.nan_to_num(gg.radial_velocity.values)*u.km/u.s)
mu = np.hypot(gg.pmra.values, gg.pmdec.values); dpc = 1000/gg.parallax.values
fresh = pd.DataFrame({"hostname": gg.hostname.values,
                      "vtan": 4.74*mu*dpc/1000.0})
has_rv = gg.radial_velocity.notna().values
gc_ = c.galactic
vgal = np.array([gc_.velocity.d_x.to_value(u.km/u.s),
                 gc_.velocity.d_y.to_value(u.km/u.s),
                 gc_.velocity.d_z.to_value(u.km/u.s)])
U, V, W = vgal
fresh["has_rv"] = has_rv
fresh["W"] = np.where(has_rv, W, np.nan)
fresh["vtot"] = np.where(has_rv,
                         np.sqrt(U**2+V**2+W**2), np.nan)
cmp = k.merge(fresh, on="hostname", suffixes=("_stored", "_fresh"))
print("recompute vs stored: rows compared:", len(cmp),
      "| vtan max|diff|: %.2e" % np.abs(cmp.vtan_stored-cmp.vtan_fresh).max())
both = cmp[cmp.W_stored.notna()]
print("| W max|diff| over RV stars: %.2e | vtot max|diff|: %.2e"
      % (np.abs(both.W_stored-both.W_fresh).max(),
         np.abs(both.vtot_stored-both.vtot_fresh).max()))
print("has_rv flags identical:", bool((cmp.has_rv_stored == cmp.has_rv_fresh).all()))
# dropped-by-parallax accounting
dropped = g[~g.hostname.isin(fresh.hostname)]
print("gaia rows dropped by parallax/pm filter:",
      dropped.hostname.tolist(), "parallax:", dropped.parallax.tolist())

# ---------- (c) host-loss accounting ----------
print("\n=== HOST ACCOUNTING ===")
shosts = set(m.hostname.unique()); ghosts = set(g.hostname)
khosts = set(k.hostname)
print("sample hosts:", len(shosts), "| gaia_hosts:", len(ghosts),
      "| kinematics:", len(khosts))
print("hosts w/o gaia match (cone miss or run incomplete):",
      sorted(shosts - ghosts)[:20], f"(n={len(shosts-ghosts)})")

# ---------- (b) separation analysis ----------
print("\n=== SEPARATION ANALYSIS ===")
# replicate 06's per-host archive position choice
hosts = m.sort_values("pl_rade").groupby("hostname", as_index=False).first()[
    ["hostname", "st_teff", "ra", "dec", "sy_plx", "sy_pmra", "sy_pmdec",
     "sy_dist"]]
j = g.merge(hosts, on="hostname", suffixes=("_gaia", "_arch"))
cg = SkyCoord(j.ra_gaia.values*u.deg, j.dec_gaia.values*u.deg)
ca = SkyCoord(j.ra_arch.values*u.deg, j.dec_arch.values*u.deg)
sep_true = cg.separation(ca).arcsec
dra = (j.ra_gaia - j.ra_arch) * np.cos(np.radians(j.dec_arch))
ddec = j.dec_gaia - j.dec_arch
sep_flat = np.hypot(dra, ddec) * 3600.0
j["sep_flat"] = sep_flat; j["sep_true"] = sep_true
print("flat-sky vs exact separation max abs diff: %.2e arcsec"
      % np.abs(sep_flat-sep_true).max())
print("sep [arcsec]: min %.4f  med %.4f  p90 %.4f  p99 %.4f  max %.4f"
      % (j.sep_true.min(), j.sep_true.median(), j.sep_true.quantile(.9),
         j.sep_true.quantile(.99), j.sep_true.max()))
for t in [0.1, 0.25, 0.5, 1.0]:
    print(f"  hosts with sep > {t}\": {(j.sep_true > t).sum()} "
          f"({100*(j.sep_true>t).mean():.1f}%)")

# epoch propagation expectation: archive pos (often J2000) vs Gaia J2016
mu_tot = np.hypot(j.pmra, j.pmdec)
dt = 2016.0 - 2000.0          # years, order-of-magnitude epoch gap
exp_shift = mu_tot * dt / 1000.0
j["mu_tot"] = mu_tot; j["exp_shift"] = exp_shift
big = j[j.sep_true > 0.25].sort_values("sep_true", ascending=False)
cols = ["hostname", "sep_true", "mu_tot", "exp_shift", "phot_g_mean_mag"]
print("\ntop-15 largest separations:")
print(big[cols].head(15).to_string(index=False))
corr = j[["sep_true", "mu_tot"]].corr().iloc[0, 1]
print(f"\nPearson corr(sep, mu_total) = {corr:.3f}")
print(f"hosts whose sep < expected epoch-shift*1.5 (consistent w/ epoch drift): "
      f"{int((j.sep_true < j.exp_shift*1.5).sum())}/{len(j)}")

# companion / wrong-star indicators: archive vs Gaia parallax & PM consistency
# (archive sy_plx/sy_pm* refer to the HOST; a matched visual companion or an
#  unrelated background source would disagree)
j["f_plx"] = np.abs(j.parallax - j.sy_plx) / j.sy_plx
j["f_pm"] = np.hypot(j.pmra - j.sy_pmra, j.pmdec - j.sy_pmdec) / np.hypot(
    j.sy_pmra, j.sy_pmdec)
print("\nparallax consistency |plx_gaia-sy_plx|/sy_plx:")
print("  med %.4f  p90 %.4f  max %.4f" % (j.f_plx.median(), j.f_plx.quantile(.9),
                                         j.f_plx.max()))
print("PM consistency |dPMvec|/|PMvec|:")
print("  med %.4f  p90 %.4f  max %.4f" % (j.f_pm.median(), j.f_pm.quantile(.9),
                                          j.f_pm.max()))
sus = j[(j.f_plx > 0.25) | (j.f_pm > 0.25)].sort_values("sep_true",
                                                        ascending=False)
cols = ["hostname", "sep_true", "phot_g_mean_mag", "sy_plx", "parallax",
        "sy_pmra", "pmra", "sy_pmdec", "pmdec"]
print(f"\nSUSPICIOUS matches (plx or PM inconsistent >25%): n={len(sus)}")
if len(sus):
    print(sus[cols].to_string(index=False))
else:
    print("none")

# absolute-mag plausibility: M dwarfs span M_G ~ +6..+13
j["MG"] = j.phot_g_mean_mag - 5*np.log10(j.sy_dist/10.0)
print("\nabsolute G of matched sources: med %.2f  range %.2f..%.2f"
      % (j.MG.median(), j.MG.min(), j.MG.max()))
weird = j[(j.MG < 5) | (j.MG > 14)]
print("implausible-M_G matches (M_G<5 or >14):", weird.hostname.tolist())

# ---------- (d) resume-read crash reproduction ----------
print("\n=== 06 RESUME-READ REPRODUCTION ===")
try:
    pd.read_csv(D + "gaia_hosts_M.csv").hostname
    print("pd.read_csv(OUT).hostname OK")
except Exception as e:
    print("06 resume path CRASHES on legacy file: %s: %s"
          % (type(e).__name__, e))
