"""
02_explore.py — First-pass look at the samples.

Questions this answers before any science:
  1. Is the radius valley visible in our M sample? Where does it sit?
  2. How many planets have host rotation periods (gyro ages)?
  3. How many have full astrometry+RV (kinematic ages)?
  4. Super-Earth vs sub-Neptune counts under a simple valley cut.
"""
import pathlib
import numpy as np
import pandas as pd

D = pathlib.Path(__file__).resolve().parents[1] / "data"

m = pd.read_csv(D / "planet_sample_M.csv")
f = pd.read_csv(D / "planet_sample_FGK.csv")
print(f"M sample: {len(m)} planets, {m.hostname.nunique()} hosts")
print(f"FGK     : {len(f)} planets, {f.hostname.nunique()} hosts")

def hist_ascii(vals, lo, hi, nb=40, width=60):
    h, edges = np.histogram(np.clip(vals, lo, hi), bins=np.linspace(lo, hi, nb + 1))
    peak = h.max()
    out = []
    for i, c in enumerate(h):
        bar = "#" * int(width * c / peak)
        out.append(f"{edges[i]:5.2f}-{edges[i+1]:5.2f} | {bar} {c}")
    return "\n".join(out)

print("\n=== Radius distribution, M sample (<4 Rearth) ===")
print(hist_ascii(m.pl_rade.values, 0.5, 4.0))

# valley floor = density minimum STRICTLY BETWEEN the two modes
# (global-min-in-window fails here because the sparse sub-Neptune tail keeps
#  falling toward 2.6+, so the window edge beats the real dip)
r = np.sort(m.pl_rade.values)
grid = np.linspace(0.9, 3.6, 200)
kde = np.exp(-0.5 * ((grid[:, None] - r[None, :]) / 0.12) ** 2).sum(axis=1)

mode_se = grid[(grid > 0.9) & (grid < 1.7)][np.argmax(kde[(grid > 0.9) & (grid < 1.7)])]
mode_sn = grid[(grid > 1.7) & (grid < 3.2)][np.argmax(kde[(grid > 1.7) & (grid < 3.2)])]
between = (grid > mode_se) & (grid < mode_sn)
vmin = grid[between][np.argmin(kde[between])]
# GUARD: the minimum must be an interior dip, not a window-edge artifact
# (at large smoothing widths the sparse sub-Neptune tail falls monotonically
# and argmin lands on the boundary -- verified failure mode at sigma=0.20)
assert mode_se + 0.15 < vmin < mode_sn - 0.15, \
    f"valley floor {vmin:.2f} hit window edge between modes {mode_se:.2f}/{mode_sn:.2f}"
rel_depth = 1 - kde[between].min() / kde[between].max()
print(f"\nmodes: SE {mode_se:.2f} / SN {mode_sn:.2f} Rearth")
print(f"valley floor (between modes): ~{vmin:.2f} Rearth  "
      f"(relative dip depth {rel_depth:.3f}; floor systematic +-0.03 over "
      f"smoothing sigma 0.08-0.15)")

for cut_lo, cut_hi, name in [(0.5, vmin, "super-Earths"), (vmin, 4.0, "sub-Neptunes")]:
    n = ((m.pl_rade > cut_lo) & (m.pl_rade <= cut_hi)).sum()
    print(f"  {name}: {n}")

print("\n=== Data coverage in M sample ===")
has_rot = m.st_rotp.notna() & (m.st_teff < 4200)
print(f"planets w/ host rotation period      : {has_rot.sum()} ({m.loc[has_rot,'hostname'].nunique()} unique hosts)")
rot_warm = has_rot & (m.st_teff >= 3500)
print(f"  ...and Teff>=3500K (gyro-reliable) : {rot_warm.sum()} ({m.loc[rot_warm,'hostname'].nunique()} unique hosts)")
astro = m.sy_plx.notna() & m.sy_pmra.notna() & m.sy_pmdec.notna()
print(f"planets w/ parallax + proper motions : {astro.sum()}")
full = astro & m.st_radv.notna()
print(f"  ...plus radial velocity            : {full.sum()}")

print("\n=== Rotation periods available by Teff bin ===")
rr = m[m.st_rotp.notna()].groupby(pd.cut(m.st_teff[m.st_rotp.notna()], [2500, 3000, 3200, 3500, 3800, 4200]),
                                  observed=True).agg(n_pl=("pl_name", "size"), n_host=("hostname", "nunique"))
print(rr)

print("\n=== Same coverage stats for FGK sample (baseline) ===")
print(f"rotation: {f.st_rotp.notna().sum()} planets ({f.loc[f.st_rotp.notna(),'hostname'].nunique()} hosts)")

# quick FGK valley check for sanity
rf = np.sort(f.pl_rade.values)
kdef = np.exp(-0.5 * ((grid[:, None] - rf[None, :]) / 0.12) ** 2).sum(axis=1)
vmask = (grid > 1.3) & (grid < 2.6)
vflo = grid[vmask][np.argmin(kdef[vmask])]
print(f"FGK valley floor: ~{vflo:.2f} Rearth")

m.to_csv(D / "_cache_M.csv", index=False)
print("\nsaved cache ok")
