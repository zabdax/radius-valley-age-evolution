"""
audit1_valley.py — Tasks 1 & 2: valley-floor robustness + KDE method check.
READ-ONLY audit; writes nothing outside results/audit_tmp2/.
"""
import pathlib
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

D = pathlib.Path("C:/Users/MIT/Downloads/Astro_research_ASGSR/mdwarf_radius_valley/data")

m = pd.read_csv(D / "planet_sample_M.csv")
f = pd.read_csv(D / "planet_sample_FGK.csv")
print(f"M sample n={len(m)}, FGK n={len(f)}")

# ---------- replicate 02_explore algorithm exactly ----------
def valley_floor(rades, sigma, lo_mode_hi=1.70, sn_lo=1.70, grid_lo=0.9, grid_hi=3.6, ngrid=200):
    r = np.sort(np.asarray(rades, float))
    grid = np.linspace(grid_lo, grid_hi, ngrid)
    kde = np.exp(-0.5 * ((grid[:, None] - r[None, :]) / sigma) ** 2).sum(axis=1)
    w_se = (grid > 0.9) & (grid < lo_mode_hi)
    w_sn = (grid > sn_lo) & (grid < 3.2)
    mode_se = grid[w_se][np.argmax(kde[w_se])]
    mode_sn = grid[w_sn][np.argmax(kde[w_sn])]
    between = (grid > mode_se) & (grid < mode_sn)
    vmin = grid[between][np.argmin(kde[between])]
    return mode_se, mode_sn, vmin

print("\n=== TASK 1a: baseline reproduction (sigma=0.12, boundary 1.70) ===")
mse, msn, vm = valley_floor(m.pl_rade.values, 0.12)
fse, fsn, vf = valley_floor(f.pl_rade.values, 0.12)
print(f"M  : modes {mse:.4f}/{msn:.4f} -> floor {vm:.4f} (reported 1.85)")
print(f"FGK: floor {vf:.4f} (reported 1.88)")

print("\n=== TASK 1b: sigma x boundary sweep ===")
print(f"{'sigma':>6} | {'b=1.65':>8} {'b=1.70':>8} {'b=1.75':>8}   (M floor)")
rows = []
for sig in [0.08, 0.10, 0.12, 0.15, 0.20]:
    line = []
    for b in [1.65, 1.70, 1.75]:
        _, _, v = valley_floor(m.pl_rade.values, sig, lo_mode_hi=b, sn_lo=b)
        line.append(v); rows.append(("M", sig, b, v))
    print(f"{sig:>6} | " + " ".join(f"{x:8.3f}" for x in line))
print()
for sig in [0.08, 0.10, 0.12, 0.15, 0.20]:
    line = []
    for b in [1.65, 1.70, 1.75]:
        _, _, v = valley_floor(f.pl_rade.values, sig, lo_mode_hi=b, sn_lo=b)
        line.append(v); rows.append(("FGK", sig, b, v))
    print(f"FGK sig={sig}: " + " ".join(f"{x:7.3f}" for x in line))

# finer look at the M KDE minimum shape at sigma=0.12: how sharp is the min?
r = np.sort(m.pl_rade.values)
grid = np.linspace(0.9, 3.6, 200)
kde = np.exp(-0.5 * ((grid[:, None] - r[None, :]) / 0.12) ** 2).sum(axis=1)
between = (grid > 1.23) & (grid < 2.08)
kb = kde[between]; gb = grid[between]
i0 = np.argmin(kb)
print("\nM-sample KDE around minimum (sigma=0.12):")
for di in range(-4, 5):
    j = i0 + di
    if 0 <= j < len(gb):
        mark = " <-- min" if di == 0 else ""
        print(f"  R={gb[j]:.4f}  dens={kb[j]:9.2f}{mark}")
print(f"relative dip depth vs neighbors +-2 grid pts: "
      f"(kde[i-2]-min)/max={((kb[i0-2]-kb[i0])/kb[i0]):.4f}, (kde[i+2]-min)/max={((kb[i0+2]-kb[i0])/kb[i0]):.4f}")

# ---------- TASK 2: truncation / range sensitivity ----------
print("\n=== TASK 2: KDE range & truncation checks ===")
def floor_with_subset(rades, sigma=0.12, subset=None, grid_lo=0.9, grid_hi=3.6, ngrid=200):
    r_all = np.sort(np.asarray(rades, float))
    r_use = np.sort(np.asarray(rades if subset is None else subset, float))
    grid = np.linspace(grid_lo, grid_hi, ngrid)
    # kernels from ALL planets (as in 02)
    kde_all = np.exp(-0.5 * ((grid[:, None] - r_all[None, :]) / sigma) ** 2).sum(axis=1)
    # kernels only from planets inside the plotted window
    kde_sub = np.exp(-0.5 * ((grid[:, None] - r_use[None, :]) / sigma) ** 2).sum(axis=1)
    def fl(k):
        w_se = (grid > 0.9) & (grid < 1.7); w_sn = (grid > 1.7) & (grid < 3.2)
        ms = grid[w_se][np.argmax(k[w_se])]; Mn = grid[w_sn][np.argmax(k[w_sn])]
        btw = (grid > ms) & (grid < Mn)
        return ms, Mn, grid[btw][np.argmin(k[btw])]
    return fl(kde_all), fl(kde_sub)

grid = np.linspace(0.9, 3.6, 200)
for name, dfx in [("M", m), ("FGK", f)]:
    rr = dfx.pl_rade.values
    inside = rr[(rr >= 0.9) & (rr <= 3.6)]
    # all-data kernels (as 02 does)
    kall = np.exp(-0.5 * ((grid[:, None] - np.sort(rr)[None, :]) / 0.12) ** 2).sum(axis=1)
    # in-window-only kernels
    ksub = np.exp(-0.5 * ((grid[:, None] - np.sort(inside)[None, :]) / 0.12) ** 2).sum(axis=1)
    def _floor(k):
        w_se = (grid > 0.9) & (grid < 1.7); w_sn = (grid > 1.7) & (grid < 3.2)
        ms = grid[w_se][np.argmax(k[w_se])]; Ms = grid[w_sn][np.argmax(k[w_sn])]
        btw = (grid > ms) & (grid < Ms)
        return grid[btw][np.argmin(k[btw])]
    v_all, v_in = _floor(kall), _floor(ksub)
    n_out = len(rr) - len(inside)
    print(f"{name}: n_outside_window={n_out}; floor all-data={v_all:.4f}, "
          f"in-window-only={v_in:.4f}, diff={100*(v_in-v_all)/v_all:+.2f}%")
    # max tail contribution of an out-of-window point AT the floor location:
    dmin = min(abs(v_all - rr.min()), abs(rr.max() - v_all))
    print(f"   nearest data edge to floor: {dmin:.2f} R_E away -> kernel weight "
          f"exp(-0.5*(d/0.12)^2) <= {np.exp(-0.5*(dmin/0.12)**2):.2e}")

# grid-range choice sensitivity
print("\nGrid range sensitivity (sigma=0.12, boundary 1.70):")
for glo, ghi, ng in [(0.9, 3.6, 200), (0.5, 4.0, 300), (1.0, 3.2, 200), (0.9, 3.6, 400)]:
    for name, dfx in [("M", m), ("F", f)]:
        _, _, v = valley_floor(dfx.pl_rade.values, 0.12, grid_lo=glo, grid_hi=ghi, ngrid=ng)
        print(f"  grid[{glo},{ghi}]n{ng} {name}: floor={v:.4f}")

# normalized-vs-unnormalized check (does missing 1/(sigma sqrt(2pi)) matter?)
print("\nNormalization constant irrelevant for argmin (constant factor per fixed sigma) -- verified algebraically.")

# ---------- TASK 1c: downstream SE/SN counts + Fisher ----------
print("\n=== TASK 1c: downstream impact of floor choice ===")
ag = pd.read_csv(D / "host_ages_M.csv")
ok_hosts = ag.loc[ag.gyro_ok, "hostname"]
mm = m.merge(ag[["hostname", "age_gyr", "gyro_ok"]], on="hostname")
mm = mm[mm.gyro_ok & (mm.age_gyr > 0.05)]
print(f"gyro_ok M hosts: {len(ok_hosts)}; usable planets after age_gyr>0.05 cut: {len(mm)} "
      f"({mm.hostname.nunique()} hosts)")

for floor in [1.75, 1.85, 1.95]:
    se = int(((m.pl_rade > 0.5) & (m.pl_rade <= floor)).sum())
    sn = int(((m.pl_rade > floor) & (m.pl_rade <= 4.0)).sum())
    sub = mm.copy(); sub["is_sn"] = sub.pl_rade > floor
    med = sub.age_gyr.median()
    y = sub[sub.age_gyr <= med]; o = sub[sub.age_gyr > med]
    tab = [[int(y.is_sn.sum()), int((~y.is_sn).sum())],
           [int(o.is_sn.sum()), int((~o.is_sn).sum())]]
    odds, p = fisher_exact(tab)
    print(f"\nfloor={floor}: full-M SE={se} SN={sn} (SN frac {sn/(se+sn):.3f})")
    print(f"  gyro subsample N={len(sub)}: young SN {y.is_sn.mean():.3f} ({tab[0][0]}/{sum(tab[0])}), "
          f"old SN {o.is_sn.mean():.3f} ({tab[1][0]}/{sum(tab[1])})")
    print(f"  Fisher odds={odds:.3f} p={p:.3f}; direction 'old<yOUNG SN frac' "
          f"(Gaidos-consistent): {o.is_sn.mean() < y.is_sn.mean()}")

# also show what fraction of planets move across the cut between floors
band = ((m.pl_rade > 1.75) & (m.pl_rade <= 1.95)).sum()
print(f"\nplanets in (1.75,1.95] reclassified when floor moves 1.75->1.95: {band} "
      f"({100*band/len(m):.1f}% of M sample)")
