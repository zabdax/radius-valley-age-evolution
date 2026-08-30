"""
08_headline_fgk_vtan.py — The powered FGK positive-control analysis.

Reproduces every number in the ROADMAP 'HEADLINE FINDING' section:
  - FGK x v_tan Fisher split (all hosts)
  - metallicity-restricted splits
  - distance/local-volume cut
  - mission (hostname-prefix) localization
Run after 06_gaia_crossmatch.py planet_sample_FGK.csv and
07_kinematic_ages.py FGK.
"""
import pathlib
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

D = pathlib.Path(__file__).resolve().parents[1] / "data"
VALLEY_F = 1.88


def fisher(y, o):
    t = [[int(y.is_sn.sum()), int((~y.is_sn).sum())],
         [int(o.is_sn.sum()), int((~o.is_sn).sum())]]
    odds, p = fisher_exact(t)
    return y.is_sn.mean(), o.is_sn.mean(), p


def main():
    pl = pd.read_csv(D / "planet_sample_FGK.csv")
    k = pd.read_csv(D / "hosts_kinematics_FGK.csv")
    df = pl.merge(k[["hostname", "vtan"]], on="hostname", how="inner")
    df["is_sn"] = df.pl_rade > VALLEY_F

    def split(d, label):
        med = d.vtan.median()
        y, o = d[d.vtan <= med], d[d.vtan > med]
        sy, so, p = fisher(y, o)
        print(f"{label:<42} N={len(d):5d}  SNy={sy:.3f} SNo={so:.3f}  p={p:.4f}")

    print("=== FGK x vtan (positive control) ===")
    split(df, "ALL FGK")
    split(df[df.st_met.abs() < 0.10], "|Fe/H| < 0.10")
    split(df[df.st_met.abs() < 0.05], "|Fe/H| < 0.05")
    split(df[df.sy_dist < 200], "dist < 200 pc")
    split(df[(df.st_met.abs() < 0.10) & (df.sy_dist < 200)],
          "solar-Fe/H & dist<200")

    print("\n=== mission localization ===")
    df["mission"] = df.hostname.str.extract(
        r"^(TOI|Kepler|KIC|K2|HIP|HD|GJ|Wolf|LHS|TRAPPIST|EPIC)",
        expand=False).fillna("other")
    for msub, g in df.groupby("mission"):
        if len(g) >= 60:
            split(g, f"mission={msub}")

    lo, hi = df[df.vtan <= df.vtan.median()], df[df.vtan > df.vtan.median()]
    print(f"\nmedian distance: low-v {lo.sy_dist.median():.0f} pc, "
          f"high-v {hi.sy_dist.median():.0f} pc")


if __name__ == "__main__":
    main()
