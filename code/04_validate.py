"""
04_validate.py — Validation gate (v2): does our pipeline reproduce known results?

Test A (Gaidos+24 reproduction): young-vs-old sub-Neptune fraction decline
        around M hosts using rotation-based gyro ages (Teff>=3800 subsample).
Test B (FGK positive control): identical test on FGK hosts with gyro ages.
Test C (kinematic proxy sanity): young-vs-old SN-fraction split on the
        kinematic velocity ranking for ALL M hosts with astrometry.

The science claim is NOT made here. This file only checks that directions and
magnitudes are consistent with published work before we build on them.
"""
import pathlib
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

D = pathlib.Path(__file__).resolve().parents[1] / "data"
VALLEY_M, VALLEY_F = 1.85, 1.88


def fisher_split(df, agecol, label):
    df = df.dropna(subset=[agecol]).copy()
    if len(df) < 8:
        print(f"[{label}] N={len(df)} too small for a split test")
        return None
    med = df[agecol].median()
    y = df[df[agecol] <= med]
    o = df[df[agecol] > med]
    tab = [[int(y.is_sn.sum()), int((~y.is_sn).sum())],
           [int(o.is_sn.sum()), int((~o.is_sn).sum())]]
    odds, p = fisher_exact(tab)
    tag = "" if len(df) >= 20 else "  ** N<20: SUGGESTIVE ONLY **"
    print(f"\n--- {label} ---{tag}")
    print(f"N={len(df)} planets / {df.hostname.nunique()} hosts | "
          f"split at {med:.2f} {agecol}")
    print(f"SN frac young: {y.is_sn.mean():.3f} ({tab[0][0]}/{tab[0][0]+tab[0][1]})")
    print(f"SN frac old  : {o.is_sn.mean():.3f} ({tab[1][0]}/{tab[1][0]+tab[1][1]})")
    print(f"Fisher odds={odds:.2f} p={p:.3f} | direction matches Gaidos: "
          f"{o.is_sn.mean() < y.is_sn.mean()}")
    return p


def base(sample_file, valley):
    pl = pd.read_csv(D / sample_file)
    pl["is_sn"] = pl.pl_rade > valley
    return pl


if __name__ == "__main__":
    # Test A: M hosts, gyro ages. NOTE: no blanket age floor -- genuinely
    # young (floor-piled) hosts are real anchors, not artifacts.
    m = base("planet_sample_M.csv", VALLEY_M)
    try:
        ag = pd.read_csv(D / "host_ages_M.csv")
        mm = m.merge(ag[["hostname", "age_gyr", "gyro_ok"]], on="hostname")
        mm = mm[mm.gyro_ok]
        pA = fisher_split(mm, "age_gyr", "A: M hosts x gyro age")
    except FileNotFoundError:
        print("host_ages_M.csv missing")

    # Test B: FGK control
    try:
        f = base("planet_sample_FGK.csv", VALLEY_F)
        agf = pd.read_csv(D / "host_ages_FGK.csv")
        ff = f.merge(agf[["hostname", "age_gyr", "gyro_ok"]], on="hostname")
        ff = ff[ff.gyro_ok]
        pB = fisher_split(ff, "age_gyr", "B: FGK control x gyro age")
    except FileNotFoundError:
        print("host_ages_FGK.csv missing")

    # Test C: M hosts, kinematic velocity ranking (vtan as ordinal time)
    try:
        k = pd.read_csv(D / "hosts_kinematics_M.csv")
        mk = m.merge(k[["hostname", "vtan", "W"]], on="hostname")
        pC = fisher_split(mk, "vtan",
                          "C: M hosts x vtan ranking (high v = older)")
        # vertical-velocity ranking MUST use |W|: signed W is symmetric about
        # the solar value and a median split would rank by sign, not energy.
        if mk.W.notna().sum() >= 8:
            mk["absW"] = mk.W.abs()
            pD = fisher_split(mk, "absW", "C2: M hosts x |W| ranking (RV subset)")
    except FileNotFoundError:
        print("hosts_kinematics_M.csv missing - run 07 first")
