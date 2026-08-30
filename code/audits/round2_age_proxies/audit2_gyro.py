"""
audit2_gyro.py — Tasks 3 & 4: host_ages_{M,FGK}.csv internal consistency +
st_rotperr1 handling in 03_gyro_ages.py.
"""
import pathlib
import numpy as np
import pandas as pd

D = pathlib.Path("C:/Users/MIT/Downloads/Astro_research_ASGSR/mdwarf_radius_valley/data")

for name, sample in [("M", "planet_sample_M.csv"), ("FGK", "planet_sample_FGK.csv")]:
    print(f"\n================ {name} ================")
    ag = pd.read_csv(D / f"host_ages_{name}.csv")
    pl = pd.read_csv(D / sample)
    ok = ag[ag.gyro_ok]
    bad = ag[~ag.gyro_ok]
    print(f"rows={len(ag)}, gyro_ok={len(ok)}, flagged={len(bad)}")

    # (a) validity windows
    t_out = ok[(ok.teff < 3800) | (ok.teff > 6200)]
    p_out = ok[ok.prot > 45.0]
    print(f"(a) gyro_ok rows outside Teff [3800,6200]: {len(t_out)}; "
          f"gyro_ok rows w/ Prot>45: {len(p_out)}")
    # flagged rows: confirm they fail ONLY for legit reasons
    if len(bad):
        why_t = ((bad.teff < 3800) | (bad.teff > 6200)).sum()
        why_p = (bad.prot > 45.0).sum()
        why_nan = (~np.isfinite(bad.teff) | ~np.isfinite(bad.prot)).sum()
        other = len(bad) - sum(int(x) for x in [why_t, why_p, why_nan])
        print(f"    flagged breakdown: teff-out={why_t}, prot>45={why_p}, "
              f"nan-input={why_nan}, other(model-level)={other}")
        if other:
            print(bad[(bad.teff >= 3800) & (bad.teff <= 6200) & (bad.prot <= 45)
                      & np.isfinite(bad.teff) & np.isfinite(bad.prot)])

    # (b) interval ordering
    v_ok = ag.age_lo.notna()
    viol_lo = ag[v_ok & (ag.age_lo > ag.age_gyr)]
    viol_hi = ag[v_ok & (ag.age_gyr > ag.age_hi)]
    nan_mismatch = ag[ag.gyro_ok != ag.age_gyr.notna()]
    print(f"(b) lo>med violations: {len(viol_lo)}; med>hi violations: {len(viol_hi)}; "
          f"gyro_ok/nan mismatches: {len(nan_mismatch)}")

    # (c)+(d) inspect usable set
    o = ok.sort_values("prot")
    print("(c/d) usable hosts sorted by Prot:")
    print(o[["hostname", "teff", "prot", "age_gyr", "age_lo", "age_hi"]].to_string(index=False))
    ceil = (o.age_hi >= 2.95).mean()
    ceil_med = (o.age_gyr >= 2.9).mean()
    print(f"at/near ceiling: age_hi>=2.95: {ceil:.0%}; age_gyr>=2.9: {ceil_med:.0%}")

    # monotonicity: within +-400 K Teff bands, is age increasing in Prot?
    viol = 0; comps = 0
    for i in range(len(o)):
        for j in range(len(o)):
            if i < j and abs(o.teff.iloc[i] - o.teff.iloc[j]) <= 400:
                comps += 1
                if o.prot.iloc[i] < o.prot.iloc[j] and not (
                        o.age_hi.iloc[i] >= o.age_lo.iloc[j]):   # intervals overlap -> not a violation
                    if o.age_gyr.iloc[i] > o.age_gyr.iloc[j]:
                        viol += 1
    print(f"(c) pairwise same-band (dTeff<=400K) Prot-ordering checks: {comps}, "
          f"clear non-overlapping violations: {viol}")

    # cross-check against planet sample: do the CSV hosts match archive rotp hosts?
    rot_hosts = set(pl.loc[pl.st_rotp.notna(), "hostname"].unique())
    csv_hosts = set(ag.hostname)
    print(f"archive rot-hosts={len(rot_hosts)} vs csv={len(csv_hosts)}; "
          f"missing={sorted(rot_hosts-csv_hosts)[:5]} extra={sorted(csv_hosts-rot_hosts)[:5]}")
    # per-host st_rotp/st_teff consistency across multiple planet rows
    g = pl.dropna(subset=["st_rotp"]).groupby("hostname").agg(
        nt=("st_teff", "nunique"), np_=("st_rotp", "nunique"))
    incons = g[(g.nt > 1) | (g.np_ > 1)]
    print(f"hosts with inconsistent st_teff/st_rotp across planet rows: {len(incons)}")

# ---------- TASK 4: archive st_rotperr1 vs assumed 0.05*Prot ----------
print("\n================ TASK 4: st_rotperr1 vs assumed 5% ================")
for name, sample in [("M", "planet_sample_M.csv"), ("FGK", "planet_sample_FGK.csv")]:
    pl = pd.read_csv(D / sample)
    ag = pd.read_csv(D / f"host_ages_{name}.csv")
    h = pl.dropna(subset=["st_rotp"]).sort_values("pl_rade").groupby(
        "hostname", as_index=False).first()[["hostname", "st_teff", "st_rotp", "st_rotperr1"]]
    h["assumed"] = 0.05 * h.st_rotp
    have = h.st_rotperr1.notna() & (h.st_rotperr1 > 0)
    print(f"\n{name}: rot hosts={len(h)}, with positive archive err={have.sum()} "
          f"(NaN/absent: {len(h)-have.sum()})")
    hh = h[have].copy()
    hh["ratio"] = hh.st_rotperr1 / hh.assumed
    wild = hh[(hh.ratio > 2) | (hh.ratio < 0.5)].sort_values("ratio")
    print(f"  ratio archive_err/(5%Prot): median={hh.ratio.median():.2f}, "
          f"min={hh.ratio.min():.2f}, max={hh.ratio.max():.2f}")
    print(wild[["hostname", "st_teff", "st_rotp", "st_rotperr1", "assumed", "ratio"]].to_string(index=False))

# ---------- identify the FGK host with inconsistent per-row star params ----------
print("\n================ inconsistent multi-row host (FGK) ================")
plf = pd.read_csv(D / "planet_sample_FGK.csv")
gf = plf.dropna(subset=["st_rotp"]).groupby("hostname").agg(
    nt=("st_teff", "nunique"), np_=("st_rotp", "nunique"))
bad_h = gf[(gf.nt > 1) | (gf.np_ > 1)].index.tolist()
for bh in bad_h:
    cols = ["pl_name", "pl_rade", "st_teff", "st_rotp", "st_rotperr1"]
    print(plf.loc[plf.hostname == bh, cols].to_string(index=False))

