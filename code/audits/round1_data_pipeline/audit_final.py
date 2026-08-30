"""Audit part 4: Test A details, gyro-validate repro, FGK spot checks,
pandas header behavior, groupby.first() risk, merge-loss identity."""
import io
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, spearmanr

D = "C:/Users/MIT/Downloads/Astro_research_ASGSR/mdwarf_radius_valley/data/"

# --- pandas names/header behavior (07.load_gaia latent hazard) ---
buf = "hostname,num\nAU Mic,1\n"
df = pd.read_csv(io.StringIO(buf), names=["hostname", "num"], skiprows=0)
print("read_csv(names=...) on headed file -> rows:", len(df),
      "first hostname cell:", repr(df.hostname.iloc[0]))

# --- host_ages_M details ---
ag = pd.read_csv(D + "host_ages_M.csv")
print("\nhost_ages_M: rows", len(ag), "| gyro_ok True:", int(ag.gyro_ok.sum()))
print(ag[ag.gyro_ok].age_gyr.describe())
ok_ages = ag[ag.gyro_ok]
print("gyro_ok ages <= 0.05:", int((ok_ages.age_gyr <= 0.05).sum()),
      "| < 0.1:", int((ok_ages.age_gyr < 0.1).sum()),
      "| min:", ok_ages.age_gyr.min())
print("AU Mic age row:", ag[ag.hostname == "AU Mic"].to_dict("records"))

# --- 07.validate gyro-Spearman reproduction ---
k = pd.read_csv(D + "hosts_kinematics_M.csv")
mg = k.merge(ok_ages, on="hostname")
okm = mg.age_gyr > 0.05
rho, p = spearmanr(np.log10(mg.age_gyr[okm]), mg.vtan[okm])
print(f"\ngyro-vs-vtan Spearman n={okm.sum()} rho={rho:.2f} p={p:.3f} "
      "(results file says rho=0.12 p=0.712 n=12)")

# --- planets lost in 04 Test-C merge ---
m = pd.read_csv(D + "planet_sample_M.csv")
lost = m[~m.hostname.isin(k.hostname)]
print("\nplanets dropped by Test C merge:",
      lost[["pl_name", "hostname"]].to_string(index=False))

# --- groupby.first() mixed-row risk in 06 ---
sub = m[["hostname", "ra", "dec"]]
nan_ra = sub.ra.isna().sum(); nan_dec = sub.dec.isna().sum()
g_n = sub.groupby("hostname").nunique()
mixed = g_n[(g_n.ra > 1) | (g_n.dec > 1)]
print(f"\nra NaN={nan_ra}, dec NaN={nan_dec} in planet_sample_M")
print("hosts with >1 distinct ra or dec:", len(mixed))
print(mixed.head())

# --- FGK spot checks (Test B + FGK vtan headline in ROADMAP) ---
VALLEY_F = 1.88
f = pd.read_csv(D + "planet_sample_FGK.csv"); f["is_sn"] = f.pl_rade > VALLEY_F
kf = pd.read_csv(D + "hosts_kinematics_FGK.csv")
agf = pd.read_csv(D + "host_ages_FGK.csv")
ff = f.merge(agf[["hostname", "age_gyr", "gyro_ok"]], on="hostname")
ffb = ff[ff.gyro_ok & (ff.age_gyr > 0.05)]
dropped_f = ff[ff.gyro_ok & ~(ff.age_gyr > 0.05)]
print(f"\nFGK Test B: N={len(ffb)} (roadmap says 128); "
      f"drops {len(dropped_f)} planets by age>0.05 filter")

def split(df, col):
    df = df.dropna(subset=[col])
    med = df[col].median()
    y, o = df[df[col] <= med], df[df[col] > med]
    tab = [[int(y.is_sn.sum()), int((~y.is_sn).sum())],
           [int(o.is_sn.sum()), int((~o.is_sn).sum())]]
    odds, p = fisher_exact(tab)
    print(f"split@{med:.2f}: SN low {y.is_sn.mean():.3f} high {o.is_sn.mean():.3f}"
          f" odds={odds:.2f} p={p:.4f}")

print("Test B (gyro age):"); split(ffb, "age_gyr")
mgf = f.merge(kf[["hostname", "vtan", "W"]], on="hostname")
print("FGK x vtan (roadmap: young 0.561 old 0.604 p=0.023):")
split(mgf, "vtan")
w = kf.W.notna()
print(f"FGK C2 signed-W:"); split(mgf[mgf.W.notna()].assign(), "W")
tmp = mgf.copy(); tmp["absW"] = tmp.W.abs()
print("FGK C2 |W|-corrected:"); split(tmp[tmp.W.notna()], "absW")

# FGK anchor matching contamination check
ANCHORS = {"HIP 67522": 17, "TOI 1227": 8, "K2-25": 650, "K2-136": 550}
print("\nFGK anchor matching (contains-token):")
for host, age in ANCHORS.items():
    row = kf[kf.hostname.str.contains(host.split()[0], case=False, na=False)]
    if len(row):
        print(f"{host:<10} token='{host.split()[0]:<8}' hits={len(row):<4} "
              f"picked='{row.iloc[0].hostname}' vtan={row.iloc[0].vtan:.1f}")
print("true HIP 67522 row:",
      kf[kf.hostname == "HIP 67522"][["vtan"]].to_dict("records"))
print("true K2-136 row:",
      kf[kf.hostname == "K2-136"][["vtan"]].to_dict("records"))
