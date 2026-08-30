"""Audit part 1: Tests A/C/C2 recompute, merge cardinality, anchor matching."""
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, spearmanr

D = "C:/Users/MIT/Downloads/Astro_research_ASGSR/mdwarf_radius_valley/data/"
VALLEY_M = 1.85

m = pd.read_csv(D + "planet_sample_M.csv")
m["is_sn"] = m.pl_rade > VALLEY_M
k = pd.read_csv(D + "hosts_kinematics_M.csv")

print("=== CARDINALITY CHECKS ===")
print("planet_sample_M rows:", len(m), "unique hosts:", m.hostname.nunique())
print("hosts_kinematics_M rows:", len(k), "unique hosts:", k.hostname.nunique(),
      "| duplicate hostnames in kinematics:", k.hostname.duplicated().sum())
print("kin hosts missing from sample:",
      set(k.hostname) - set(m.hostname))
mk_all = m.merge(k[["hostname", "vtan", "W", "vtot", "has_rv"]], on="hostname")
print("merged rows (sample x kin):", len(mk_all), " (expect =", len(m), ")")
print("has_rv True rows in kinematics:", k.has_rv.sum(), "/", len(k))

def fisher_split(df, agecol, label):
    df = df.dropna(subset=[agecol]).copy()
    if len(df) < 20:
        print(f"[{label}] N={len(df)} too small for a split test")
        return None
    med = df[agecol].median()
    y = df[df[agecol] <= med]
    o = df[df[agecol] > med]
    tab = [[int(y.is_sn.sum()), int((~y.is_sn).sum())],
           [int(o.is_sn.sum()), int((~o.is_sn).sum())]]
    odds, p = fisher_exact(tab)
    print(f"\n--- {label} ---")
    print(f"N={len(df)} planets / {df.hostname.nunique()} hosts | split at {med:.2f} {agecol}")
    print(f"SN frac low : {y.is_sn.mean():.3f} ({tab[0][0]}/{tab[0][0]+tab[0][1]})")
    print(f"SN frac high: {o.is_sn.mean():.3f} ({tab[1][0]}/{tab[1][0]+tab[1][1]})")
    print(f"Fisher odds={odds:.2f} p={p:.4f}")
    return p

# ---- Test C: vtan ranking (reproduce pipeline) ----
pC = fisher_split(mk_all, "vtan", "Test C reproduce: vtan signed split (pipeline)")

# ---- Test C2 AS RUN (signed W) ----
pC2_signed = fisher_split(mk_all, "W", "Test C2 AS RUN: RAW SIGNED W split")

# ---- Test C2 CORRECTED (|W|) ----
mk_all["absW"] = mk_all.W.abs()
pC2_abs = fisher_split(mk_all, "absW", "Test C2 CORRECTED: |W| split")

# extra: what does the signed-W median split actually separate?
w = mk_all.W.dropna()
print("\nW stats: n=%d median=%.2f mean=%.2f std=%.2f min=%.2f max=%.2f"
      % (len(w), w.median(), w.mean(), w.std(), w.min(), w.max()))
medW = mk_all.W.median()
low = mk_all[mk_all.W <= medW]; high = mk_all[mk_all.W > medW]
print("signed-W low side |W| median: %.2f ; high side |W| median: %.2f"
      % (low.W.abs().median(), high.W.abs().median()))
print("=> signed split does NOT order vertical energy: sides' mean |W| = %.2f vs %.2f"
      % (low.W.abs().mean(), high.W.abs().mean()))
# also: does signed-W split equal |W| split at all?
agree = ((mk_all.W <= medW) == (mk_all.W.abs() <= mk_all.absW.median())).mean()
print("fraction of stars where signed-W young-side == |W| young-side: %.2f" % agree)

# ---- Test A both ways ----
ag = pd.read_csv(D + "host_ages_M.csv")
print("\n=== host_ages_M ===")
print("rows:", len(ag), "dup hostnames:", ag.hostname.duplicated().sum())
mm = m.merge(ag[["hostname", "age_gyr", "gyro_ok"]], on="hostname")
mm_f = mm[mm.gyro_ok & (mm.age_gyr > 0.05)]
pA_filtered = fisher_split(mm_f, "age_gyr", "A AS RUN: gyro age > 0.05 filter")
mm_all = mm[mm.gyro_ok]
pA_nofilter = fisher_split(mm_all, "age_gyr", "A NO-FILTER: all gyro_ok ages incl floor")
dropped = mm[mm.gyro_ok & ~(mm.age_gyr > 0.05)]
print("planets dropped by age>0.05 filter:", len(dropped),
      " hosts:", dropped.hostname.nunique())
print("dropped SN frac: %.3f (n=%d)" % (dropped.is_sn.mean(), len(dropped)))
print("dropped hosts:", sorted(dropped.hostname.unique()))

# ---- 07 validate() anchor matching reproduction ----
import sys
sys.path.insert(0, "../code")
ANCHORS = {
    "AU Mic": 22, "HIP 67522": 17, "DS Tuc": 40,
    "TOI 1227": 8, "K2-25": 650, "K2-136": 550,
}
print("\n=== ANCHOR MATCHING REPRODUCTION (07.validate) ===")
known = []
for host, age_myr in ANCHORS.items():
    row = k[k.hostname.str.contains(host.split()[0], case=False, na=False)]
    if len(row):
        r = row.iloc[0]
        known.append((host, age_myr, r.vtan))
        print(f"{host:<12} token='{host.split()[0]:<9}' matched_rows={len(row):<3} "
              f"iloc[0].hostname='{r.hostname}' vtan={r.vtan:.1f}")
    else:
        print(f"{host:<12} token='{host.split()[0]:<9}' NO MATCH")
if len(known) >= 3:
    rho, p = spearmanr([np.log10(a) for _, a, _ in known], [v for _, _, v in known])
    print(f"Spearman(log lit_age, vtan) n={len(known)} rho={rho:.2f} p={p:.3f}")

# correct anchor matching (exact-ish)
print("\ncorrected exact-name matching:")
known2 = []
for host, age_myr in ANCHORS.items():
    hits = k[k.hostname.str.contains(host.replace(" ", "").replace("-", "").lower(),
                                     case=False, regex=False)]
    hits = k[k.hostname.str.lower().str.startswith(host.split()[0].lower()[:3])]
print(k[k.hostname.str.contains("AU Mic|HIP 67522|DS Tuc|TOI-1227|K2-25|K2-136",
                                case=False, regex=True)][["hostname", "vtan", "W", "has_rv"]])
