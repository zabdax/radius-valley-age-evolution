"""
audit5_mcquillan.py — Task 6: verify 05's parser drops units/dashes rows,
verify Kepler-1646 match genuineness.
"""
import pathlib
import numpy as np
import pandas as pd

D = pathlib.Path(__file__).resolve().parents[3] / "data"

# ---- replicate 05 parser exactly ----
rows = []
with open(D / "mcquillan_cool.csv", encoding="utf-8", errors="replace") as f:
    header = None
    raw_noncomment = 0
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        if header is None:
            header = [p.strip() for p in parts]
            continue
        rows.append(parts)
        raw_noncomment += 1
mcq = pd.DataFrame(rows, columns=header)
print(f"header parsed as: {header}")
print(f"data rows before numeric coercion: {len(mcq)}")
print("first 3 raw rows:")
print(mcq.head(3).to_string(index=False))
for c in ["Teff", "Prot", "e_Prot", "_RA", "_DE"]:
    mcq[c] = pd.to_numeric(mcq[c], errors="coerce")
n_before_drop = len(mcq)
mcq = mcq.dropna(subset=["Prot", "_RA", "_DE"])
print(f"rows dropped by coercion+dropna: {n_before_drop - len(mcq)} "
      f"(expect 2: units row + dashes row)")
print(f"final catalog rows: {len(mcq)}; Teff range ({mcq.Teff.min():.0f},{mcq.Teff.max():.0f}); "
      f"Prot range ({mcq.Prot.min():.3f},{mcq.Prot.max():.3f})")
# confirm no leftover junk
junk = mcq[(mcq.Prot <= 0.1) | (mcq.Prot > 70)]
print(f"implausible remaining Prot rows: {len(junk)}")

# does KIC column survive as string junk anywhere?
if "KIC" in mcq.columns:
    bad_kic = pd.to_numeric(mcq["KIC"], errors="coerce").isna().sum()
    print(f"non-numeric KIC left after dropna: {bad_kic}")

# ---- Kepler-1646 match verification ----
m = pd.read_csv(D / "planet_sample_M.csv")
h = m[m.hostname == "Kepler-1646"].iloc[0]
ra_h, dec_h = h.ra, h.dec
print(f"\nKepler-1646 NEA position: RA={ra_h:.7f} Dec={dec_h:.7f} Teff_NEA={h.st_teff}")

c_ra, c_dec = float(mcq._RA.values[0]), None  # placeholder
# nearest McQuillan source to host position
d2 = ((mcq._RA - ra_h) * np.cos(np.radians(dec_h))) ** 2 + (mcq._DE - dec_h) ** 2
i = int(d2.idxmin())
row = mcq.loc[i]
dra = (row._RA - ra_h) * np.cos(np.radians(dec_h)) * 3600
ddec = (row._DE - dec_h) * 3600
sep = np.hypot(dra, ddec)
print(f"nearest McQuillan row: KIC={int(row.KIC)} Teff={row.Teff:.0f} Prot={row.Prot:.3f}"
      f"+-{row.e_Prot:.3f} sep={sep:.3f} arcsec (dRA*cos={dra:.3f}, dDec={ddec:.3f})")

hosts = pd.read_csv(D / "hosts_M_with_mcquillan.csv")
kh = hosts[hosts.hostname == "Kepler-1646"].iloc[0]
print(f"05 output record: mcq_prot={kh.mcq_prot}, mcq_sep={kh.mcq_sep:.4f}, "
      f"st_rotp={kh.st_rotp} (NaN -> genuinely new)")

# second-nearest neighbor for confusion check
d2s = d2.drop(i)
j = int(d2s.idxmin())
r2 = mcq.loc[j]
sep2 = np.hypot((r2._RA - ra_h) * np.cos(np.radians(dec_h)) * 3600,
                (r2._DE - dec_h) * 3600)
print(f"second-nearest McQuillan source: KIC={int(r2.KIC)} sep={sep2:.1f} arcsec")

# plausibility: Prot 3.514 d at Teff 3299 K.
# Rossby number check: R_o = Prot / tau_c (Guenther+22-ish tau ~ 60-70 d at 3300K)
tau = 70.0
Ro = 3.514 / tau
print(f"\nplausibility: Ro = Prot/tau(~{tau:.0f}d at 3299K) = {Ro:.3f} "
      f"(saturated fast-rotator regime Ro<~0.13 typical of young field M dwarfs)")
sat = mcq[(mcq.Teff < 3400) & (mcq.Prot < 5)]
print(f"McQuillan context: {len(sat)} catalog stars with Teff<3400 & Prot<5d "
      f"-> short-P cool rotators are common in this catalog")
frac_short = (mcq[mcq.Teff < 3500].Prot < 10).mean()
print(f"fraction of Teff<3500 catalog stars with Prot<10d: {frac_short:.0%}")
