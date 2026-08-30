"""
05_literature_rotations.py — Recover rotation periods for hosts missing them.

Step 1: cross-match M-dwarf planet hosts against McQuillan et al. 2014
        (Kepler rotation periods, J/ApJS/211/24/table1) by position.
Reports how many cool (<3800 K) hosts gain a Prot -> those become
gyro-eligible once paired with an M-dwarf-appropriate age relation.
"""
import pathlib
import subprocess
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u

D = pathlib.Path(__file__).resolve().parents[1] / "data"
MCQ = D / "mcquillan_cool.csv"


def ensure_mcquillan():
    """Fetch McQuillan+2014 cool rotators (J/ApJS/211/24/table1) via the
    classic VizieR asu-tsv service if not already on disk. The TAP mirror
    403s on this query; asu-tsv is the reliable route."""
    if MCQ.exists() and MCQ.stat().st_size > 10000:
        return
    url = "https://vizier.cds.unistra.fr/viz-bin/asu-tsv"
    cmd = ["curl", "-sL", "-G", "--max-time", "240", url,
           "--data-urlencode", "-source=J/ApJS/211/24/table1",
           "--data-urlencode", "-out=KIC,Teff,Prot,e_Prot,_RA,_DE",
           "--data-urlencode", "-out.max=50000",
           "--data-urlencode", "Teff=<4200",
           "--data-urlencode", "Prot=>0",
           "-o", str(MCQ)]
    subprocess.run(cmd, check=True)
    print(f"downloaded {MCQ} ({MCQ.stat().st_size} bytes)")


ensure_mcquillan()

# --- load McQuillan TSV (skip comment header block) ---
rows = []
with open(D / "mcquillan_cool.csv", encoding="utf-8", errors="replace") as f:
    header = None
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        if header is None:
            header = [p.strip() for p in parts]
            continue
        rows.append(parts)
mcq = pd.DataFrame(rows, columns=header)
for c in ["Teff", "Prot", "e_Prot", "_RA", "_DE"]:
    mcq[c] = pd.to_numeric(mcq[c], errors="coerce")
mcq = mcq.dropna(subset=["Prot", "_RA", "_DE"])
print(f"McQuillan cool rotators loaded: {len(mcq)}")

# --- load planet sample, one row per host ---
m = pd.read_csv(D / "planet_sample_M.csv")
hosts = m.sort_values("pl_rade").groupby("hostname", as_index=False).first()[
    ["hostname", "st_teff", "st_rotp", "ra", "dec"]]
print(f"M hosts total: {len(hosts)}, already have archive Prot: {hosts.st_rotp.notna().sum()}")

c_hosts = SkyCoord(ra=hosts.ra.values * u.deg, dec=hosts.dec.values * u.deg)
c_mcq = SkyCoord(ra=mcq["_RA"].values * u.deg, dec=mcq["_DE"].values * u.deg)

idx, d2d, _ = c_hosts.match_to_catalog_sky(c_mcq)
sep = d2d.arcsec
hit = sep < 5.0

hosts["mcq_prot"] = np.where(hit, mcq.Prot.values[idx], np.nan)
hosts["mcq_sep"] = np.where(hit, sep, np.nan)

need = hosts.st_rotp.isna()
gained = need & hit
print(f"cross-matches <5 arcsec      : {hit.sum()}")
print(f"new Prot recovered for hosts : {gained.sum()} "
      f"({(gained & (hosts.st_teff<3800)).sum()} of them Teff<3800)")
print("\nRecovered hosts:")
cols = ["hostname", "st_teff", "st_rotp", "mcq_prot", "mcq_sep"]
print(hosts.loc[gained, cols].to_string(index=False))

# sanity: agreement where both exist
both = hosts.st_rotp.notna() & hit
if both.sum():
    ratio = (hosts.mcq_prot[both] / hosts.st_rotp[both]).round(2)
    print(f"\nagreement check (archive/mcq) n={both.sum()}: median ratio "
          f"{np.median(ratio):.2f}")

hosts.to_csv(D / "hosts_M_with_mcquillan.csv", index=False)
