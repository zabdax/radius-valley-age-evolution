"""
audit4_tap_counts.py — Task 5: fresh independent COUNT(*) queries against
NEA TAP + local duplicate checks.
"""
import pathlib
import urllib.parse
import urllib.request
import pandas as pd

BASE = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
D = pathlib.Path("C:/Users/MIT/Downloads/Astro_research_ASGSR/mdwarf_radius_valley/data")

CUTS = ("default_flag=1 AND pl_controv_flag=0 "
        "AND discoverymethod='Transit' "
        "AND pl_rade>0.5 AND pl_rade<4 AND pl_orbper<100")

QUERIES = {
    "M":   f"SELECT COUNT(*) AS n FROM ps WHERE {CUTS} AND st_teff<4200",
    "FGK": f"SELECT COUNT(*) AS n FROM ps WHERE {CUTS} AND st_teff>=4200 AND st_teff<=7000",
}

for name, q in QUERIES.items():
    url = BASE + "?" + urllib.parse.urlencode({"query": q, "format": "csv"})
    with urllib.request.urlopen(url, timeout=120) as r:
        body = r.read().decode()
    print(f"TAP fresh count {name}: {body.strip().splitlines()[-1]}")

for fname, expected in [("planet_sample_M.csv", 419), ("planet_sample_FGK.csv", 2880)]:
    df = pd.read_csv(D / fname)
    dup_pl = df.pl_name.duplicated().sum()
    dup_full = df.duplicated().sum()
    print(f"{fname}: rows={len(df)} (expected {expected}), "
          f"dup pl_name={dup_pl}, fully-dup rows={dup_full}, "
          f"unique hosts={df.hostname.nunique()}, teff range=({df.st_teff.min()},{df.st_teff.max()})")
    if dup_pl:
        print(df[df.pl_name.duplicated(keep=False)].sort_values("pl_name")[
            ["pl_name", "hostname", "pl_rade", "st_teff"]].to_string(index=False))
