"""
01_download.py — Pull the planet samples from the NASA Exoplanet Archive (TAP).

Samples:
  M  : confirmed transiting planets, host Teff < 4200 K, 0.5 < Rp < 4 Rearth, P < 100 d
  FGK: same cuts, 4200 <= Teff <= 7000 K   (comparison baseline)

Both use default_flag=1 (one row per planet) and pl_controv_flag=0.
Columns include astrometry for later kinematic-age work:
  ra/dec, sy_plx(+err), sy_pmra/sy_pmdec, sy_dist, st_radv
and rotation: st_rotp(+err).
"""
import urllib.parse
import urllib.request
import pathlib

BASE = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
OUT = pathlib.Path(__file__).resolve().parents[1] / "data"

COLS = ("pl_name,hostname,discoverymethod,"
        "pl_rade,pl_radeerr1,pl_radeerr2,"
        "pl_orbper,pl_insol,"
        "st_teff,st_mass,st_rad,st_met,"
        "st_rotp,st_rotperr1,"
        "st_radv,ra,dec,sy_plx,sy_plxerr1,sy_pmra,sy_pmdec,sy_dist")

CUTS = ("default_flag=1 AND pl_controv_flag=0 "
        "AND discoverymethod='Transit' "
        "AND pl_rade>0.5 AND pl_rade<4 AND pl_orbper<100")

QUERIES = {
    "planet_sample_M.csv":   f"SELECT {COLS} FROM ps WHERE {CUTS} AND st_teff<4200",
    "planet_sample_FGK.csv": f"SELECT {COLS} FROM ps WHERE {CUTS} AND st_teff>=4200 AND st_teff<=7000",
}

def fetch(fname, query):
    url = BASE + "?" + urllib.parse.urlencode({"query": query, "format": "csv"})
    print(f"downloading {fname} ...")
    with urllib.request.urlopen(url, timeout=120) as r, open(OUT / fname, "wb") as f:
        f.write(r.read())
    n = sum(1 for _ in open(OUT / fname)) - 1
    print(f"  -> {n} planets saved")

for fname, q in QUERIES.items():
    fetch(fname, q)
print("done")
