"""
06_gaia_crossmatch.py — Match every M-dwarf host to Gaia DR3.

Purpose: fetch radial_velocity, BP-RP photometry, and confirm astrometry
for each host -> inputs for kinematic ages (and any future gyro work).

Method: sequential cone searches (1.5 arcsec) against the local copy of
positions pulled from the NASA Exoplanet Archive. Results saved
incrementally; re-runs skip already-matched hosts.
"""
import pathlib
import sys
import time
import urllib.parse
import urllib.request
import pandas as pd

D = pathlib.Path(__file__).resolve().parents[1] / "data"
SAMPLE = sys.argv[1] if len(sys.argv) > 1 else "planet_sample_M.csv"
TAG = "M" if "_M" in SAMPLE else "FGK"
OUT = D / f"gaia_hosts_{TAG}.csv"

QUERY = """
SELECT TOP 10 source_id, ra, dec, parallax, pmra, pmdec, radial_velocity,
       phot_bp_mean_mag, phot_rp_mean_mag, phot_g_mean_mag
FROM gaiadr3.gaia_source
WHERE 1=CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra:.6f}, {dec:.6f}, {r}))
"""
# NOTE: no `ORDER BY COORDDISTANCE` here -- ESA TAP rejects it with HTTP 400
# (costs the whole run). We pick the nearest candidate locally instead.


def gaia_cone(ra, dec, radius_arcsec=1.5):
    q = QUERY.format(ra=ra, dec=dec, r=radius_arcsec / 3600.0)
    url = ("https://gea.esac.esa.int/tap-server/tap/sync?"
           + urllib.parse.urlencode({"REQUEST": "doQuery", "LANG": "ADQL",
                                     "FORMAT": "csv", "QUERY": q}))
    with urllib.request.urlopen(url, timeout=90) as r:
        text = r.read().decode()
    lines = [l for l in text.splitlines() if l.strip()]
    if len(lines) < 2:
        return None
    # nearest row by plain angular separation (arcsec-scale, flat-sky ok)
    import math
    best, best_d = None, 1e9
    for line in lines[1:]:
        parts = line.split(",")
        cra, cdec = float(parts[1]), float(parts[2])
        d = math.hypot((cra - ra) * math.cos(math.radians(dec)),
                       cdec - dec) * 3600.0
        if d < best_d:
            best_d, best = d, line
    return best


def main():
    m = pd.read_csv(D / SAMPLE)
    hosts = m.sort_values("pl_rade").groupby("hostname", as_index=False).first()[
        ["hostname", "st_teff", "ra", "dec"]]

    names = ["hostname", "source_id", "ra", "dec", "parallax", "pmra",
             "pmdec", "radial_velocity", "phot_bp_mean_mag",
             "phot_rp_mean_mag", "phot_g_mean_mag"]
    reject_log = D / f"gaia_rejects_{TAG}.txt"
    done = set()
    if OUT.exists():
        # tolerate legacy integer-header files written by the first version
        first = open(OUT).readline().strip()
        skip = 1 if first.split(",")[0] == "0" else 0
        done = set(pd.read_csv(OUT, names=names, skiprows=skip).hostname)
        print(f"resuming: {len(done)} hosts already matched")

    def _f(x):
        try:
            return float(x)
        except (TypeError, ValueError):
            return float("nan")

    n_new = 0
    backoff = [2, 8, 30, 90]          # seconds; ESA TAP throttles fast cadences
    for i, h in hosts.iterrows():
        if h.hostname in done:
            continue
        row = None
        for attempt, wait in enumerate([0] + backoff):
            time.sleep(wait)
            try:
                row = gaia_cone(float(h.ra), float(h.dec))
                break
            except Exception as e:
                if attempt == len(backoff):
                    print(f"SKIP {h.hostname} after {len(backoff)} retries: {e}",
                          flush=True)
        if row is not None:
            parts = row.split(",")
            # row layout: source_id, ra, dec, parallax, pmra, pmdec, rv, ...
            plx, pmra_, pmdec_ = _f(parts[3]), _f(parts[4]), _f(parts[5])
            if not (plx > 0 and pmra_ == pmra_ and pmdec_ == pmdec_):
                # Gaia source without astrometry: spurious detection or
                # companion -- not usable for kinematics. Log and move on.
                msg = f"{h.hostname}\t{parts[0]}\tplx={parts[3]} pm={parts[4]},{parts[5]}"
                with open(reject_log, "a") as fh:
                    fh.write(msg + "\n")
                print(f"REJECT {msg}", flush=True)
                done.add(h.hostname)
                continue
            pd.DataFrame([[h.hostname] + parts], columns=names).to_csv(
                OUT, mode="a", header=not OUT.exists(), index=False)
            done.add(h.hostname)
            n_new += 1
        if i % 20 == 0:
            print(f"[{i+1}/{len(hosts)}] matched so far: {len(done)}", flush=True)
        time.sleep(1.2)               # stay under the rate limit

    print(f"done. new matches: {n_new}, total: {len(done)}/{len(hosts)}")


if __name__ == "__main__":
    main()
