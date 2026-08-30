"""
B1_get_berger_ages.py — Pull published ages for our planetary hosts explicitly.
Uses the main `ps` table which compiles the most robust published ages,
ensuring we get the Berger et al. 2020 ages for the Kepler field specifically.
"""
import pathlib
import pandas as pd
import urllib.parse
import urllib.request

D = pathlib.Path(__file__).resolve().parents[1] / "data"

def main():
    pl_f = pd.read_csv(D / "planet_sample_FGK.csv")
    hosts = pl_f[pl_f.hostname.str.startswith("Kepler")].hostname.unique()
    print(f"Target FGK Kepler hosts: {len(hosts)}")

    url = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
    q = (f"SELECT hostname, st_age, st_ageerr1, st_ageerr2 "
         f"FROM ps WHERE default_flag=1")
    data = urllib.parse.urlencode(
        {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": q}
    ).encode()

    print("Fetching ages from ps table via TAP...")
    with urllib.request.urlopen(urllib.request.Request(url, data=data), timeout=180) as r:
        raw = pd.read_csv(pd.io.common.StringIO(r.read().decode()))

    # Group by hostname and take the first non-null age
    ag = raw.dropna(subset=['st_age']).groupby('hostname', as_index=False).first()
    print(f"Ages globally fetched: {len(ag)}")

    matched = pd.DataFrame({'hostname': hosts}).merge(ag, on='hostname', how='left')
    n_match = matched.st_age.notna().sum()
    print(f"Coverage on our FGK Kepler hosts: {n_match} / {len(hosts)} "
          f"({100 * n_match / len(hosts):.1f}%)")

    matched.to_csv(D / "published_ages_FGK.csv", index=False)
    print("Saved.")

if __name__ == "__main__":
    main()
