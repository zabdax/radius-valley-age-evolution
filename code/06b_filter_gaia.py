"""
06b_filter_gaia.py — Offline quality gate for existing Gaia crossmatch files.

Applies the same acceptance rules as 06's live gate (finite parallax > 0,
finite proper motions) to already-downloaded tables, WITHOUT re-querying:
  data/gaia_hosts_M.csv   (legacy integer-header version from _attic has all
                           303 rows; the post-gate regeneration wrongly
                           dropped 187 hosts due to a column-index bug)
  data/gaia_hosts_FGK.csv (written before any gate existed)

Writes clean gaia_hosts_{TAG}.csv with proper headers + rejected rows to
gaia_rejects_{TAG}.txt.
"""
import pathlib
import pandas as pd

D = pathlib.Path(__file__).resolve().parents[1] / "data"
NAMES = ["hostname", "source_id", "ra", "dec", "parallax", "pmra",
         "pmdec", "radial_velocity", "phot_bp_mean_mag",
         "phot_rp_mean_mag", "phot_g_mean_mag"]


def load_any(tag):
    p = D / f"gaia_hosts_{tag}.csv"
    first = open(p).readline().strip()
    skip = 1 if first.split(",")[0] in ("0", "hostname") else 0
    df = pd.read_csv(p, names=NAMES, skiprows=skip)
    for c in NAMES[1:]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def main(tag):
    df = load_any(tag)
    n0 = len(df)
    ok = (df.parallax.notna() & (df.parallax > 0)
          & df.pmra.notna() & df.pmdec.notna())
    bad = df[~ok]
    good = df[ok]
    out = D / f"gaia_hosts_{tag}.csv"
    good.to_csv(out, index=False)          # proper header, clean rows
    rej = D / f"gaia_rejects_{tag}.txt"
    with open(rej, "w") as fh:
        for _, r in bad.iterrows():
            fh.write(f"{r.hostname}\t{r.source_id}\tplx={r.parallax} "
                     f"pm={r.pmra},{r.pmdec}\n")
    print(f"[{tag}] {n0} -> {len(good)} accepted, {len(bad)} rejected "
          f"(logged to {rej.name})")
    print(f"  rejected: {', '.join(bad.hostname.head(10))}"
          + (" ..." if len(bad) > 10 else ""))


if __name__ == "__main__":
    main("M")
    main("FGK")
