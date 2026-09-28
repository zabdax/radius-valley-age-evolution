#!/usr/bin/env python3
"""fetch_gaia_quality_v2.py — same as fetch_gaia_quality.py but via raw
TAP sync POST (requests, timeout=300s) with bigger chunks. Writes
data/gaia_quality_DR3_v2.csv. Whichever fetch finishes first wins; the
other file is then deduplicated/discarded in analysis.
"""
import io
import pathlib
import time

import numpy as np
import pandas as pd
import requests

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / "data"

SYNC_URL = "https://gea.esac.esa.int/tap-server/tap/sync"
CHUNK = 400


def main():
    ids = []
    for tag in ("M", "FGK"):
        g = pd.read_csv(DATA / f"gaia_hosts_{tag}.csv", usecols=["source_id"])
        ids.append(g)
    ids = pd.concat(ids, ignore_index=True)
    uniq = ids["source_id"].astype(np.int64).unique().tolist()
    print(f"hosts: {len(ids)} unique: {len(uniq)}", flush=True)
    cols = ("source_id, ruwe, parallax, parallax_error, pmra, pmdec, "
            "pmra_error, pmdec_error, astrometric_excess_noise")
    parts = []
    for i in range(0, len(uniq), CHUNK):
        ch = uniq[i:i + CHUNK]
        inlist = ",".join(str(int(x)) for x in ch)
        q = (f"SELECT {cols} FROM gaiadr3.gaia_source "
             f"WHERE source_id IN ({inlist})")
        for attempt in range(5):
            try:
                r = requests.post(SYNC_URL,
                                  data={"REQUEST": "doQuery", "LANG": "ADQL",
                                        "FORMAT": "csv", "QUERY": q},
                                  timeout=300)
                r.raise_for_status()
                tab = pd.read_csv(io.StringIO(r.text))
                break
            except Exception as e:
                print(f"  chunk {i//CHUNK+1} attempt {attempt+1}: {e}",
                      flush=True)
                time.sleep(10 * (attempt + 1))
        else:
            raise RuntimeError(f"chunk {i} failed 5x")
        parts.append(tab)
        print(f"  chunk {i//CHUNK+1}/{(len(uniq)+CHUNK-1)//CHUNK}: "
              f"{len(tab)} rows", flush=True)
    qtab = pd.concat(parts, ignore_index=True).drop_duplicates("source_id")
    qtab.to_csv(DATA / "gaia_quality_DR3_v2.csv", index=False)
    print(f"retrieved {len(qtab)} rows", flush=True)
    print(qtab[["ruwe", "parallax", "parallax_error"]].describe().to_string())


if __name__ == "__main__":
    main()
