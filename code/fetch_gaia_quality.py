#!/usr/bin/env python3
"""fetch_gaia_quality.py — pull RUWE + parallax errors for gate hosts.

Reads source_id lists from data/gaia_hosts_{M,FGK}.csv, queries Gaia DR3
TAP (gea.esac.esa.int) in chunks, writes data/gaia_quality_DR3.csv with:
  source_id, ruwe, parallax, parallax_error, pmra, pmdec,
  pmra_error, pmdec_error, astrometric_excess_noise
Frozen archive tables are NOT modified.
"""
import pathlib
import time

import numpy as np
import pandas as pd
import pyvo

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / "data"

TAP_URL = "https://gea.esac.esa.int/tap-server/tap"
CHUNK = 200


def main():
    ids = []
    for tag in ("M", "FGK"):
        g = pd.read_csv(DATA / f"gaia_hosts_{tag}.csv", usecols=["source_id"])
        g["tag"] = tag
        ids.append(g)
    ids = pd.concat(ids, ignore_index=True)
    ids["source_id"] = ids["source_id"].astype(np.int64)
    print(f"hosts: {len(ids)} unique source_ids: {ids.source_id.nunique()}")
    svc = pyvo.dal.TAPService(TAP_URL)
    cols = ("source_id, ruwe, parallax, parallax_error, pmra, pmdec, "
            "pmra_error, pmdec_error, astrometric_excess_noise")
    parts = []
    uniq = ids["source_id"].unique().tolist()
    for i in range(0, len(uniq), CHUNK):
        ch = uniq[i:i + CHUNK]
        inlist = ",".join(str(int(x)) for x in ch)
        q = (f"SELECT {cols} FROM gaiadr3.gaia_source "
             f"WHERE source_id IN ({inlist})")
        for attempt in range(4):
            try:
                job = svc.submit_job(q)
                job.run()
                job.wait(phases=["COMPLETED", "ERROR", "ABORTED"], timeout=120)
                if job.phase != "COMPLETED":
                    raise RuntimeError(f"TAP phase {job.phase}: "
                                       f"{job.errors}")
                tab = job.fetch_result().to_table().to_pandas()
                break
            except Exception as e:
                print(f"  chunk {i//CHUNK+1} attempt {attempt+1} failed: {e}")
                time.sleep(5 * (attempt + 1))
        else:
            raise RuntimeError(f"chunk starting {i} failed 4x")
        parts.append(tab)
        print(f"  chunk {i//CHUNK+1}/{(len(uniq)+CHUNK-1)//CHUNK}: "
              f"{len(tab)} rows", flush=True)
    qtab = pd.concat(parts, ignore_index=True).drop_duplicates("source_id")
    print(f"retrieved {len(qtab)} rows; missing "
          f"{len(set(uniq)) - len(set(qtab.source_id))}")
    qtab.to_csv(DATA / "gaia_quality_DR3.csv", index=False)
    print("wrote data/gaia_quality_DR3.csv")
    print(qtab[["ruwe", "parallax", "parallax_error"]].describe().to_string())


if __name__ == "__main__":
    main()
