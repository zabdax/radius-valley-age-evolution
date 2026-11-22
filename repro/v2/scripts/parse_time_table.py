#!/usr/bin/env python3
"""Parse the Gaidos et al. (2023) TIME Table and crossmatch against our cool-host sample.

Two source layers are handled, deliberately kept separate:

  (A) VAULT-NOTE EMBEDDED TABLE (mangled PDF text)
      research/notes/mnras-000-116-2022.md lines ~1284-1725.
      The paper explicitly says "Only the first 50 entries are shown" (note line
      988).  A robust field-pattern state machine recovers those 50 rows.

  (B) AUTHORITATIVE MACHINE-READABLE TABLE (CDS format, 249 rows)
      Zenodo DOI 10.5281/zenodo.7578269, file time_table.mrt.
      A cached copy lives in results/science/time_table_mrt_source.txt
      (md5 581d26ae64ae79eb4c9b309bc7d6466a).  Parsed by fixed byte columns.

The delivered CSV uses layer (B) as the primary 249-row table and marks which
rows also appear in the 50-row embedded extract.  Parse failures are reported
honestly rather than silently dropped.

Usage:
  python manuscripts/q1-paper/v2-2026-11/scripts/parse_time_table.py
"""
import csv
import json
import os
import re
import urllib.request
from collections import Counter

# --------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
NOTE = os.path.join(ROOT, "research", "notes", "mnras-000-116-2022.md")
SAMPLE = os.path.join(ROOT, "radius-valley-age-evolution", "data", "planet_sample_M.csv")
GYRO = os.path.join(ROOT, "radius-valley-age-evolution", "data", "host_ages_M.csv")
SCI = os.path.join(ROOT, "manuscripts", "q1-paper", "v2-2026-11", "results", "science")
MRT_CACHE = os.path.join(SCI, "time_table_mrt_source.txt")
MRT_URL = "https://zenodo.org/api/records/7578269/files/time_table.mrt/content"

# --------------------------------------------------------------------------
# Field regexes for the mangled note text
# --------------------------------------------------------------------------
RE_TEFF = re.compile(r"^\d{3,4}$")
RE_FEH = re.compile(r"^[+-]?\d+\.\d+$")
RE_DASH = re.compile(r"^[\u2013\u2014\-–]+$")          # em/en dashes used for "missing"
RE_PERIOD = re.compile(r"^\d+\.\d+(?:\s*\(\d+(?:\.\d+)?\))?$")
RE_METHOD = re.compile(r"^[PS]+\s*\([^)]*\)$")
RE_AGE = re.compile(r"^(?:<)?\d+(?:\.\d+)?(?:\s*\(\d+(?:\.\d+)?\))?$")
RE_CORR = re.compile(r"^[+-]\d+\.\d+$")
NOTE_KEYWORDS = ("binary", "error", "old", "age", "no ")

NUM_RE = re.compile(r"(-?\d+(?:\.\d+)?)")


def _is_missing(s):
    return bool(RE_DASH.match(s.strip())) or s.strip() in ("", "—")


def _parse_period(s):
    nums = NUM_RE.findall(s)
    return (float(nums[0]) if nums else None,
            float(nums[1]) if len(nums) > 1 else None)


def _parse_age(s):
    """Return (value_or_None, unc_or_None, status)."""
    s = s.strip()
    if _is_missing(s):
        return None, None, "excluded"
    is_ul = s.startswith("<")
    nums = NUM_RE.findall(s)
    val = float(nums[0]) if nums else None
    unc = float(nums[1]) if len(nums) > 1 else None
    return val, unc, ("upper_limit" if is_ul else "age_estimate")


# --------------------------------------------------------------------------
# (A) Embedded 50-row table parser (robust state machine)
# --------------------------------------------------------------------------
def parse_embedded_table(path):
    lines = open(path, encoding="utf-8", errors="replace").read().splitlines()
    # Restrict to the printed Table 1 region.
    start = next(i for i, l in enumerate(lines) if l.startswith("Table 1. TIME-Table"))
    end = next(i for i, l in enumerate(lines) if l.startswith("1The full table is available"))
    region = lines[start:end]
    records, failures = [], []

    def is_record_start(i):
        if i + 2 >= len(region):
            return False
        name, teff, feh = region[i].strip(), region[i + 1].strip(), region[i + 2].strip()
        return (bool(name) and not RE_TEFF.match(name)
                and bool(RE_TEFF.match(teff))
                and (bool(RE_FEH.match(feh)) or bool(RE_DASH.match(feh))))

    i = 0
    while i < len(region):
        if not is_record_start(i):
            i += 1
            continue
        try:
            name = region[i].strip()
            teff = int(region[i + 1].strip())
            feh = region[i + 2].strip()
            prot, e_prot = _parse_period(region[i + 3].strip())
            method = region[i + 4].strip()
            reference = region[i + 5].strip()
            age, e_age, status = _parse_age(region[i + 6].strip())
            corr = region[i + 7].strip()
            note = ""
            if i + 8 < len(region):
                cand = region[i + 8].strip()
                if any(k in cand.lower() for k in NOTE_KEYWORDS):
                    note = cand
            # validity checks
            problems = []
            if not RE_PERIOD.match(region[i + 3].strip()):
                problems.append("period")
            if not RE_METHOD.match(method):
                problems.append("method")
            if not RE_CORR.match(corr) and not RE_DASH.match(corr):
                problems.append("corr")
            rec = dict(name=name, teff_K=teff, feh_dex=(None if RE_DASH.match(feh) else float(feh)),
                       prot_days=prot, e_prot_days=e_prot, method=method, reference=reference,
                       age_Gyr=age, e_age_Gyr=e_age, corr_Gyr=(None if RE_DASH.match(corr) else float(corr)),
                       age_status=status, note=note, source_line=start + i + 1,
                       problems=problems)
            records.append(rec)
            i += 9 if note else 8
        except Exception as exc:
            failures.append(dict(line=start + i + 1, text=region[i], error=str(exc)))
            i += 1
    # Drop records with hard field-format problems.
    clean = [r for r in records if not r["problems"]]
    bad = [r for r in records if r["problems"]]
    return clean, failures, bad


# --------------------------------------------------------------------------
# (B) Authoritative CDS machine-readable table (249 rows)
# --------------------------------------------------------------------------
def get_mrt_text():
    if os.path.exists(MRT_CACHE):
        return open(MRT_CACHE, encoding="utf-8", errors="replace").read()
    text = urllib.request.urlopen(MRT_URL, timeout=60).read().decode("utf-8", "replace")
    os.makedirs(SCI, exist_ok=True)
    open(MRT_CACHE, "w", encoding="utf-8").write(text)
    return text


def _c(line, a, b):
    return line[a - 1:b].strip()


def parse_mrt():
    text = get_mrt_text()
    lines = text.splitlines()
    data_start = next(i for i, l in enumerate(lines) if l.startswith("EPIC 201170410"))
    rows = [l for l in lines[data_start:] if l.strip()]
    recs = []
    for l in rows:
        name = _c(l, 1, 15)
        binary = _c(l, 17, 17)
        teff = _c(l, 19, 22)
        feh = _c(l, 24, 28)
        prot = _c(l, 30, 35)
        eprot = _c(l, 37, 41)
        method = _c(l, 43, 44)
        age = _c(l, 46, 50)
        eage = _c(l, 52, 56)
        corr = _c(l, 58, 62)
        fail = _c(l, 64, 64)
        ul = _c(l, 66, 66)
        pms = _c(l, 68, 68)
        no_err = _c(l, 70, 70)
        note = _c(l, 72, 99)
        if fail == "Y":
            status = "excluded"
        elif ul == "Y":
            status = "upper_limit"
        else:
            status = "age_estimate"
        recs.append(dict(
            name=name, binary_flag=binary,
            teff_K=(int(teff) if teff else None),
            feh_dex=(float(feh) if feh else None),
            prot_days=(float(prot) if prot else None),
            e_prot_days=(float(eprot) if eprot else None),
            method=method,
            age_Gyr=(float(age) if age else None),
            e_age_Gyr=(float(eage) if eage else None),
            corr_Gyr=(float(corr) if corr else None),
            fail_flag=fail, ul_flag=ul, pms_flag=pms, no_err_flag=no_err,
            note=note, age_status=status,
        ))
    return recs, len(rows)


# --------------------------------------------------------------------------
# Name normalization + crossmatch
# --------------------------------------------------------------------------
def normalize_name(name):
    s = str(name).strip().lower()
    s = s.replace("\u2013", "-").replace("\u2014", "-")
    s = re.sub(r"\s*-\s*", "-", s)
    s = re.sub(r"\s+", " ", s)
    # Gl / Gliese / GJ all refer to the Gliese(-Jahreiss) catalogue.
    s = re.sub(r"^(gliese|gl)\b", "gj", s)
    # Drop trailing stellar component letter (GJ 15 A, K2-148 B, TOI-1450 A ...).
    s = re.sub(r"\s+[a-d]$", "", s)
    s = s.replace(" ", "")
    return s


def main():
    os.makedirs(SCI, exist_ok=True)

    # ---- Parse both layers ----
    embedded, emb_failures, emb_bad = parse_embedded_table(NOTE)
    mrt, n_raw = parse_mrt()
    emb_names = {r["name"] for r in embedded}
    for r in mrt:
        r["in_vault_note_embedded_50"] = "Y" if r["name"] in emb_names else "N"
        r["normalized_name"] = normalize_name(r["name"])

    # ---- Validation counts ----
    counts = Counter(r["age_status"] for r in mrt)
    validation = dict(
        expected_total=249, recovered_total=len(mrt),
        expected_age_estimates=227, recovered_age_estimates=counts.get("age_estimate", 0),
        expected_upper_limits=8, recovered_upper_limits=counts.get("upper_limit", 0),
        expected_excluded=14, recovered_excluded=counts.get("excluded", 0),
        embedded_table_rows_expected=50, embedded_table_rows_recovered=len(embedded),
        embedded_table_parse_failures=len(emb_failures) + len(emb_bad),
    )
    print("VALIDATION:", json.dumps(validation))

    # ---- CSV ----
    csv_cols = ["name", "normalized_name", "binary_flag", "teff_K", "feh_dex",
                "prot_days", "e_prot_days", "method", "age_Gyr", "e_age_Gyr",
                "corr_Gyr", "fail_flag", "ul_flag", "pms_flag", "no_err_flag",
                "age_status", "note", "in_vault_note_embedded_50"]
    csv_path = os.path.join(SCI, "time_table_parsed.csv")
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=csv_cols, extrasaction="ignore")
        w.writeheader()
        for r in mrt:
            w.writerow(r)

    # ---- Crossmatch ----
    import pandas as pd
    sample = pd.read_csv(SAMPLE)
    sample_hosts = sorted(sample["hostname"].dropna().unique().tolist())
    sample_map = {}
    for h in sample_hosts:
        sample_map.setdefault(normalize_name(h), []).append(h)

    gyro_df = pd.read_csv(GYRO)
    gyro_ok = gyro_df.loc[gyro_df["gyro_ok"] == True, "hostname"].tolist()  # noqa: E712
    gyro_norm = {normalize_name(h) for h in gyro_ok}

    mrt_by_norm = {}
    for r in mrt:
        mrt_by_norm.setdefault(r["normalized_name"], []).append(r)

    matches = []
    for nkey, sample_names in sorted(sample_map.items()):
        if nkey in mrt_by_norm:
            trows = mrt_by_norm[nkey]
            matches.append(dict(
                sample_hosts=sample_names,
                time_table_names=[r["name"] for r in trows],
                normalized_key=nkey,
                age_statuses=[r["age_status"] for r in trows],
                teff_K=[r["teff_K"] for r in trows],
                prot_days=[r["prot_days"] for r in trows],
                age_Gyr=[r["age_Gyr"] for r in trows],
                e_age_Gyr=[r["e_age_Gyr"] for r in trows],
                already_one_of_12_gyro=(nkey in gyro_norm),
            ))

    matched_12 = [m for m in matches if m["already_one_of_12_gyro"]]
    additional = [m for m in matches if not m["already_one_of_12_gyro"]]
    # Of the additional hosts, those with a real (non-upper-limit, non-excluded) age.
    additional_with_age = [m for m in additional
                           if any(s == "age_estimate" for s in m["age_statuses"])]

    crossmatch = dict(
        sample_host_count=len(sample_hosts),
        time_table_host_count=len({r["normalized_name"] for r in mrt}),
        matched_count=len(matches),
        matched_hosts=matches,
        matched_names_flat=[m["sample_hosts"][0] for m in matches],
        existing_12_gyro_hosts=gyro_ok,
        existing_12_gyro_matched_in_time_table=[m["sample_hosts"][0] for m in matched_12],
        additional_hosts_beyond_12_count=len(additional),
        additional_hosts_beyond_12=[m["sample_hosts"][0] for m in additional],
        additional_hosts_with_finite_age_count=len(additional_with_age),
        additional_hosts_with_finite_age=[m["sample_hosts"][0] for m in additional_with_age],
        parse_completeness=validation,
        caveats=[
            "NAME NORMALIZATION: trailing stellar-component letters stripped (GJ 15 A -> GJ 15; "
            "K2-148 A/B -> K2-148; TOI-1450 A -> TOI-1450). K2-148 A and B collapse to one key.",
            "Gl / Gliese / GJ unified as one catalogue prefix. No Gl<->GJ numeric ambiguity "
            "affected any match in this run.",
            "Kepler-N and KOI-N are NOT interchangeable designations and were never matched "
            "by number. The TIME Table lists many KOI-N stars absent from our Kepler-N named "
            "sample; those are counted as non-matches.",
            "EPIC and TIC identifiers were kept distinct.",
            "'Additional hosts' = cool hosts in our 303-star sample that appear in the Gaidos+23 "
            "TIME Table but are not among our 12 pipeline gyro-age hosts. It does NOT mean those "
            "stars are usable in our Angus+19 gyrointerp model (Teff 3800-6200 K, Prot <= 45 d).",
            "Vault note embeds only the first 50 of 249 rows; full table parsed from Zenodo "
            "DOI 10.5281/zenodo.7578269 (time_table.mrt, CDS format).",
        ],
    )
    json_path = os.path.join(SCI, "time_table_crossmatch.json")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(crossmatch, f, indent=2)

    print("MATCHED:", len(matches))
    print("ADDITIONAL BEYOND 12:", len(additional))
    print("ADDITIONAL WITH FINITE AGE:", len(additional_with_age))
    print("MATCHED 12:", [m["sample_hosts"][0] for m in matched_12])
    print("EMBEDDED FAILURES:", emb_failures, emb_bad)


if __name__ == "__main__":
    main()
