# Data

All files in this folder are frozen snapshots of pipeline inputs and outputs.
They are committed so that every number in the manuscript is reproducible
without re-querying remote services (archive contents can drift over time).
`checksums.sha256` contains the SHA-256 hash of every file listed below —
verify with `sha256sum -c checksums.sha256` (Linux/macOS) or
`Get-FileHash -Algorithm SHA256` (Windows).

## Core analysis inputs (used by the calibrated gate and power grid)

| File | Contents | Produced by |
|---|---|---|
| `planet_sample_FGK.csv` | FGK control planet sample (2,880 planets, 0.5-4.0 R_e, P<100 d) | `code/01_download.py` |
| `planet_sample_M.csv` | M dwarf planet sample (419 planets) | `code/01_download.py` |
| `hosts_kinematics_FGK.csv` | Gaia DR3 v_tan / U,V,W for 2,123 FGK hosts | `code/06_gaia_crossmatch.py` + `code/07_kinematic_ages.py` |
| `hosts_kinematics_M.csv` | Gaia DR3 kinematics for 301 M dwarf hosts | same |
| `completeness_weights_FGK.csv` | Inverse-detection-efficiency weights (Kepler DR25 CDPP, trimmed) | `code/12_completeness.py` |

## Supporting pipeline data

- `dr25_cdpp.csv`, `dr25_cdpp_matched.csv` — Kepler DR25 6-h CDPP stellar
  table and the host match (used by `code/12_completeness.py`)
- `gaia_hosts_FGK.csv`, `gaia_hosts_M.csv` — raw Gaia DR3 crossmatch tables
  (`code/06_gaia_crossmatch.py`); `gaia_rejects_*.txt` record quality-gate
  rejections (both empty: the offline gate rejected hosts recorded elsewhere)
- `mags_FGK.csv`, `mags_M.csv` — host magnitudes for the IDEM completeness
  proxy; `magtest.csv` — probe output
- `dr25_stellar_ages.csv`, `host_ages_FGK.csv`, `host_ages_M.csv` —
  gyrochronological ages (`code/03_gyro_ages.py`, Angus et al. 2019 via
  `gyrointerp`)
- `mcquillan_cool.csv`, `koi_names.csv`, `koi_ages.csv`, `phaseB1_hosts.txt` —
  rotation-period and KOI crossmatch tables
  (`code/05_literature_rotations.py`, `code/04_validate.py`)
- `published_ages_FGK.csv`, `berger_h.tsv`, `berger_probe.tsv`,
  `ks_test.csv`, `sh_test.csv` — inputs/outputs of the isochrone-based
  validation attempt (`code/B1_get_berger_ages.py`,
  `code/B2_isochrone_gate.py`) that was **not completed** and is reported as
  a limitation in the manuscript

## Not redistributed

- `data/berger/` (local only) — downloaded supplementary material of
  Berger et al. (2020, AJ 160, 108; GKSPC planet/host catalog) used by the
  isochrone attempt. Third-party material is not relicensed here; obtain it
  from the publisher's supplementary resources if you need to re-run
  `code/B2_isochrone_gate.py`.
- Superseded caches (`_attic/` in the pre-release working folder) were
  removed entirely.