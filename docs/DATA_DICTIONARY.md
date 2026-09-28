# Data dictionary (frozen inputs)

All files in `data/` are committed snapshots; verify with `data/checksums.sha256`.

## Core gate/power inputs

| File | Rows (approx) | Contents |
|---|---|---|
| `planet_sample_FGK.csv` | 2,880 planets | FGK control sample (transit, default_flag=1, 0.5–4 R⊕, P<100 d) |
| `planet_sample_M.csv` | 419 planets | M-dwarf sample, `st_teff < 4200 K` |
| `hosts_kinematics_FGK.csv` | 2,123 hosts | Gaia DR3 `v_tan`, U, V, W |
| `hosts_kinematics_M.csv` | 301 hosts | same for M |
| `completeness_weights_FGK.csv` | 2,880 | Inverse-detection weights (DR25 CDPP, trimmed) |
| `completeness_weights_M.csv` | 419 | M weights |
| `completeness_weights_M_sharp.csv` | 419 | Sharpened M sensitivity variant |

## Supporting pipeline tables

- `dr25_cdpp.csv`, `dr25_cdpp_matched.csv` — Kepler DR25 6-h CDPP + host match
- `gaia_hosts_FGK.csv`, `gaia_hosts_M.csv` — raw DR3 crossmatch; `gaia_rejects_*.txt` (empty: rejections logged elsewhere)
- `gaia_quality_DR3.csv`, `gaia_quality_DR3_v2.csv` — RUWE / parallax-S/N quality
- `mags_FGK.csv`, `mags_M.csv`, `magtest.csv` — host magnitudes / probe
- `dr25_stellar_ages.csv`, `host_ages_FGK.csv` (79 usable), `host_ages_M.csv` (12 usable) — gyro ages
- `mcquillan_cool.csv`, `koi_names.csv`, `koi_ages.csv`, `phaseB1_hosts.txt` — rotation/KOI crossmatches
- `published_ages_FGK.csv`, `berger_h.tsv`, `berger_probe.tsv`, `ks_test.csv`, `sh_test.csv` — isochrone-attempt I/O (incomplete by design)

## Authoritative outputs (`results/highrep/`)

| File | Contents |
|---|---|
| `task2_gate_rep1000_replicates_raw.csv` | 2,000 replicate fits (1,000 null + 1,000 signal) |
| `task2_calibration_kepler_only_trim0.95_rep1000.json` | Gate summary: beta_real, envelopes, percentiles, verdict |
| `power_grid_rep500_replicates_raw.csv` | 7,500 cell replicates (15 cells × 500) |
| `power_analysis_rep500.json` | Per-cell detection/coverage table |
| `manifest.json`, `accelerator_manifest.json` | Run provenance |

Derived data behind the RNAAS figure lives in `manuscript/`-adjacent CSVs in the
full project archive; the open-source figure CSVs regenerate via
`figures/make_overview.py` and `code/regen_manuscript_figs.py`.
