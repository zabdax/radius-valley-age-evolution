# Reproducibility

All numbers in the current (v2) manuscript regenerate from committed files.
No network access is needed except `code/01_download.py`, which performs the
original live TAP query.

**Read [REPO_STATUS.md](../REPO_STATUS.md) first.** It identifies which results
are authoritative and which older files are historical.

## 1. Verify the committed artifacts

```bash
sha256sum -c data/checksums.sha256                       # 31 frozen input files
cd repro/v2 && sha256sum -c provenance/BUNDLE_SHA256.txt  # 85 v2 bundle files
```

On Windows PowerShell, use `Get-FileHash -Algorithm SHA256` for the same check.

## 2. Reproduce the v2 results (recommended)

The self-contained bundle in `repro/v2/` holds the exact scripts and results
behind the current manuscript:

```text
repro/v2/scripts/          executed analysis + figure scripts (portable paths)
repro/v2/results/science/  authoritative JSON/CSV outputs (*_FINAL = final)
repro/v2/results/audit/    Gaia cuts, demographics, multiplicity budget
repro/v2/figures/          seven publication figures
repro/v2/provenance/       checksums, provenance lock, Zenodo manifests
```

Path resolution is centralised in `repro/v2/scripts/_paths.py`; set the
`RV_REPO_ROOT` environment variable if the repository is relocated. The frozen
inputs stay in top-level `data/` and shared analysis machinery in `code/`.

Key v2 outputs and their drivers:

| Result | Script | Output |
|---|---|---|
| TIME Table crossmatch | `scripts/parse_time_table.py` | `results/science/time_table_crossmatch.json` |
| TIME Table age-ordering check | `scripts/timetable_spearman.py` | `results/science/timetable_spearman.json` |
| Period mixing | `scripts/period_mixing_analysis.py` | `results/science/period_classifier_results.csv` |
| Full power-vs-effect curve | `scripts/complete_power_curve.py` | `results/science/power_vs_effect_summary_FINAL.json` |
| Gaia production cuts | `scripts/gaia_production_cuts.py` | `results/audit/gaia_cut_summary.json` |
| Publication figures | `scripts/regen_fig*.py`, `make_new_figures.py` | `figures/fig1`--`fig7` |

## 3. Rebuild the original pipeline (historical, still runnable)

```bash
pip install -r requirements.txt
python code/01_download.py        # network: NASA Exoplanet Archive TAP
python code/06_gaia_crossmatch.py
python code/06b_filter_gaia.py
python code/07_kinematic_ages.py
python code/03_gyro_ages.py
python code/04_validate.py
python code/12_completeness.py
python code/09_hierarchical_model.py
python code/13_calibration_fgk.py kepler_only 0.95   # original 120-rep run
python code/11_power_analysis.py                     # original grid
python code/14_one_per_host_check.py
python code/audit_highrep.py                         # audit results/highrep/
python code/regen_manuscript_figs.py                 # earlier Fig 2 + Fig 4
python figures/make_overview.py                      # hero overview figure
python figures/make_gate_animation.py                # gate convergence GIF
python figures/make_power_effect_animation.py        # power-vs-effect GIF
```

Run order follows the numeric prefixes in `code/`. All scripts use
repository-relative paths.

## 4. Authoritative raw results

- `results/highrep/task2_gate_rep1000_replicates_raw.csv` (2,000 replicates)
- `results/highrep/task2_calibration_kepler_only_trim0.95_rep1000.json`
- `results/highrep/power_grid_rep500_replicates_raw.csv` (7,500 replicates)
- `results/highrep/power_analysis_rep500.json`
- `results/colab_audit/AUDIT_REPORT.md` (independent audit of the earlier run)

## 5. High-replicate Colab run

See `colab/README.md` and `docs/USAGE_COLAB.md`. The as-executed variant
`colab/colab_highrep_gate_fast_v3.py` (SHA-256 `fed4824d.d53a69e59`) drives the
reference `colab_highrep_gate.py` as a library with an identical RNG order.

## 6. Known limitations (disclosed)

- The bootstrap grid reuses roughly 301 cool hosts; host-covariate diversity
  saturates, so sigma proportional to N^(-1/2) is a premise, not a result.
- The cool sample (Teff < 4200 K) spans late-K to late-M, not isolated M3+.
- Isochrone-age validation was attempted but not completed (`code/B1`, `B2`).
- The Gaia bound-pair crossmatch could not be completed: endpoint failures and
  a 1.42 GB catalog at ~0.07 MB/s. No binary fraction is claimed (see
  `repro/v2/results/audit/P3_bound_pair_crossmatch.md`).
