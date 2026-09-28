# Reproducibility

All numbers in the manuscript regenerate from committed files. No network needed
except `code/01_download.py` (live TAP query).

## Quick verify (no install beyond Python + stdlib)

```bash
sha256sum -c data/checksums.sha256        # Linux/macOS
# Windows PowerShell:
# Get-FileHash -Algorithm SHA256 -Path data/*.csv
```

## Full pipeline (Python >= 3.10)

```bash
pip install -r requirements.txt   # numpy, pandas, scipy, matplotlib, astropy, gyrointerp
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
python code/audit_highrep.py                         # audits results/highrep/
python code/regen_manuscript_figs.py                 # Fig 2 + Fig 4 from highrep
python figures/make_overview.py                      # hero overview.png/pdf
```

## Authoritative results (already committed)

- `results/highrep/task2_gate_rep1000_replicates_raw.csv` (2,000 reps)
- `results/highrep/task2_calibration_kepler_only_trim0.95_rep1000.json`
- `results/highrep/power_grid_rep500_replicates_raw.csv` (7,500 reps)
- `results/highrep/power_analysis_rep500.json`
- `results/colab_audit/AUDIT_REPORT.md` (independent audit)

## High-replicate run (Colab, ~5–7 h gate)

See `colab/README.md` + `docs/USAGE_COLAB.md`. The as-executed variant
`colab/colab_highrep_gate_fast_v3.py` (SHA-256 `fed4824d…d53a69e59`) drives the
reference `colab_highrep_gate.py` as a library with identical rng order.

## Known limitations (disclosed)

- Bootstrap grid upsamples ~301 M hosts; host-covariate diversity saturates.
- σ ∝ N<sup>−1/2</sup> extrapolation is a premise, not a grid demonstration.
- M sample (T<sub>eff</sub> < 4200 K) spans late-K–late-M, not isolated M3+.
- Isochrone validation incomplete (`code/B1/B2`).
