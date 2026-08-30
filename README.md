# The Radius Valley's Age Evolution Around M Dwarfs: An Inconclusive Verdict and a Quantified Path Forward

Code, data, and results for a calibrated statistical gate that asks whether
current archival samples can discriminate a literature-sized age-evolution
effect in radius-valley demographics around M dwarfs from a null effect.

**Central result.** Applied to a completeness-weighted FGK Kepler-only
control (2,135 planets, ESS = 1,645) with 1,000 Monte Carlo replicates per
hypothesis, the gate returns **FAIL (inconclusive)**: the real age-proxy
coefficient (beta_age/IQR = -0.0296 +/- 0.1061) lies at the 45.6th percentile
of the null-injected distribution and the 81.2nd percentile of the
literature-sized signal-injected distribution — both hypotheses remain
admissible. A wrong-sign binned age signal is identified as a **period-mixing
artifact** of this dataset. For the M-dwarf question, a 15-cell power grid
(500 replicates per cell) shows that at 8x the current inventory
(~3,337 planets) a literature-sized effect is **detected in only 23.2%**
(95% empirical range 19.6-27.2%) of realizations — no tested scale approaches
conventional 80% power, so 8x is a floor, not a sufficient target. A simple
Gaussian independent-host extrapolation places 80% power at order
1e4-2e4 planets (point estimate ~1.7e4), a model-dependent guide rather than
an empirically demonstrated requirement.

The manuscript is in [`manuscript/`](manuscript/) (AASTeX v7.01,
submission-ready).

## Repository layout

```
├── code/                     numbered pipeline scripts (run order = number order)
│   ├── 01_download.py        NASA Exoplanet Archive (TAP) sample construction
│   ├── 02_explore.py         exploratory checks
│   ├── 03_gyro_ages.py       gyrochronological ages (gyrointerp / Angus+2019)
│   ├── 04_validate.py        sample validation
│   ├── 05_literature_rotations.py   McQuillan+2014 rotation crossmatch
│   ├── 06_gaia_crossmatch.py Gaia DR3 crossmatch (+ 06b quality filter)
│   ├── 07_kinematic_ages.py  v_tan and U,V,W kinematic age proxies
│   ├── 08_headline_fgk_vtan.py      FGK headline fit
│   ├── 09_hierarchical_model.py     mixture model + host-robust sandwich SE
│   ├── 11_power_analysis.py         15-cell Monte Carlo power grid
│   ├── 12_completeness.py           DR25-CDPP completeness weighting (+ 12b sharpened M variant)
│   ├── 13_calibration_fgk.py        calibrated gate (null/signal injection)
│   ├── 14_one_per_host_check.py     within-host clustering robustness check
│   ├── B1_get_berger_ages.py, B2_isochrone_gate.py   isochrone validation attempt (incomplete)
│   ├── generate_figures.py          original 120/150-replicate figure pipeline
│   ├── regen_manuscript_figs.py     authoritative high-replicate Fig 2 / Fig 4
│   ├── audit_highrep.py             re-runnable audit of results/highrep/
│   └── audits/                      the three adversarial audit rounds (see audits/README.md)
├── colab/                    high-replicate Colab experiment (1,000/arm gate, 500/cell grid)
├── data/                     frozen input/output CSVs + provenance and SHA-256 checksums
├── figures/                  manuscript figures (Fig 2 / Fig 4 = high-replicate versions)
├── results/                  pipeline outputs; results/highrep/ = authoritative raw replicate results
├── results/colab_audit/      independent audit of the high-replicate run (AUDIT_REPORT.md)
├── manuscript/               AASTeX source, figures, revision log, compile instructions
└── docs/                     Colab usage guide, project notes, legacy draft
```

## Requirements

Python >= 3.10 with:

```
numpy, pandas, scipy, matplotlib, astropy, gyrointerp
```

(`pip install -r requirements.txt`; `gyrointerp` implements the Angus et al.
2019 gyrochronology model. Network access is needed only for
`code/01_download.py`, which queries the NASA Exoplanet Archive TAP service.)

## Reproducing the analysis

All scripts are run from the repository root and write to `results/` via
paths relative to the repo; no absolute paths are used.

1. **Sample construction and validation** (requires network):
   `python code/01_download.py` -> `code/06_gaia_crossmatch.py` ->
   `code/06b_filter_gaia.py` -> `code/07_kinematic_ages.py` ->
   `code/03_gyro_ages.py` -> `code/04_validate.py`.
   The frozen outputs of every step are already committed under `data/`, so
   steps 2-7 work offline against the exact snapshots used in the paper.
2. **Completeness weighting**: `python code/12_completeness.py`
   (FGK Kepler-only DR25-CDPP weights; `12b_sharpen_weights.py` is the
   sharpened M-dwarf sensitivity variant).
3. **Hierarchical mixture model fits**: `python code/09_hierarchical_model.py`
   (add `SELFTEST` for the synthetic recovery self-test).
4. **Calibrated gate (original 120-replicate run)**:
   `python code/13_calibration_fgk.py kepler_only 0.95`.
5. **Power grid (original 120/150-replicate run)**:
   `python code/11_power_analysis.py`.
6. **High-replicate experiment (authoritative)**: see
   [`colab/README.md`](colab/README.md) — 1,000 replicates per gate arm and
   500 per grid cell on Google Colab. The archived outputs of this run are in
   [`results/highrep/`](results/highrep/) and are the numbers used in the
   manuscript.
7. **Figures**: `python code/regen_manuscript_figs.py` regenerates the
   authoritative Figure 2 and Figure 4 from `results/highrep/`;
   `python code/generate_figures.py` regenerates the original-run versions of
   all four figures (archived for the record).
8. **Robustness and audit**: `python code/14_one_per_host_check.py`
   (within-host clustering check); `python code/audit_highrep.py`
   (integrity + statistics audit of `results/highrep/`).

## Key authoritative results (high-replicate run)

| Quantity | Value |
|---|---|
| Real beta_age/IQR (FGK control) | -0.0296 +/- 0.1061 |
| Null-injected median / 95% envelope (n=1000) | -0.0186 / [-0.2476, +0.2246] |
| Real percentile in null (binomial 95%) | 45.6% [42.5, 48.7] |
| P(null >= real) | 0.544 |
| Signal-injected median / 95% envelope (n=1000) | -0.1365 / [-0.3898, +0.0836] |
| Real percentile in signal | 81.2% [78.6, 83.6] |
| Gate verdict | FAIL (inconclusive) |
| 8x SWEET detection (n=500) | 23.2% [19.6, 27.2] |
| 8x strong detection (n=500) | 60.6% [56.2, 64.9] |
| Pooled null FPR / pooled coverage | 6.68% / 91.3% |
| Extrapolated 80%-power scale (model-dependent) | ~1.7e4 planets [1.4, 2.2]e4 |

`python code/audit_highrep.py` recomputes the gate statistics, false-positive
rates, coverage, and power table directly from the raw replicate CSVs and
prints the full comparison against the summary JSONs (which agree exactly).

## Known caveats (also documented in the manuscript)

- Planet count is not independent age information: stellar age is a
  host-level quantity, and the M-dwarf age information is dominated by the
  kinematic proxy v_tan (population-level, not per-star ages); only 12 of 303
  M-dwarf hosts yield usable gyrochronological ages.
- The bootstrap grid upsamples the same ~301 M-dwarf hosts, so host-covariate
  diversity saturates at large scale; the sigma proportional to N^(-1/2)
  extrapolation is a premise, not a grid demonstration.
- Analytical sandwich intervals used inside the power diagnostics are mildly
  anti-conservative (pooled coverage 91.3%, pooled null FPR 6.68%); the
  primary gate verdict is based on empirical injected distributions and is
  unaffected.
- The M-dwarf sample (T_eff < 4200 K) spans late K through late M and does
  not isolate the approximately M3+ regime where the cited valley
  disappearance occurs; spectral-type subdivision is future work.

See the manuscript and `results/colab_audit/AUDIT_REPORT.md` for the complete
discussion.

## Reproducibility limitations (disclosed)

- The executed Colab orchestration variant is archived:
  [`colab/colab_highrep_gate_fast_v3.py`](colab/colab_highrep_gate_fast_v3.py)
  (SHA-256 `fed4824d...d53a69e59`; full hash in `colab/README.md`). It drives
  the reference implementation
  [`colab/colab_highrep_gate.py`](colab/colab_highrep_gate.py) as a library —
  serial synthetic generation in the original rng order, fits dispatched to a
  fork-based process pool, and the per-gate-replicate Hessian (unused by the
  gate statistics) skipped for speed. The audit's fingerprint evidence
  (null-arm first-120 statistics reproduce the original run to <=6e-5;
  identical calibration constants and planet-count sequences) is in
  `results/colab_audit/AUDIT_REPORT.md`.
- `data/berger/` (Berger et al. 2020 supplementary catalog, used only by the
  incomplete isochrone attempt) is third-party material and is not
  redistributed; see `data/README.md`.

## Citation

See [`CITATION.cff`](CITATION.cff). If you use this code, data, or results,
please cite the accompanying manuscript (AASTeX source in
[`manuscript/`](manuscript/)) and this repository.

## License

Code: MIT (see [`LICENSE`](LICENSE)). Curated data files follow the access
conditions of the NASA Exoplanet Archive, Gaia DR3, and Kepler DR25 services.