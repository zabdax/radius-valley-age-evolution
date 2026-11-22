# v2 reproducibility bundle

This directory contains the artifacts required to reproduce the **v2 (Q1)**
scientific results of the radius-valley age-evolution study: the analysis
scripts actually executed, the authoritative result tables and JSON summaries,
the publication figures, and the provenance records.

It deliberately contains **no manuscript source**: `main.tex`, `references.bib`,
the compiled PDF, AASTeX style files, and the paper-writing/process documents
are not distributed here.

## Layout

- `scripts/` — v2 analysis and figure-regeneration scripts. These are the exact
  drivers used for the reported results. Path constants were rewritten to be
  portable through `scripts/_paths.py`; set the `RV_REPO_ROOT` environment
  variable if the repository is relocated.
- `results/science/` — authoritative numerical outputs (effect-size derivation,
  TIME Table crossmatch and Spearman check, period-mixing tests, full 10-point
  power-versus-effect curve and checkpoints). `*_FINAL` files are authoritative.
- `results/audit/` — data-level audit outputs (Gaia production cuts, host
  demographics, multiplicity budget, archived-run verification).
- `figures/` — the seven publication figures.
- `provenance/` — baseline checksums, provenance lock, Zenodo release manifests,
  and `BUNDLE_SHA256.txt` covering every file in this bundle.

## Verification

```bash
cd repro/v2
sha256sum -c provenance/BUNDLE_SHA256.txt
```

The frozen input CSVs live in the repository's top-level `data/` directory and
are covered by `data/checksums.sha256`. The gate and mixture machinery they
depend on lives in top-level `code/`.

## Notes

The TIME Table gyrochronological ages are used only for ordinal proxy
validation (`timetable_spearman.json`; Spearman rho = 0.372, p = 0.0023,
n = 65). They are not used in the primary calibrated gate because they fall
outside the Angus+19 validity domain.
