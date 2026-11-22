# Priority 0 — Provenance lock (read-only)

Date: 2026-11-21. Workspace: `manuscripts/q1-paper/v2-2026-11/`.

## Policy executed

- All frozen code, data and prior results under `radius-valley-age-evolution/` are read-only. Nothing there was modified for this build. (New first-party artifacts added earlier in the run, `results/proxy_validation_fgk_recomputed.json` and its script under the research run temp folder, are additions, not edits.)
- `paper-workspace/` and `manuscripts/revised-package/` were not edited.
- No raw data were copied into the manuscript folder. Only the manuscript baseline (`.tex`, `.bib`, class, bst, figures) was staged into the versioned folder.
- `BASELINE_SHA256.txt` records SHA-256 hashes of the 10 staged baseline files.

## Verification performed

- `radius-valley-age-evolution/data/checksums.sha256`: 31 entries, 31 verified, 0 mismatches or missing files.
- Stream-fingerprint evidence (from `results/colab_audit/AUDIT_REPORT.md`): the 1,000-replicate run's null-arm first-120 subset reproduces the original 120-replicate run (median −0.013282 vs −0.013343, difference 6.1e−5; q025 difference 2.4e−7; q975 difference 2.3e−6; percentile rank 0.4333 identical; replicate-level draws match an independent local implementation). The generative stream is verified unchanged.

## Python environment

Scientific Python located at `C:\Users\MIT\AppData\Local\Programs\Python\Python312\python.exe` with numpy, pandas, scipy, matplotlib, astropy and gyrointerp available. The hyperresearch venv Python lacks these. All new computations use the scientific Python.
