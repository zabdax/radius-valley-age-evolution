# High-replicate Monte Carlo experiment on Google Colab

`colab_highrep_gate.py` re-runs the calibrated-gate Monte Carlo (and the
15-cell power grid) at **1,000 replicates per gate arm** and **500 replicates
per power-grid cell** — the counts used in the published manuscript — using
the identical model specification and generative process as
`../code/13_calibration_fgk.py` and `../code/11_power_analysis.py`.
It adds no new methodology; it only increases the replicate count and dumps
raw replicate-level results.

## Usage (Colab)

1. Runtime > Change runtime type > **CPU** (GPU does not help).
2. Build and upload the data bundle:
   - run `python make_data_bundle.py` locally (creates `colab_gate_data.zip`,
     ~0.4 MB from the five frozen CSVs in `../data/`), then upload it to
     `/content/`; or set `COLAB_GATE_DATA` to a Drive folder containing the
     five CSVs.
3. Upload `colab_highrep_gate.py` and run `!python colab_highrep_gate.py`.
   Optional: `--smoke` first (tiny plumbing test, ~10 min).
4. Outputs land in `/content/colab_gate_output/` (also zipped): the two raw
   replicate CSVs, the two summary JSONs, run manifests, and resumable
   checkpoints. If the session dies, re-upload the output zip and re-run —
   the script resumes from the last checkpoint with the exact rng streams.
5. Copy the four result files back into `../results/highrep/` (this is the
   authoritative results set used by the manuscript and by
   `../code/audit_highrep.py`).

The script prints a measured runtime estimate before committing to the full
run (~5-7 h for the gate; the grid budget is configurable and the gate always
completes first). Full configuration options are documented in the script
docstring and in `../docs/USAGE_COLAB.md`.

## Reproducibility note

The as-executed orchestration variant is archived in this folder:
`colab_highrep_gate_fast_v3.py`
(SHA-256 `fed4824dd0fde7c2cc13077f19ba866577ea551d8dfaf1d01806d56d53a69e59`).
It imports the reference implementation (`colab_highrep_gate.py`) as a
library, overrides only the replicate counts and output paths, generates
synthetic data serially in the original rng order, and dispatches the fits
to a fork-based process pool; for gate replicates it skips the per-replicate
Hessian evaluation, which the gate statistics never use. Fingerprint evidence
that the generative process was unchanged (identical null-arm first-120
statistics, identical calibration constants and planet count sequences) is
archived in `../results/colab_audit/AUDIT_REPORT.md`. The reference
implementation (`colab_highrep_gate.py`) remains the documented, serial
version of the same protocol.