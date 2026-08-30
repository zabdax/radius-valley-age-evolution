# Colab high-replicate Monte Carlo — usage note

**Script:** `colab_highrep_gate.py` (standalone; also runs locally).
**Data upload:** `colab_gate_data.zip` (~350 KB, already built in this
folder) — contains the 5 frozen CSVs: `planet_sample_FGK.csv`,
`hosts_kinematics_FGK.csv`, `completeness_weights_FGK.csv`,
`planet_sample_M.csv`, `hosts_kinematics_M.csv`.

## What it runs
The project's existing calibrated-gate Monte Carlo with the model
specification **unchanged** (per the Issue C investigation):
working-independence mixture (`fixed_su=True`, sigma_u=0), host-level
cluster-robust sandwich SEs, host-bootstrap injections with one shared
host effect per drawn host (sigma_u=1.5). Gate: FGK Kepler-only, trim
95th pct (N=2,135, ESS~1,645), arms null / signal (B_INJ=-0.14 per SD),
seed 20260826 — identical generative process to `code/13_calibration_fgk.py`.
Grid: 15 cells, per-cell seeds and rng order identical to
`code/11_power_analysis.py`. Only the replicate counts increase.
No new methodology; no conclusions; does not touch the pipeline or PDF.

## Colab steps
1. Runtime > Change runtime type > **CPU** (GPU does not help).
2. Upload `colab_gate_data.zip` to `/content/`.
3. Upload (or paste) `colab_highrep_gate.py`, then run:
   `!python colab_highrep_gate.py`
   Optional first: `!python colab_highrep_gate.py --smoke` (~10 min
   plumbing test, writes to `colab_gate_output_smoke/`).
4. Edit the CONFIG block to change counts/budget if desired:
   `GATE_REPS=1000` (per arm), `RUN_GRID=True`, `GRID_REPS=500` (per
   cell), `GRID_BUDGET_HOURS=5.0` (grid-only budget; the gate is never
   shortened), `CHECKPOINT_EVERY=25`.

## Expected runtime (measured warm-up; Colab CPU is similar)
- Gate: ~9–13 s per replicate at N=2,135 -> **2 x 1000 reps ~ 5–7 h**.
  This is the priority and always completes first.
- Grid: ~6–7 ms per planet (synth+fit+sandwich) -> full 15 cells x 500
  reps ~ 14–18 h total, so the default 5 h budget completes roughly the
  0.5x–4x cells and part of 8x. Unfinished cells are skipped cleanly;
  **re-run the script in a follow-up session (re-upload
  `colab_gate_output.zip` next to the script or unzip it into
  `/content/`) and it RESUMES from checkpoints with the exact rng
  streams** until all cells reach `GRID_REPS`.
- The script prints its own measured estimate before committing.

## What comes out (`/content/colab_gate_output/`, zipped as
`colab_gate_output.zip` beside it)
- `task2_calibration_kepler_only_trim0.95_rep1000.json` — same schema as
  `results/task2_calibration_kepler_only_trim0.95.json` (superset with
  provenance keys). Drop-in.
- `task2_gate_rep1000_replicates_raw.csv` — RAW rows:
  `arm,rep,beta_inj,alpha,beta_hat_per_iqr,n_planets`.
- `power_analysis_rep500.json` — same list-of-cells schema as
  `results/power_analysis.json` (`reps` records the actual count per
  cell). Drop-in.
- `power_grid_rep500_replicates_raw.csv` — RAW rows: `cell,scale,effect,
  beta_inj,rep,beta_hat_per_sd,se_per_sd,detect,cover,n_planets`.
- `checkpoints/` (resume state) and `manifest.json` (run provenance,
  including any grid cells skipped by the budget).

## Import back into the project
Copy the two `*_replicates_raw.csv` and the two summary JSONs into
`mdwarf_radius_valley/results/`. The summary JSONs are supersets of the
existing formats and drop in without reformatting; rename (drop the
`_repNNN` suffix) if they should supersede the local versions.
