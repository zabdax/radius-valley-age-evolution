# P3 — Exact host demographics of the frozen M-dwarf planet sample

**Date:** 2026-11-21. **Status:** recomputed from read-only frozen data. All numbers in
`P3_host_demographics.csv` are machine-readable; this note explains method and caveats.

## Source and method

- Input: `radius-valley-age-evolution/data/planet_sample_M.csv` (419 planet rows; **303 unique
  `hostname` values**), plus `data/gaia_hosts_M.csv` for the crossmatch inventory only.
- **Host-level** statistics: one value per unique `hostname` (`.groupby("hostname").first()`).
  This avoids double-counting multi-planet systems. `st_teff` is present for all 303 hosts;
  `st_mass` is usable for **282 of 303** hosts (21 have blank/non-usable mass).
- **Planet-level** statistics: every one of the 419 rows counts once. These are the numbers the
  calibrated gate actually weights, so they are reported alongside the host fractions. Two
  denominators are given for planet mass: usable-mass planets (385) and all 419 planets.
- Quartiles: with numpy's default linear interpolation the host `st_teff` quartiles are
  **3417.0 / 3866.5 K**; the task-stated **3415 / 3867 K** corresponds to discrete lower and
  upper hinges (`method="lower"` gives 3415 K, `method="higher"` gives 3867 K). Both conventions
  are tabulated in the CSV so the published choice is explicit. The median is **3683 K** under
  either convention. Range: 2566–4195 K; mass range 0.0898–0.808 M_sun.

## Recomputed values

| Quantity | Host-level | Planet-level (gate-weighted) |
|---|---|---|
| Sample size | 303 hosts | 419 planets |
| Median `st_teff` | 3683 K | — |
| Teff Q1 / Q3 (linear) | 3417.0 / 3866.5 K | — |
| Teff Q1 / Q3 (discrete hinges) | 3415 / 3867 K | — |
| T_eff < 3400 K | 71/303 = **23.43%** | 100/419 = **23.87%** |
| T_eff < 3000 K | 4/303 = **1.32%** | 13/419 = **3.10%** |
| M_* < 0.4 M_sun | 101/282 = **35.82%** | 139/385 = 36.10%; 139/419 = **33.17%** of all planets |
| M_* < 0.3 M_sun | 58/282 = **20.57%** | 83/385 = 21.56%; 83/419 = **19.81%** of all planets |

The task-stated figures (303 hosts, 3683 K median, 23.4% <3400 K, 1.3% <3000 K, 282 usable
masses, 35.8% <0.4 M_sun, 20.6% <0.3 M_sun) are reproduced to the quoted precision. The only
non-trivial convention choice is the quartile interpolation, documented above.

## Interpretation caveat — MUST be stated in the manuscript

The sample **straddles the ~0.3–0.4 M_sun transition** where the radius valley is observed to
fade or is contested. By the host counts, ~36% of hosts with usable mass are below 0.4 M_sun and
~21% are below 0.3 M_sun (planet-weighted: ~33% and ~20% of all planets). Roughly two-thirds of
the sample lies above the transition. Therefore:

- Do **not** claim the whole M-dwarf sample is "valley-present." The valley morphology in this
  mass range is mass-dependent and still debated; a single pooled gate is a deliberate
  analysis choice, not a statement that every host has a valley.
- The FGK control sample is the part of the analysis that carries the +4.14 pp residual
  period-mixing AME; the demographics here characterize the M-dwarf science sample.

## Provenance

- Script: `manuscripts/q1-paper/v2-2026-11/scripts/build_P3_audit_tables.py`
- Output table: `manuscripts/q1-paper/v2-2026-11/results/audit/P3_host_demographics.csv`
- No frozen file was modified.
