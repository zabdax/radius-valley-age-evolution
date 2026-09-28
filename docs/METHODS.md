# Methods (summary)

This is a **calibrated statistical gate**, not a single fit. Every age-evolution
estimate must survive a null/signal injection test before it can be claimed.

## 1. Samples

- **M dwarfs:** 419 planets / 303 hosts, `st_teff < 4200 K`, Kepler transits,
  `0.5 < R_p < 4 R_earth`, `P < 100 d` (`code/01_download.py`, NASA Exoplanet Archive TAP).
- **FGK control:** 2,880 planets / 2,150 hosts, `4200 < st_teff < 7000 K`, same cuts.
- **Gate subset:** FGK Kepler-only, completeness-weighted, N = 2,135 planets,
  ESS = 1,645 (DR25 CDPP trimmed weights, `code/12_completeness.py`).

## 2. Age proxies

| Proxy | M dwarfs | FGK | Code |
|---|---|---|---|
| Kinematic `v_tan` (+ U,V,W) | 301/303 hosts | 2,123/2,150 | `07_kinematic_ages.py` |
| Gyrochronology (Angus+2019 via `gyrointerp`) | 12 usable of 66 with rotation | 79 usable of 88 | `03_gyro_ages.py` |
| Literature rotation crossmatch | McQuillan+2014 | same | `05_literature_rotations.py` |
| Isochrone (Berger+2020 attempt) | incomplete, reported as limitation | same | `B1/B2` |

Planet count is **not** independent age information: age is host-level, and the
M-dwarf information is dominated by the population-level `v_tan` proxy.

## 3. Model

Working-independence mixture likelihood over radius-valley membership with
host-level cluster-robust (sandwich) SEs (`code/09_hierarchical_model.py`).
Mixture components pinned for stratified variants. Coefficient reported per
IQR of `v_tan`: `beta_age/IQR`.

## 4. Calibrated gate (authoritative: 1,000 reps/arm)

1. Fit real data → `beta_real = -0.0296 ± 0.1061`.
2. Inject null (`beta = 0`) and literature-sized signal into 1,000 synthetic
   replicates each, refit identically.
3. PASS requires real to sit **inside** signal-95% **and outside** null-95%.
4. Result: 45.6th %ile of null, 81.2nd %ile of signal → **FAIL (inconclusive)**.
   Both hypotheses remain admissible. A wrong-sign binned age signal is a
   **period-mixing artifact** (old hosts at longer P where sub-Neptunes dominate).

## 5. Power grid (500 reps/cell, 15 cells)

Bootstrap grid over sample-size × effect (null / SWEET literature-sized /
strong 2×). At 8× current inventory (~3,337 planets) literature-sized detection
is **23.2% [19.6, 27.2]**, strong **60.6% [56.2, 64.9]** — no tested scale
reaches 80%. Gaussian extrapolation places 80% power at ~1.7e4 planets
([1.4, 2.2]e4), a **model-dependent guide**, not a demonstrated requirement.

## 6. Robustness

One-per-host check, stratified (vmag/dist) gate, RUWE/parallax-S/N gate,
arm-delta and beta-distribution appendices, independent-host confirmation at
production scale (M 1000/arm, FGK 400/arm) — verdicts unchanged. Analytical
sandwich intervals are mildly anti-conservative (pooled FPR 6.68%, coverage
91.3%); the gate verdict uses empirical injected distributions and is unaffected.
