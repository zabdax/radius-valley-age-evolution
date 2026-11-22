# P3 — Production Gaia quality cuts: verification ledger, recomputation, disclosure numbers

**Working root:** `F:\Research_Radius-Valley-Age-Evolution`
**Frozen project (READ ONLY):** `radius-valley-age-evolution/`
**Writes confined to:** `manuscripts/q1-paper/v2-2026-11/results/audit/` and `.../scripts/`
**Reproduce with:** `C:\Users\MIT\AppData\Local\Programs\Python\Python312\python.exe manuscripts\q1-paper\v2-2026-11\scripts\gaia_production_cuts.py`
**Cut rule (pre-specified in frozen `code/appendix_C_ruwe.py`):** keep `RUWE ≤ 1.4 AND parallax_SNR = parallax/parallax_error > 5`.

---

## 1. Verification ledger — archived post-hoc RUWE/SNR sensitivity result

Sources (all archive-sourced unless marked RECOMPUTED):
- `radius-valley-age-evolution/results/appendix/ruwe_clean_summary.json`
- `radius-valley-age-evolution/results/appendix/ruwe_{M,FGK}_{null,signal}_raw.csv` (200 reps/arm each)
- `radius-valley-age-evolution/results/appendix/confirm_summary.json` (baseline no-cut arm separations)
- `radius-valley-age-evolution/results/appendix/ruwe_gate_all_raw.csv`
- `radius-valley-age-evolution/data/gaia_quality_DR3_v2.csv` and `gaia_quality_DR3.csv`

**Quality-table identity (RECOMPUTED):** `gaia_quality_DR3.csv` (the file the archived appendix actually read) and `gaia_quality_DR3_v2.csv` (the frozen v2 table) have identical columns and 2,420 rows; `source_id`, `ruwe`, `parallax`, `parallax_error` agree to machine precision (max abs diff = 0). The two fetches are interchangeable, so all numbers below are valid for v2.

| Quantity | Cool / M sample | Control / FGK sample | Provenance |
|---|---|---|---|
| Planets before cut | 417 | 2,135 | archive + RECOMPUTED (exact match) |
| Hosts before cut | 301 | 1,588 | archive + RECOMPUTED |
| Planets dropped by `RUWE>1.4` | **35 / 417** (8.4%) | **178 / 2,135** (8.3%) | archive + RECOMPUTED |
| Planets with `SNR≤5` | 0 | 24 (all 24 also `RUWE>1.4`; union = 178) | archive + RECOMPUTED |
| Planets after cut | 382 | 1,957 | archive + RECOMPUTED |
| Hosts dropped | 23 | 143 | archive + RECOMPUTED |
| Median RUWE / median SNR (pre) | 1.0881 / 647.3 | 1.0040 / 73.6 | archive + RECOMPUTED |
| Arm separation, **no cut → after cut** | **−0.19 → −0.25** | **−1.01 → −0.84** | baseline from `confirm_summary.json` meta (`production_legacy`: M_arm −0.19, FGK_arm −1.01) **archive-sourced**; after-cut values RECOMPUTED from raw CSVs |
| After-cut separation (exact) | −0.2545327933 | −0.8368765312 | RECOMPUTED from `ruwe_{sample}_{null,signal}_raw.csv`; matches archive to <1e-12 |
| `real_pct_null` after cut | 0.505 | 0.480 | RECOMPUTED from raw null-arm draws vs archived `real_beta_cleaned`; exact match |
| Verdict | `FAIL(inconclusive)` | `FAIL(inconclusive)` | **archive-sourced**; unchanged from baseline `FAIL` — cut does not rescue or reverse either conclusion |

Arm statistics recomputed from the raw replicate CSVs (median / SD / 2.5% / 97.5% of `beta_hat`), all matching the archive to <1e-12:
- M null: −0.006083 / 0.714328 / −1.697386 / +1.264273; M signal: −0.187903 / 0.746590 / −1.807587 / +0.980250.
- FGK null: −0.044124 / 0.143986 / −0.335670 / +0.238183; FGK signal: −0.164623 / 0.155162 / −0.469515 / +0.145263.

**Verdict on the archived claim:** fully verified. The post-hoc cut drops 35/417 cool and 178/2135 control planets, shifts arm separation from −0.19 to −0.25 (cool) and −1.01 to −0.84 (control), and leaves both verdicts as inconclusive failures.

---

## 2. Production-cut version (cut applied BEFORE the gate) — recomputed

Script: `manuscripts/q1-paper/v2-2026-11/scripts/gaia_production_cuts.py`. It loads the frozen base via `s15.load_control` (Kepler-only filter for FGK, completeness weights, 95th-percentile weight trim), merges `source_id` → Gaia quality, applies the cut, then rebuilds the design on the *cleaned* frame, pins that frame's own mixture components, and refits the real pooled β with the fixed-`su` cluster-robust sandwich SE — i.e. the same estimator/mode/weighting as the archived run.

**ESS definition (identical to frozen `15_stratified_gate.py` line 469):** `ESS = (Σ w)² / Σ w²` on the 95%-trimmed completeness weights.

| Quantity | Cool / M | Control / FGK |
|---|---|---|
| N (planets) before → after cut | 417 → **382** | 2,135 → **1,957** |
| Hosts before → after cut | 301 → **278** | 1,588 → **1,445** |
| ESS before (archive baseline) → after (RECOMPUTED) | 139.51 → **129.40** (92.7% retained) | 1,644.99 → **1,512.54** (92.0% retained) |
| Real pooled β/SE before (RECOMPUTED) | +0.045763 / 0.320965 | −0.029596 / 0.106103 |
| Real pooled β/SE after cut (RECOMPUTED) | **−0.002802 / 0.360170** | **−0.048681 / 0.117196** |
| Archived `real_beta_cleaned` / `real_se_cleaned` | −0.002801691 / 0.360170365 | −0.048680738 / 0.117196406 |
| Agreement | exact (<1e-12) | exact (<1e-12) |

Interpretation: the production cut costs ~8% of planets and ~7–8% of ESS. The real pooled β stays consistent with zero in both samples (|β| < 0.05, well within its SE); this is the same qualitative result as the no-cut production fit. The blank residual signal is therefore **not** an artifact of bad astrometry in the excluded tail.

**What is RECOMPUTED vs CARRIED:**
- RECOMPUTED here: post-hoc quality counts; all arm statistics and separations from the archived raw CSVs; production-cut N/hosts/ESS; production-cut real β and SE; parallax-error quantiles.
- CARRIED FROM ARCHIVE: baseline (no-cut) ESS values (139.5074 / 1644.9874) and baseline production arm separations (−0.19 / −1.01) from `confirm_summary.json`; the null/signal arm Monte Carlo draws (and hence verdict labels) when the script is run at default `--gate-reps 0`. The after-cut arms are *not* a new simulation — they are the archived 200-rep/post-hoc ensemble, which already applied the cut to the base before `synth()`. The script supports `--gate-reps 200` to regenerate those arms from scratch with the production cut (multiprocessing is Windows-guarded); that option was **not** run for this ledger at the orchestrator's request to deliver promptly.

---

## 3. Parallax uncertainty propagation (for footnote/disclosure)

Source: `radius-valley-age-evolution/data/gaia_quality_DR3_v2.csv` (2,420 unique host `source_id`s; the frozen `gaia_hosts_*.csv` carry `parallax` but no error, so `parallax_error` is taken from the quality table). The merged parallax columns agree exactly (max abs diff = 0); no non-positive parallaxes. With `d[kpc] = 1/π[mas]`, `σ_d/d = σ_π/π = parallax_error/parallax`.

| Subset | n | median σ_d/d | 95th pct | 5th pct | max |
|---|---|---|---|---|---|
| All quality-table hosts | 2,420 | **0.00877** | **0.04623** | 0.00071 | 3.504 |
| M hosts | 301 | 0.00149 | 0.01343 | 0.00030 | 0.122 |
| FGK hosts | 2,123 | 0.01014 | 0.04943 | 0.00110 | 3.504 |
| Pass `SNR>5` | 2,400 | 0.00864 | 0.04318 | 0.00070 | 0.186 |
| **Fail `SNR>5` (≤5)** | **20** | **0.34959** | **2.42378** | 0.21993 | 3.504 |

All 20 SNR failures are FGK hosts. **Orchestrator cross-check** (`P3_parallax_orchestrator_check.md`) agrees: median 0.0088 (all) / 0.0082 (both cuts applied), p95 0.0462 (all) / 0.0354 (both cuts), failures median 0.350 / p95 2.42. (Small pass-subset differences reflect whether only SNR or both RUWE and SNR are applied.)

**Disclosure statement:** inverse-parallax distances are adequate at the adopted sample median — fractional distance error ≈ 0.9% for all hosts, ≈ 0.15% for the M-dwarf hosts and ≈ 1.0% for the FGK hosts, i.e. sub-percent to ≈1%, with 95th percentile ≈ 4.6% overall. They are *not* adequate for the 20 low-SNR hosts (median 35%, 95th percentile 242%, where `d = 1/π` is effectively meaningless), and those hosts are exactly what the production `SNR>5` cut removes; no zero-point correction was applied. Safe to state that distance errors are negligible at the sample median and that the low-SNR tail is excluded.

---

## 4. Caveats

1. The 2.5%/97.5% arm quantiles from 200 replicates are Monte-Carlo-limited (binomial SE ≈ 1% in percentile level); the *archive* separations −0.19/−1.01 came from larger ensembles (400/1000 reps) and are not exactly like-for-like with the 200-rep after-cut suite — this is already how the archived appendix reported them.
2. `RUWE>1.4` and `SNR≤5` overlap in the FGK sample (24 low-SNR planets all also have high RUWE), so the total drop (178) is the union, not the sum.
3. Dropped hosts are kinematically unremarkable in FGK (median vtan 34.7 kept vs 34.7 dropped) and slightly colder in M (33.7 vs 27.9 km/s); the cut is not selectively removing the high-vtan tail that drives the proposed age signal.
4. The quality table fetched in both `DR3.csv` and `DR3_v2.csv` has no missing quality rows for the sample hosts (`planets_no_source_id = 0`, `planets_no_quality_row = 0`), so no silent attrition beyond the explicit cut.
5. No full gate rerun with the cut applied at the primary sample size was performed here; arm draws and verdict labels are carried from the archived post-hoc suite (which used the same cut-before-`synth()` procedure). Counts, ESS, and real β are genuine recomputations from the frozen tables.

## 5. Output files
- `manuscripts/q1-paper/v2-2026-11/results/audit/P3_gaia_production_cuts.md` (this file)
- `manuscripts/q1-paper/v2-2026-11/results/audit/gaia_cut_summary.json` (machine-readable ledger; written by the script)
- `manuscripts/q1-paper/v2-2026-11/scripts/gaia_production_cuts.py` (reusable script, Windows-guarded)
