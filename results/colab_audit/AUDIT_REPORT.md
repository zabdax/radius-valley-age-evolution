# HIGH-REPLICATE MONTE CARLO AUDIT REPORT
**Project:** The Radius Valley's Age Evolution Around M Dwarfs — calibrated-gate pipeline
**Auditor:** Monte Carlo Results Auditor | **Date:** 2026-08-30
**Inputs:** `colab_gate_output_fast_v3.zip` and `mdwarf_radius_valley_highrep_fast_v3-20260830T080159Z-1-001.zip` (byte-identical contents, SHA-256 verified).
**Ground truth:** raw replicate CSVs; JSONs used only as cross-checks.
**Audit code:** `results/colab_audit/audit_v3.py` (full console log `audit_out.txt`, machine-readable `audit_summary.json`).
**Purpose:** hand-off to the manuscript-revision agent. Evidence-first; old manuscript not assumed correct.

---

## 1. EXECUTIVE SUMMARY

The high-replicate run (1,000/arm gate; 500/cell grid, all 15 cells complete) **confirms the manuscript's core conclusions**: the FGK gate remains **FAIL (inconclusive)**, the 8x scale remains far below 80% power (updated **23.2%**, up from 18.7%), and "8x is a floor, not a sufficient target" stands. Three things change materially:

1. The previously reported FPR **non-monotonicity is confirmed as finite-replicate noise**, but a **mild, persistent, scale-independent FPR inflation** (pooled 6.68% vs 5.0%, p = 1e-4) replaces the old "recovers toward nominal with scale" narrative.
2. Analytic-interval **coverage is significantly below nominal** (pooled 91.3% vs 95%).
3. The extrapolated requirement updates from ~2.3e4 to **~1.7e4 planets** [1.4, 2.2]e4 — still "of order 2e4" — with a newly quantified caveat: the grid's own SEs **saturate with host diversity** (upsampled cells reuse the same 301 hosts), so sigma ∝ N^(-1/2) scaling is a premise the grid cannot itself demonstrate.

**Classification: B — the run strengthens the manuscript but requires several substantive (numerical/narrative) revisions.**

---

## 2. DATA INTEGRITY — PASS

**Files:** all five outputs hash-identical across both archives. `manifest.json`: executed script `colab_highrep_gate_fast_v3.py`, gate 1000/arm, grid 500/cell, workers=2 (process pool), batch=8, `grid_skipped_cells: []`, 6.03 h elapsed, `model_spec_unchanged: true`.

**Gate CSV (2,000 rows):** null 1000/1000, signal 1000/1000; rep IDs contiguous 0-999; zero duplicates; zero NaN/Inf; single alpha per arm (-0.708348); n_planets 2135.5 +/- 29.0 [2039-2233] — the host-bootstrap signature of the frozen N=2,135 sample.

**Grid CSV (7,500 rows):** all 15 cells 500/500 (budget not binding); contiguous; zero duplicates; zero NaN/Inf; detect/cover all in {0,1}; per-cell mean n_planets tracks scale x 417 to <0.3% (208.2 / 417.3 / 834.5 / 1667.4 / 3337.8).

**JSON cross-check:** every gate and grid summary statistic recomputed from raw agrees with the JSONs to exactly 0.0e+00.

**Fidelity caveats (reported, not repaired):**
- The executed `colab_highrep_gate_fast_v3.py` was **not included** in the archives — line-level code audit impossible; mitigated statistically (below). Remediate before submission by archiving the script with its hash.
- **Stream fingerprint:** the new run's *null-arm first-120 subset* reproduces the original 120-replicate run essentially exactly (median -0.013282 vs -0.013343, d=6.1e-5; q025 d=2.4e-7; q975 d=2.3e-6; percentile rank 0.4333 **identical**); rep-level draws match an independent local validation run to 3.3e-4 with identical alpha and n_planets sequences. The signal-arm first-120 differs (median -0.146 vs -0.178) **by design**: its rng stream position legitimately begins after 1,000 rather than 120 null synths. Fits differ from local runs only at optimizer/floating-point tolerance (<=3.3e-4 vs replicate SD ~0.12). Conclusion: generative process unchanged.

---

## 3. GATE RESULTS (recomputed from raw, n=1000/arm; JSON agrees to 0.0e+00)

Real-data coefficient: **beta/IQR = -0.0296 +/- 0.1061** (JSON; agrees with independent local fit to 4e-7 / 3e-8).

| Quantity | Null-injected (n=1000) | Signal-injected (n=1000) |
|---|---|---|
| median beta-hat | -0.0186 (bootstrap 95% [-0.0256, -0.0113]) | -0.1365 (bootstrap 95% [-0.1473, -0.1270]) |
| q025 / q975 | -0.2476 / +0.2246 | -0.3898 / +0.0836 |
| mean +/- SD | -0.0140 +/- 0.1169 | -0.1412 +/- 0.1248 |
| MAD | 0.0759 | 0.0822 |
| min / max | -0.355 / +0.361 | -0.611 / +0.269 |
| real beta percentile rank | **45.6%** (binom95 [42.5, 48.7]) | **81.2%** (binom95 [78.6, 83.6]) |

- **A. Null compatibility:** yes — P(null-injected >= real) = **0.544**; real coefficient at the 45.6th percentile.
- **B. Signal compatibility:** yes — real beta inside the signal 95% envelope ([-0.3898, +0.0836]), at its 81.2th percentile.
- **C. Verdict under the exact existing rule:** distinguishable-from-null = False -> **FAIL (inconclusive)**. Unchanged.
- **D. Stability vs the 120-replicate run:** null-arm changes negligible (old values inside new CIs). Signal arm: median -0.178 -> -0.137, q025 -0.355 -> -0.390, q975 +0.068 -> +0.084, real percentile 90.8% -> 81.2%. The signal-median shift is ~2.8 sigma of combined MC error — **modest**, reflecting the small original sample plus the necessarily different stream segment; no conclusion changes, but the manuscript's signal-envelope numbers must be replaced.

---

## 4. MONTE CARLO PRECISION (distinct from model/statistical uncertainty)

- Gate percentile ranks: MC error now +/-2.5pp (null: binom95 [42.5, 48.7]) vs +/-4.5pp at n=120 (~1.8x precision gain).
- Gate medians: bootstrap SE ~0.0037 (null) / 0.0052 (signal) vs ~0.0128/0.0143 at n=120.
- Proportions (n=500/cell): binomial SE ~1.0-2.2pp depending on rate; e.g. 8x SWEET 23.2% +/- 1.9pp.
- All astrophysical/model uncertainty (attenuation, sample noise, covariate distribution) is unchanged by replication and remains dominant; replication only shrinks the Monte Carlo layer. The two layers are strictly separated throughout.

---

## 5. NULL FPR CALIBRATION — narrative must change

| Scale | FPR (n=500) | CP 95% CI | vs nominal 5% | binomial p |
|---|---|---|---|---|
| 0.5x | 7.2% (36/500) | [5.1, 9.8] | marginally high | 0.0196 |
| 1x | 8.2% (41/500) | [6.0, 11.0] | **high** | 0.0015 |
| 2x | 5.2% (26/500) | [3.4, 7.5] | consistent | 0.447 |
| 4x | 5.6% (28/500) | [3.8, 8.0] | consistent | 0.296 |
| 8x | 7.2% (36/500) | [5.1, 9.8] | marginally high | 0.0196 |

- **The 8x anomaly is gone:** new 8x (36/500) is identical to 0.5x (Fisher p = 1.000) and indistinguishable from every other scale (p = 0.24-1.00; vs pooled others p = 0.617). The old 8x = 7.5% vs 4x = 2.5% contrast was **finite-replicate noise** — the old 4x cell was the fluctuation (new 4x = 5.6%). The Issue-A anticipation ("Monte Carlo noise at 120 replicates") is confirmed.
- **But "recovers toward nominal" is dead:** pooled FPR = **167/2500 = 6.68%**, binomial P(X >= 167 | n=2500, p=0.05) = **1e-4**. The inflation is mild (~1.7pp), roughly scale-independent, and persists at 8x. It is consistent with the documented anti-conservative sandwich SEs — a *CI calibration* property, **not** a gate problem (the gate uses empirical null distributions, not analytic CIs).

---

## 6. COVERAGE — mild persistent anti-conservatism

Pooled: **6850/7500 = 91.3%** vs nominal 95% (P(X <= 6850) < 1e-4). New per-cell range: **85.2-95.4%** (old claim "83-98%" no longer applies). Strong-effect cells 85.2-91.0% (all CP-upper < 0.95, significantly anti-conservative); SWEET cells 90.8-95.4%; null cells 92.8-94.8%. Worst cell: 8x strong (85.2% [81.8, 88.2]). Root cause shared with the FPR inflation: analytic sandwich CIs are mildly too narrow (~2-4pp) at a roughly constant level, worst for the strongest effects. The manuscript's "empirical coverage spans 83-98%" and "broadly well-calibrated once sample asymptotic properties stabilize" must be revised — the miscalibration does **not** vanish with scale.

---

## 7. POWER GRID (recomputed from raw; n=500 per cell)

| Scale | ~N planets (mean raw) | Null FPR | SWEET detection | Strong detection |
|---|---|---|---|---|
| 0.5x | 208 | 7.2% [5.1, 9.8] | 9.0% [6.6, 11.9] | 16.2% [13.1, 19.7] |
| 1x | 417 | 8.2% [6.0, 11.0] | 9.4% [7.0, 12.3] | 27.0% [23.2, 31.1] |
| 2x | 834 | 5.2% [3.4, 7.5] | 13.6% [10.7, 16.9] | 38.8% [34.5, 43.2] |
| 4x | 1667 | 5.6% [3.8, 8.0] | 18.8% [15.5, 22.5] | 50.2% [45.7, 54.7] |
| 8x | 3337 | 7.2% [5.1, 9.8] | **23.2% [19.6, 27.2]** | **60.6% [56.2, 64.9]** |

SWEET and strong curves are monotone; null FPR shows no scale trend beyond the uniform mild inflation. **New finding (SE-scaling audit):** median replicate SE (SWEET cells) = 0.221 / 0.162 / 0.127 / 0.105 / 0.093 — ratio to the 8x cell = 2.37 / 1.73 / 1.36 / 1.12 / 1.00, versus the sqrt-scaling expectation 4.00 / 2.83 / 2.00 / 1.41 / 1.00. Upsampling reuses the same 301 real M-dwarf hosts (~8 copies each at 8x), so *host-covariate diversity saturates*: the 8x grid cell behaves like ~3x in independent-host information. (Below 1x, sqrt-scaling holds — 0.5x matches — confirming the mechanism.) Consequence: the grid **understates** the power a genuinely independent-host sample of the same size would have; the host-reuse channel pushes required-N *down*, while the existing unmodeled-survey-differences caveat pushes *up*.

---

## 8. ~3,300-PLANET INTERPRETATION — stands, numbers update

1. Updated 8x SWEET power: **23.2% [19.6, 27.2]** (was 18.7%).
2. Updated 8x strong power: **60.6% [56.2, 64.9]** (was 64.0%).
3. "8x is a floor, not a sufficient target" — **remains valid**; 23.2% is still far below 80%, and no curve reaches the 80% line (max 60.6%).
4. 3,300 should **remain** in the abstract/discussion — it is still the largest tested scale and the concrete anchor.
5. Wording changes: "~19%" -> "~23%" everywhere; Section 5.3's "first exceeds ~15%" -> "reaches only 23.2%"; Figure-4 caption value updates to 23.2% (its structural claim "no curve reaches the dashed 80% line" remains true).

---

## 9. ~2e4 EXTRAPOLATION — defensible as labeled, with updated numbers and one new caveat

Gaussian logic check: power p implies standardized effect x = Phi^-1(p) + 1.96; 80% power needs x80 = 2.80; N proportional to x^-2. Internally consistent. Updated inputs:

- p8 = 0.232 -> x = 1.228, factor = 5.21 -> **N80 ~= 17,400 (~42x current)**;
- CP interval on p8 [0.196, 0.272] -> factor 4.30-6.45 -> **N80 in [1.43e4, 2.15e4]**;
- old p8 = 0.187 gave factor 6.86 -> N80 ~= 2.29e4. The old point estimate sits at the upper edge of the new interval; "of order 2e4" survives as an order-of-magnitude statement, with the point estimate better written as **~1.7e4 (range ~1.4-2.2e4)** and the "additional factor ~5-12" parenthetical updated to **~4.3-6.5**.
- **Critical caveat (new, must be added):** the sigma proportional to N^(-1/2) premise is contradicted *within the grid itself* (SE saturation, Section 7) because upsampled cells cannot create host diversity. The extrapolation implicitly assumes a future sample whose host diversity keeps growing — plausible for a real PLATO-era sample, but undemonstrated here. EMPIRICAL ("23.2% observed in the simulated 8x cell") must not be conflated with EXTRAPOLATED ("~1.7e4 under a Gaussian, independent-host scaling assumption").
- The Section 4.3/5.3 "figures are lower bounds" direction claim is contradicted by the host-diversity channel (which favors *lower* required-N) and must be replaced by the two-channel caveat (host-diversity saturation down; unmodeled survey properties up).

---

## 10. NUMERICAL / OPTIMIZATION DIAGNOSTICS — clean

Across 9,500 replicate fits (2,000 gate + 7,500 grid): zero NaN/Inf; zero exact-duplicate beta-hat values; |beta| > 1 in exactly 2 of 7,500 grid fits (both at 0.5x, the smallest cell — expected heavy-tail behavior); no SE < 0.02; SE > 1 in one cell (0.5x null); SE distributions plausible and shrink monotonically with scale (8x: [0.068, 0.162]); gate |beta| <= 0.611. No evidence of non-convergence, boundary-hitting, or singular Hessians observable at replicate level (component parameters are not exported — the one unobservable). Fit-level platform differences (Colab vs local) are <=3.3e-4 per replicate — two orders of magnitude below replicate SD.

---

## 11. OLD vs HIGH-REP COMPARISON

| Quantity | Old | New high-rep | Difference | Interpretation |
|---|---|---|---|---|
| Null median | -0.0133 | -0.0186 | -0.0052 | negligible (old inside new boot CI) |
| Null q025 | -0.2481 | -0.2476 | +0.0005 | negligible |
| Null q975 | +0.2214 | +0.2246 | +0.0032 | negligible |
| Real pct in null | 43.3% | 45.6% | +2.3pp | negligible (within binomial CI) |
| Signal median | -0.1784 | -0.1365 | +0.0419 | modest (~2.8 sigma; superseded) |
| Signal q025 | -0.3551 | -0.3898 | -0.0347 | modest |
| Signal q975 | +0.0681 | +0.0836 | +0.0156 | modest |
| Real pct in signal | 90.8% | 81.2% | -9.6pp | modest; conclusion unchanged |
| Gate verdict | FAIL (inconclusive) | **FAIL (inconclusive)** | — | unchanged |
| 1x SWEET power | 9.3% (150) | 9.4% (500) | +0.1pp | negligible |
| 2x SWEET | 14.7% | 13.6% | -1.1pp | negligible (within CI) |
| 4x SWEET | 14.7% | 18.8% | +4.1pp | modest (old was low outlier) |
| **8x SWEET** | **18.7%** | **23.2%** | **+4.5pp** | material — update everywhere |
| 8x null FPR | 7.5% | 7.2% | -0.3pp | negligible; anomaly resolved |
| 8x strong | 64.0% | 60.6% | -3.4pp | modest (64.0 inside new CI) |
| 1x null FPR | 11.7% | 8.2% | -3.5pp | old value was high outlier |
| Pooled null FPR | — | 6.68% (p=1e-4) | — | new finding: persistent mild inflation |
| Pooled coverage | "83-98%" | 91.3% (85.2-95.4%) | — | new finding: below nominal |

---

## 12. MANUSCRIPT CLAIM AUDIT

- **Abstract** — MODIFY: "~19%" -> "~23%" (single instance); everything else stands.
- **Introduction** — KEEP (no power numbers).
- **Methods 3.3** — KEEP; ADD one sentence: replicate counts updated to 1,000/arm and 500/cell (Colab, process-pool fitting, serial generation; generative stream fingerprint-verified against the original run).
- **Results 4.1** — MODIFY: replace the gate table (null median -0.019, CI [-0.248, +0.225], pct 45.6%; signal median -0.137, CI [-0.390, +0.084], pct 81.2%); verdict line unchanged.
- **Results 4.3** — MODIFY substantially: new power table (n=500/cell); **rewrite the FPR paragraph** — delete "drop toward nominal... recovering from 11.7% to 2.5%" and the 8x-vs-4x non-monotonicity discussion; state the 8x anomaly is resolved as MC noise (7.2% vs 0.5x's 7.2%, Fisher p=1.00) and that a **mild persistent inflation** exists (pooled 6.68%, p=1e-4; 1x 8.2%, p=0.0015), attributed to the anti-conservative sandwich CIs and explicitly non-affecting the gate; UPDATE coverage ("83-98%" -> pooled 91.3%, range 85.2-95.4%, significantly below nominal, gate unaffected); UPDATE the concluding paragraph: 18.7% -> 23.2%, factor range ~4.3-6.5, N80 ~1.7e4.
- **Discussion 5.1** — MODIFY: "~19%" -> "~23%".
- **Discussion 5.3** — MODIFY: 23.2%; N80 ~1.7e4 [1.4, 2.2]e4; REPLACE the "lower bounds" direction claim with the two-channel caveat; keep "8x achievable with PLATO-type yields but leaves the question unresolved."
- **Conclusion** — MODIFY: "~19%" -> "~23%"; "of order 2e4" may stay with the model-dependent label.
- **Figure 4 caption** — MODIFY: "18.7%" -> "23.2%"; "no curve reaches the dashed 80% line" remains true (max 60.6%).
- **Tables** — REPLACE the 4.3 power table and 4.1 gate table with the Section 15 values.

## 13. TOP 5 REFEREE RISKS

1. **"Your null FPR never reaches nominal."** Evidence: pooled 6.68% (167/2500, p=1e-4); 1x = 8.2% (p=0.0015); 8x = 7.2% (p=0.0196). Severity: **moderate** (gate untouched — it is distribution-based). Response: rewrite the 4.3 calibration paragraph; note the direction is anti-conservative, which is conservative for the inconclusive verdict.
2. **"Your analytic CIs under-cover by ~4pp."** Evidence: pooled coverage 91.3% vs 95% (p<1e-4); all strong cells significantly anti-conservative. Severity: **moderate**. Response: same root cause; update text; state explicitly that PASS/FAIL uses empirical injection distributions.
3. **"Your 8x cell is not an 8x sample."** Evidence: SE ratios vs 8x = 2.37/1.73/1.36/1.12/1.00 vs sqrt-scaling 4.00/2.83/2.00/1.41/1.00 (301 hosts reused ~8x). Severity: **moderate** — undercuts the sigma proportional to N^-1/2 extrapolation *as demonstrated* and reverses the "lower bounds" direction. Response: add the caveat; keep N80 explicitly model-dependent; present both caveat channels with opposite signs.
4. **"Your signal envelope moved; why trust either?"** Evidence: signal median -0.178 (120r) -> -0.137 (1000r), ~2.8 sigma; real percentile 90.8% -> 81.2%. Severity: **minor-to-moderate**. Response: the 1,000-replicate estimate supersedes; the old prefix was a different, equally valid stream segment; verdict invariant under both.
5. **"The executed script differs from the archived one and was parallelized."** Evidence: manifest names `colab_highrep_gate_fast_v3.py` (2-worker pool); script not archived. Severity: **minor** (fingerprint: null-arm first-120 reproduces the original to <=6e-5 with identical percentile). Response: archive the executed script with its hash; cite the fingerprint check in the data-availability statement.

---

## 14. FINAL SCIENTIFIC VERDICT

- **Gate:** remains **FAIL (inconclusive)**. Updated: beta = -0.0296 +/- 0.1061 per IQR; null median -0.0186, envelope [-0.2476, +0.2246], real percentile 45.6% [42.5, 48.7], P(null >= real) = 0.544; signal median -0.1365, envelope [-0.3898, +0.0836], real percentile 81.2% [78.6, 83.6].
- **Calibration:** null distribution well-behaved and compatible with the real coefficient; analytic FPR mildly but persistently inflated (pooled 6.68%, p=1e-4) — a CI property, not a gate property; no 8x-specific anomaly.
- **Power:** SWEET 9.0 / 9.4 / 13.6 / 18.8 / **23.2%**; strong 16.2 / 27.0 / 38.8 / 50.2 / **60.6%** at 0.5x-8x.
- **Sample-size interpretation:** ~3,300 correctly described as a floor (23.2% << 80%); ~2e4 defensible **only** as a labeled model-dependent extrapolation — updated point ~1.7e4 [1.4, 2.2]e4, host-diversity saturation caveat added.
- **Classification: B — "High-replicate run strengthens the manuscript but requires several substantive revisions."** Every headline conclusion (inconclusive gate, floor-not-target, order-of-magnitude path forward) is confirmed with tighter precision, but the FPR-recovery narrative, the coverage claim, the 8x power value, the signal-envelope numbers, and the extrapolation's uncertainty/caveat framing all change enough that publishing the old text would be incorrect.

## 15. EXACT NUMBERS TO INSERT INTO THE MANUSCRIPT

| # | Value | n | Source | Type |
|---|---|---|---|---|
| 1 | Real beta = **-0.0296 +/- 0.1061** per IQR | — | task2 JSON (matches local fit to 4e-7) | JSON (cross-checked) |
| 2 | Null median **-0.0186**; q025/q975 **-0.2476 / +0.2246** | 1000 | gate raw CSV | RAW |
| 3 | Real percentile in null **45.6%** (binom95 42.5-48.7); P(null >= real) **0.544** | 1000 | gate raw CSV | RAW |
| 4 | Signal median **-0.1365**; q025/q975 **-0.3898 / +0.0836** | 1000 | gate raw CSV | RAW |
| 5 | Real percentile in signal **81.2%** (binom95 78.6-83.6) | 1000 | gate raw CSV | RAW |
| 6 | Verdict **FAIL (inconclusive)** | — | gate raw CSV (recomputed) | RAW |
| 7 | Power (null/SWEET/strong): 0.5x **7.2/9.0/16.2%**; 1x **8.2/9.4/27.0%**; 2x **5.2/13.6/38.8%**; 4x **5.6/18.8/50.2%**; 8x **7.2/23.2/60.6%** | 500/cell | grid raw CSV | RAW |
| 8 | 8x SWEET CP95 **[19.6, 27.2]%**; 8x strong CP95 **[56.2, 64.9]%** | 500 | grid raw CSV | RAW |
| 9 | Pooled null FPR **6.68%** (167/2500; p=1e-4); 1x **8.2%** (p=0.0015) | 2500 | grid raw CSV | RAW |
| 10 | Pooled coverage **91.3%** (6850/7500; p<1e-4 vs 95%); range **85.2-95.4%** | 7500 | grid raw CSV | RAW |
| 11 | SE-saturation ratios (SWEET, vs 8x): **2.37/1.73/1.36/1.12/1.00** (sqrt-scaling: 4.00/2.83/2.00/1.41/1.00) | 500/cell | grid raw CSV | RAW |
| 12 | Extrapolated N80 ~= **1.7e4 planets (~42x current)**, range **[1.4, 2.2]e4** (~34-52x) | — | derived from #8 via Phi-inverse scaling | EXTRAPOLATED |
| 13 | Mean synthetic n_planets per cell: 208.2 / 417.3 / 834.5 / 1667.4 / 3337.8 (matches scale x 417) | 500/cell | grid raw CSV | RAW |

**Cannot be concluded from the provided files:** (a) line-level audit of `colab_highrep_gate_fast_v3.py` (absent from both archives — request for the reproducibility record); (b) per-host planets-per-realization distributions and mixture-component/boundary diagnostics (not exported at replicate level); (c) Hessian conditioning per fit (not exported). None of these blocks any manuscript number above; (a) should be remedied before submission.