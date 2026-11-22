# Period mixing: classifier, within-period test, permutation control

## Adopted control (reproduced exactly)
- Kepler-prefix, Gaia-inner-joined, 95th-percentile weight-trimmed FGK sample.
- 2135 planets / 1588 hosts; ESS = 1645; weight cap = 1.119.
- Age proxy v_tan, split at the sample median (33.51 km/s); SN if R_p > R_gap(P).

## Reproduction of first-party facts
- Fixed 1.88 Re weighted two-proportion SN excess (old-young): +5.01 pp, se 2.31 pp, z=+2.16, p=0.0304 (stored +5 pp, z=2.16, p=0.0304). The paper's '+6.29 pp with z=2.16' mixes the logistic AME point estimate with the binned z-statistic.
- Weighted logistic (finite-metallicity complete case, N=1479): AME(Fe/H only) = +6.27 pp [stored 6.29], beta_FeH=-0.894 [stored -0.894]; adding log10(P): AME=+4.14 pp [stored 4.16], beta_logP=+2.383 [stored 2.383]. Point estimates reproduce exactly; our delta-method AME SEs are larger than the stored +-2.90/+-2.61 pp (those were model/cluster SEs), so the stored z/p are quoted for the significance claims.
- Median periods: old 10.69 d vs young 10.05 d (stored 10.69 / 10.05).
- Metallicity terciles (within-tercile v_tan median split): pooled dSN=+4.24 pp, z=+1.53, p=0.126; Q=0.69, p_het=0.709 (stored z=1.50, p=0.134, Q=0.69). All null.

## Task 1: period-dependent classifier
- R_gap(P) = R_10 (P/10d)^m with R_10 = 1.88 Re at 10 d (the fixed-cut value adopted in the paper's baseline); SN if R_p > R_gap.
- m bracket -0.10, -0.075, -0.05, -0.025, 0.0; unadjusted = HT weighted two-proportion at host v_tan median; adjusted = weighted logistic AME with Fe/H and log10(P) (labelled UPPER BOUND).

| m | R10 | R_gap_1d | R_gap_100d | dSN_unadj_pp | unadj_se_pp | unadj_z | unadj_p | AME_adj_pp | adj_se_pp | adj_z | adj_p | beta_logP | adj_label |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| -0.100 | 1.880 | 2.367 | 1.493 | 6.457 | 2.352 | 2.745 | 0.006 | 4.452 | 3.303 | 1.348 | 0.178 | 3.270 | UPPER BOUND (residual period controlled) |
| -0.075 | 1.880 | 2.234 | 1.582 | 6.528 | 2.350 | 2.778 | 0.005 | 4.833 | 3.348 | 1.444 | 0.149 | 3.099 | UPPER BOUND (residual period controlled) |
| -0.050 | 1.880 | 2.109 | 1.676 | 6.539 | 2.340 | 2.794 | 0.005 | 4.727 | 3.395 | 1.392 | 0.164 | 2.920 | UPPER BOUND (residual period controlled) |
| -0.025 | 1.880 | 1.991 | 1.775 | 6.678 | 2.334 | 2.861 | 0.004 | 5.254 | 3.467 | 1.515 | 0.130 | 2.655 | UPPER BOUND (residual period controlled) |
| 0.000 | 1.880 | 1.880 | 1.880 | 5.007 | 2.313 | 2.164 | 0.030 | 4.142 | 3.535 | 1.172 | 0.241 | 2.383 | UPPER BOUND (residual period controlled) |

- Data-driven m fit: both a flat-SN-fraction criterion and a valley-emptiness criterion hit the search-grid boundary, so m is NOT statistically identifiable in this sample; the analysis therefore uses the documented bracket (returned -0.300 at the boundary). Literature: FGK/photoevaporation valley slope dlogR/dlogP ~ -0.09 to -0.12 (Van Eylen+18, Fulton+18); cool hosts are flatter or opposite: Gaidos+24 -0.03, Luque & Pale 22 +0.02 +/- 0.04, Cloutier & Menou 20 opposite insolation-sign dependence.

## Task 2: within-period age test
- 4 equal-width log10(P) bins; within each, host split at the global v_tan median; fixed 1.88 Re classifier.

| bin | logP_lo | logP_hi | P_lo_d | P_hi_d | N | N_hosts | medP_d | SN_old | SN_young | dSN_pp | se_pp | z | p | sig_5pct |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | -0.450 | 0.162 | 0.355 | 1.453 | 93 | 93 | 0.935 | 0.023 | 0.013 | 0.961 | 2.604 | 0.369 | 0.712 | False |
| 2 | 0.162 | 0.774 | 1.453 | 5.948 | 581 | 525 | 3.789 | 0.227 | 0.216 | 1.148 | 3.307 | 0.347 | 0.728 | False |
| 3 | 0.774 | 1.386 | 5.948 | 24.344 | 1000 | 853 | 11.897 | 0.478 | 0.413 | 6.502 | 3.463 | 1.878 | 0.060 | False |
| 4 | 1.386 | 1.998 | 24.344 | 99.640 | 461 | 430 | 40.716 | 0.705 | 0.713 | -0.796 | 5.444 | -0.146 | 0.884 | False |

- IVW pooled dSN=+2.14 pp (z=+1.28); sign positive in 3/4 bins. Significant bins at p<0.05: 0.

## Task 3: age-permutation negative control
- Permuted host-level v_tan among hosts (planets move with their host) 2000 times, recomputing the period-adjusted AME.
- Observed +4.14 pp; permutation null mean -0.07 +/- 2.55 pp.
- Fraction >= observed: 0.0450 (90/2000); 5% gate: PASS.
- The rate is only marginally below 5% (Monte Carlo SE ~0.5 pp at 2000 permutations), so this is weak/borderline evidence, not a robust detection; the residual +4.14 pp is also null under a two-sided test.
- Two-sided fraction |null| >= |observed|: 0.1040.

## Verdict
Both. The +5.0 pp binned excess at a FIXED 1.88 Re cut is reproducible (z=2.16, p=0.030), but it is substantially a period-mixing artifact: kinematically old hosts carry longer-period planets (10.69 vs 10.05 d) and the SN fraction rises steeply with period (beta_logP ~ +2.38), so a fixed cut manufactures an apparent age trend. Conditioning on logP reduces the AME from +6.27 to +4.14 pp. The residual is NOT zero and is in the tail of the host-label permutation null, so the honest statement is that period mixing explains part of the signal and the remainder is at most weak evidence for intrinsic valley evolution; the paper must not claim a clean detection.

## Files
- period_classifier_results.csv, within_period_test.csv, permutation_control.json, period_mixing_summary.md
- Script: scripts/period_mixing_analysis.py
