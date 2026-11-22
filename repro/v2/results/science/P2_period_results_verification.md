# Orchestrator verification of P2 period-classifier outputs (preliminary, pending agent report)

Cross-checks performed against frozen inputs:
- My independent reproduction of the fixed-cut binned excess gave ΔSN = +0.0506 (young 0.3945 → old 0.4451) and median periods 10.05 d (young) vs 10.688 d (old), matching `results/phaseA_task2b.json` and the published narrative.
- The agent's `period_classifier_results.csv` at m=0 gives unadjusted dSN = 5.007 pp (z=2.16, p=0.0304) — consistent with the independent reproduction.
- The agent's adjusted AME at m=0 = 4.142 pp, β_logP = 2.383 — matches the frozen `phaseA_task2b.json` values (4.16 pp, 2.383) to rounding.
- Across m ∈ [−0.10, 0.0], adjusted AME = 4.45–5.25 pp, all p > 0.12: the period-dependent classifier does NOT restore significance. Adjusted values are labeled UPPER BOUND.
- Within-period test: only the 5.9–24.3 d bin shows a positive excess (6.50 pp, p=0.060); no bin is significant at 5% and the longest-period bin reverses sign (−0.80 pp).
- Permutation control (`permutation_control.json`): 2,000 host-level permutations, observed AME 4.142 pp, null mean −0.075 pp, null SD 2.55 pp, 90/2000 permutations ≥ observed → one-sided p = 0.045 (passes <5%), two-sided 0.104. This is a marginal result: report as "at the conventional threshold," not decisive.

Note: the regression AME (4.16 pp) and the binned ΔSN (5.0 pp) are different estimands; do not conflate.
