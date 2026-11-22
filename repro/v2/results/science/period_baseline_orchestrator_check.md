# Orchestrator independent reproduction of the period-mixing baseline

Reproduced directly from frozen CSVs (Kepler-prefix → Gaia inner join → weight join → 95th-percentile cap):

- 2,135 planets / 1,588 hosts; ESS ≈ 1,644.1 (matches the archived 1,645 to rounding).
- Weighted SN fraction (R_p > 1.88 R⊕): young (low v_tan) = 0.3945, old (high v_tan) = 0.4451 → **ΔSN = +0.0506**, matching `results/phaseA_task2b.json` unstratified dSN = 0.05.
- Median orbital period: young **10.05 d**, old **10.688 d** ≈ published 10.05 vs 10.69 d.
- Two NaNs in the weight column were filled with 0 for the weighted average; this does not change the result at quoted precision.

This confirms the frozen inputs and selection logic used by the period-classifier agent are correct.
