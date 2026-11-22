# Analytic cross-check of the power-vs-effect curve (orchestrator)

The gate PASS rule has two parts: (1) real coefficient distinguishable from the null arm, P(null ≥ real) ≤ 0.05; (2) real inside the signal arm's 95% envelope. Part (1) depends ONLY on the real coefficient and the null distribution, both of which are independent of β_inj. For the cool sample the real coefficient sits at the 55.3rd percentile of the null, so P(null ≥ real) ≈ 0.447 > 0.05 at EVERY injected effect size. Therefore the verdict is FAIL at all points of the curve, by construction. This is a structural result, not a Monte Carlo artifact.

Part (2) does depend on β_inj through the signal-arm median. From the frozen production gate (`results/highrep/task2_calibration_kepler_only_trim0.95_rep1000.json`): cool null median −0.0513 (SD 0.6720), signal median −0.1812 at β_inj = −0.14 → shift −0.1299, response slope ≈ 0.928 per unit β. Linear response implies:

| β_inj (per SD v_tan) | approx signal median | approx arm separation (null SD) | verdict |
|---|---|---|---|
| 0 | −0.051 | 0.00 | FAIL |
| −0.1419 | −0.183 | −0.20 | FAIL |
| −0.2206 | −0.256 | −0.30 | FAIL |
| −0.3089 | −0.338 | −0.43 | FAIL |
| −0.4412 | −0.461 | −0.61 | FAIL |
| −0.8824 | −0.870 | −1.22 | FAIL |

Interpretation: at a Gaidos-sized injection the two injected ensembles become well separated (~1.2 null SD), so the gate WOULD be sensitive to such an effect, but the real data show no age signal relative to the null (55th percentile). The honest statement is that the archive is power-limited for a Kamulali-sized effect, while a Gaidos-sized effect would be detectable in principle but is not present in the real coefficient. This does not by itself prove the ordinal age axis is binding; it shows the observed coefficient is null-like and the gate's sensitivity improves with effect size. Cross-check against the Monte Carlo curve when it completes.
