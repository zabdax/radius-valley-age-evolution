# Orchestrator independent check of effect-size conversions

Formulas: SN fraction f = 1/(1+SE/SN) for SE/SN ratios, or f = (SN/SE)/(1+SN/SE) for SN/SE ratios. logit(f)=ln(f/(1−f)). β per SD of the age variable = Δlogit / group-mean separation (1.60 SD for a median split).

| Source | published ratio | f_SN young → old | Δlogit | β/SD (1.60 sep) |
|---|---|---|---|---|
| Kamulali+26 | SE/SN 0.51 → 0.64 | 0.6623 → 0.6098 | −0.2267 | −0.1417 ≈ −0.14 (our primary) |
| Berger+20 | SE/SN 0.61 → 1.00 | 0.6211 → 0.5000 | −0.4943 | −0.3089 |
| Gaidos+24 | SN/SE 2.38 → 0.58 (ages 2.6 → 5.4 Gyr, 3.86 Gyr median split) | 0.7041 → 0.3671 | −1.4118 | −0.8824 on TRUE age |

Because our gate injects in v_tan units, the Gaidos effect on true age must also be shown attenuated by the age-proxy correlation: ×ρ=0.25 → −0.2206/SD; ×ρ=0.5 → −0.4412/SD. The cool-host ρ=0.25 is not significant (p=0.43), so the raw −0.88 is an intentionally optimistic upper bracket, not a preferred effect estimate.

Caveat: Berger's 1 Gyr split and Kamulali's 3 Gyr split are not median splits, so the 1.60-SD divisor is an approximation for those two; the Gaidos 3.86 Gyr median split is the exact case.

Computation verified with Python 3.12 math; cross-check against the science agent's `effect_size_derivation.json` when it arrives.
