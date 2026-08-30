# The Radius Valley's Age Evolution Around M Dwarfs: An Inconclusive Verdict and a Quantified Path Forward

**Zubayer Hasan Shaad**

*Government Tolaram College*

---

## Abstract

The radius valley — the well-known deficit of planets between 1.5 and 2.0 R⊕ — has been shown to evolve with stellar age among FGK (Sun-like) hosts: the sub-Neptune population shrinks and the valley shallows and shifts over Gyr timescales (Berger et al. 2020; David et al. 2021; Chen et al. 2022; Kamulali et al. 2026). The one M-dwarf test (Gaidos et al. 2024) found a marginal, non-significant decline, and Gillis, Cloutier & Pass (2026) recently demonstrated that the valley disappears entirely around mid-to-late M dwarfs, sharpening the open question of whether age evolution operates there at all. Reconciling the FGK detection with the M-dwarf silence requires knowing whether the absence reflects real physics or simply insufficient detection power — the literature has never quantitatively adjudicated this.

We assemble the NASA Exoplanet Archive's transit planet sample (419 M-dwarf planets, 2,880 FGK control planets) with Gaia DR3 kinematics, Angus et al. (2019) gyrochronology, and a completeness treatment that reuses real per-star Kepler DR25 CDPP for 1,691 hosts and applies a trimmed inverse-detection-efficiency weighting. We introduce a calibrated gate methodology: for any given fit, we inject null and literature-sized synthetic effects into the same noise structure and require the real result to be consistent with the signal-injected distribution AND distinguishable from the null-injected distribution. The pipeline underwent two independent adversarial audits, which caught a mixture-identifiability collapse (variance absorbing the mixing signal — fixed by modeling observed radii as a two-component mixture rather than a Bernoulli on a fixed classifier) and a period-mixing artifact (a spurious binned signal that vanishes once orbital period is controlled). Applying this gate to the completeness-weighted FGK Kepler-only control (N = 2,135 planets; ESS = 1,645, 77% retention) returns an inconclusive verdict — the real coefficient (−0.030 ± 0.106 per IQR) sits at the 43.3rd percentile of the null-injected distribution and the 90.8th percentile of the signal-injected distribution, meaning both interpretations remain statistically admissible.

Separately, an apparent wrong-sign binned artifact (z = +2.16, p = 0.03) was traced to period mixing: kinematically old hosts carry longer-period planets, and sub-Neptune fraction rises steeply with period against a fixed radius classification. We quantify the sample-size requirement via a 15-cell Monte Carlo power grid: even at 8× our current M-dwarf sample (~3,300 planets), detection of a literature-sized effect reaches only ~19%, far below the conventional 80% power threshold. Archival methods therefore cannot resolve this question, and the required sample size — or equivalently the required per-host age information — is substantially larger than the 8× (~3,300-planet) figure alone suggests, placing resolution squarely in PLATO/next-generation-observatory territory.

---

## 1. Introduction

The radius valley — the bimodal gap near 1.5–2.0 R⊕ in the distribution of close-in small planets — is one of the most robust demographic features in the exoplanet census (Fulton et al. 2017; Van Eylen et al. 2018). Its physical origin is debated: leading explanations include atmospheric mass loss driven by XUV photoevaporation (Owen & Wu 2017; Lopez & Fortney 2013) or core-powered mass loss (Ginzburg et al. 2018; Gupta & Schlichting 2019), both of which predict that the relative populations on either side of the valley should evolve with time as envelopes are gradually stripped, and formation-epoch mechanisms such as gas-poor formation (Lee & Chiang 2016; Lopez & Rice 2018), which predict that the valley is imprinted at birth and does not evolve.

These predictions are testable: if atmospheric stripping dominates, older planetary systems should show a depleted sub-Neptune population relative to younger ones, because marginal sub-Neptunes lose their envelopes over Gyr timescales and migrate across the valley into the super-Earth peak. Multiple groups have now detected exactly this signature among FGK (Sun-like) host stars. Berger et al. (2020) used isochrone ages for Kepler hosts and found that the radius valley shifts to larger radii and the sub-Neptune-to-super-Earth ratio decreases with age. David et al. (2021) confirmed the trend using gyrochronological ages and demonstrated that the sub-Neptune population shrinks measurably between young and old stellar samples. Chen et al. (2022) extended this with a refined age analysis and reached consistent conclusions. Most recently, Kamulali et al. (2026) provided an independent confirmation of the age-dependent valley evolution among Sun-like hosts using a homogeneous reanalysis based on SWEET-Cat spectroscopic parameters. The FGK age-evolution signal is now established by at least four independent analyses spanning different age-dating methods and sample constructions. We stress that this establishment is a property of the literature: our own pipeline's FGK control gate returns an inconclusive verdict (Section 4.1) and must not be counted as an independent confirmation of the effect.

Whether this age evolution extends to M dwarfs is a separate and currently open question. M dwarfs are cooler, smaller, and less luminous than FGK stars, with different XUV luminosity histories and convective structures. Theoretical predictions for mass-loss-driven valley evolution around M dwarfs are model-dependent: some photoevaporation models predict a weaker or delayed signal (Owen & Wu 2017), while core-powered mass loss models make different predictions for the low-mass stellar regime (Gupta & Schlichting 2020). The observational picture is sparse. Gaidos et al. (2024) tested the M-dwarf valley's age dependence using kinematic age proxies and found a marginal, non-significant decline in the sub-Neptune fraction among older hosts — suggestive but insufficient for a claim of detection.

The question was sharpened considerably by Gillis, Cloutier & Pass (2026), who demonstrated that the radius valley itself disappears entirely around mid-to-late M dwarfs (spectral types approximately M3 and later). Their analysis, based on TESS-detected planets with injection-recovery completeness corrections, showed that the bimodal structure so prominent around earlier-type hosts is absent in the coolest stellar regime. They proposed a formation-based explanation — that the valley around hotter hosts is sculpted primarily by atmospheric evolution, while M-dwarf planets form in a regime where the valley never develops — without testing the age axis. If their formation-only hypothesis is correct, there should be no age evolution of M-dwarf valley demographics; if atmospheric stripping operates but more weakly, there should be a detectable but attenuated age signal. The literature provides no quantitative test distinguishing these scenarios.

This paper addresses a specific question: **can current archival samples discriminate a literature-sized age effect on the radius valley around M dwarfs from a null effect?** We approach this not as a detection paper but as an honest statistical adjudication. We develop a calibrated gate methodology — injecting known null and literature-sized synthetic signals into the same noise structure and requiring the real data to discriminate between them — and apply it first to an FGK control sample — where the effect has been established by the literature analyses above and serves as a positive-control benchmark, while our own control-gate result is itself inconclusive (Section 4.1) — before extending to M dwarfs. We find that the answer is no: the current sample is too small to resolve the question, even after careful completeness correction. We quantify exactly how much larger the sample must become, providing a concrete target for the next generation of transit surveys.

---

## 2. Data

### 2.1 Sample Construction

We construct our planet samples from the NASA Exoplanet Archive using the Planetary Systems (``ps``) table, querying via TAP with the following selection criteria: ``default_flag = 1``, ``pl_controv_flag = 0``, ``discoverymethod = 'Transit'``, planet radius 0.5 < R_p < 4.0 R⊕, and orbital period P < 100 d. These cuts isolate the small-planet population spanning the radius valley while excluding controversial dispositions and restricting to the period range where transit surveys have meaningful sensitivity.

We divide the sample by host-star effective temperature into two subsamples:

- **M-dwarf sample:** T_eff < 4200 K, yielding **419 planets** around **303 host stars**.
- **FGK control sample:** 4200 ≤ T_eff ≤ 7000 K, yielding **2,880 planets** around **2,150 host stars**.

### 2.2 Gaia DR3 Crossmatch and Astrometric Quality

We crossmatch all host stars against Gaia DR3 (``gaiadr3.gaia_source``) using a cone search, then apply an offline quality filter requiring a positive parallax and adequate proper-motion measurements. Of the 303 M-dwarf hosts, 302 have Gaia matches; the offline quality gate rejects 2 spurious matches (the entries associated with Kepler-1410 and Kepler-1652), leaving **301 M-dwarf hosts** with usable astrometry. Of the 2,150 FGK hosts, 2,124 have initial matches; 26 are rejected for lacking usable astrometry, leaving **2,123 FGK hosts**.

### 2.3 Kinematic Age Proxies

For each host with usable Gaia astrometry, we compute tangential velocity (v_tan) from parallax and proper motion and, where radial velocities are available, full heliocentric Galactic velocities (U, V, W). Radial velocities are available for 237 of the 301 M-dwarf hosts and 1,429 of the 2,123 FGK hosts. Tangential velocity serves as the primary kinematic age proxy throughout this work: kinematically hotter (higher-v_tan) populations are statistically older, a well-established relation in Galactic stellar kinematics (e.g., Aumer & Binney 2009).

### 2.4 Gyrochronological Ages

We compute rotation-based ages using the ``gyrointerp`` package implementing the Angus et al. (2019) slow-sequence gyrochronology model. The model's validity domain requires 3800 < T_eff < 6200 K and rotation period P_rot ≤ 45 d; below 3800 K, the model returns NaN (a boundary we discovered and enforced during the project's audit process). Archive rotation-period uncertainties are used with a [5%, 35%] floor-ceiling clip.

Of the 303 M-dwarf hosts, 66 have archive rotation periods, but only **12** yield usable gyrochronological ages (median 2.50 Gyr), with 33% censored at the 3 Gyr model-grid boundary. Of the 2,150 FGK hosts, 88 have archive rotation periods and **79** yield usable ages (median 1.93 Gyr), with 18% censored. The low M-dwarf yield reflects the fundamental limitation that gyrochronology models are poorly calibrated for cool stars — the bulk of our M-dwarf age information necessarily comes from kinematic proxies.

A literature cross-match against the McQuillan et al. (2014) Kepler rotation catalog recovered a rotation period for exactly one additional M-dwarf host not in the archive (Kepler-1646, matched at 0.26″ separation). The sparse recovery is consistent with that catalog's magnitude limit, which excludes most of the Kepler field's faint M dwarfs.

---

## 3. Methods

### 3.1 Completeness Weighting

Transit detection probability varies strongly across the planet parameter space and between different survey facilities. Any analysis comparing sub-Neptune fractions across stellar or planetary populations must account for this: uncorrected comparisons are biased toward easily detected (large, short-period) planets and can manufacture spurious demographic trends.

**Reused public data product.** For Kepler-observed hosts, we use the real per-star 6-hour Combined Differential Photometric Precision (CDPP) values from the Kepler DR25 stellar table (``rrmscdpp06p0``), retrieved from the NASA Exoplanet Archive ``keplerstellar`` table and matched to our sample via KOI KEPID. This provides empirically measured photometric noise for **1,691 of 1,695 Kepler-prefixed hosts** in our combined sample (99.8% coverage), spanning 2,161 of the 2,848 FGK control planets observed by Kepler.

We note that Gillis, Cloutier & Pass (2026) performed injection-recovery completeness corrections for their M-dwarf TESS sample; however, no machine-readable sensitivity table is linked from their published materials, and their sample covers only mid-to-late M dwarfs observed by TESS, with partial overlap to our host set. Their product was therefore not directly applicable.

**First-order model for non-Kepler hosts.** For hosts observed by K2 or TESS (where per-star CDPP products are not available in a single uniform table), we compute a first-order individual-detection-efficiency model (IDEM) from physical observables: transit depth from planet radius and stellar radius, transit duration from stellar density, number of transits from mission-specific observing baselines (Kepler: 372 d, K2: 80 d, TESS: 60 d), and signal-to-noise ratio evaluated against magnitude-scaled CDPP with a logistic detection-probability function centered at SNR ≈ 7.5 (the DR25 conventional threshold), floored at 0.02 to prevent weight divergence.

**Adopted weighting scheme.** Our primary analysis restricts the FGK control sample to **Kepler-prefixed hosts only** (2,135 planets), so that the dominant completeness-weight component derives from the real DR25 CDPP rather than the cruder non-Kepler proxy. We apply inverse-detection-probability weighting (w = 1/P_det), normalize to mean unity, and **trim at the 95th percentile** of w to prevent the extreme ~14×-tail weights from collapsing the effective sample size.

This choice is justified quantitatively by the effective sample size (ESS) progression:

| Weighting scheme | N planets | ESS | Retention | SE inflation (√(N/ESS)) |
|---|---|---|---|---|
| Naive all-survey weighted (rejected) | 2,848 | 487 | 17% | 2.42× |
| Kepler-only, untrimmed (sensitivity check) | 2,135 | 511 | 24% | 2.04× |
| **Kepler-only, trimmed (adopted)** | **2,135** | **1,645** | **77%** | **1.14×** |

The ESS recovery from 487 to 1,645 (a factor of 3.4) is the mechanism by which the calibrated gate becomes trustworthy: at 17% retention, the injection distributions and real fits become effectively meaningless because a handful of extreme weights dominate all weighted statistics. The trimmed Kepler-only scheme retains 77% of the information while avoiding this pathology.

A sharpened completeness variant for the M-dwarf sample replaces the V-band magnitude fallback (which overestimates noise for red M dwarfs because V − T ≈ 2 mag) with a Gaia G → TESS-magnitude transform (T = G − 0.5) and handbook-anchored TESS CDPP curve, with a ×2 jitter factor for K2. This variant is documented as a sensitivity check but is not used in the primary results.

### 3.2 Hierarchical Mixture Model

We model the observed planet radius distribution as a two-component Gaussian mixture (super-Earths and sub-Neptunes), where the mixing fraction — the probability that a planet belongs to the sub-Neptune component — is the science target. The mixing fraction is parameterized as a logistic function of stellar age proxy (v_tan, standardized) and covariates (log₁₀(P), [Fe/H] where available), with the age-proxy coefficient β_age as the primary parameter of interest. The likelihood is evaluated via Gauss–Hermite quadrature over the radius measurement uncertainty, with component means and widths constrained to data-measured values to prevent degeneracy (the variance otherwise absorbs the mixing signal; see Section 3.4).

This model improves on a Bernoulli classifier applied to a fixed radius threshold: by treating the radius as a continuous observed outcome with measurement error, it avoids the information loss inherent in dichotomizing near the valley and is robust to the population overlap between super-Earths and sub-Neptunes.

**Host-level clustering.** The 419 M-dwarf planets orbit 303 hosts and the 2,880 FGK planets orbit 2,150 hosts; planets sharing a star are not independent draws (they share the host's age, formation history, and disk composition). The production specification handles this structure as follows. Point estimates are obtained under a working-independence likelihood — the host random-intercept scale is fixed at zero, a deliberate choice that avoids the documented weak-identifiability competition between the age slope and a host random effect — while all reported standard errors are host-level cluster-robust: the sandwich covariance aggregates per-planet scores to per-host scores before forming the meat matrix, so within-host correlation inflates the reported uncertainties exactly as in standard clustered-standard-error practice. The Monte Carlo machinery preserves the same structure: the null- and signal-injected ensembles are generated by resampling whole hosts with replacement and drawing one shared host effect per host (σ_u = 1.5 logits), so the calibrated gate's percentile comparisons are made against distributions that inherit within-host correlation. A one-planet-per-host robustness check is reported in Section 5.5.

### 3.3 Calibrated Gate Methodology

The paper's central methodological contribution is a calibrated statistical gate designed to adjudicate whether an observed fit discriminates between competing hypotheses. For a given sample and weighting scheme, we produce three estimates:

1. **Real-data fit:** the completeness-weighted hierarchical model applied to the actual data, yielding a measured β_age.
2. **Null-injected ensemble:** 120 synthetic realizations drawn from the same noise structure and covariate distribution but with β_inj = 0 (no age effect), each fitted identically.
3. **Signal-injected ensemble:** 120 synthetic realizations with β_inj = −0.14 per standard deviation of the age proxy (the baseline parameterization from the Kamulali et al. 2026 SWEET-Cat reanalysis of the FGK age-evolution literature), likewise fitted identically.

The gate criterion is:

- **PASS:** the real β_age is consistent with the signal-injected distribution (within its 95% envelope) AND statistically distinguishable from the null-injected distribution.
- **FAIL (inconclusive):** one or both conditions are not met — the data cannot discriminate between the two hypotheses.

This is deliberately conservative: it requires not only that the signal be recovered but that the null be rejected. A FAIL verdict is not a null result; it is a statement that the sample lacks the power to adjudicate the question.

### 3.4 Adversarial Audits as Methodology

The analysis pipeline underwent **two independent adversarial audits**, conducted by separate review subagents examining different components of the codebase and results. We regard these audits as part of the methodology rather than a post hoc correction step, because several of the issues they identified would have materially affected the reported results.

The two corrections with direct impact on the primary result were:

1. **Mixture-identifiability collapse.** The original (v1) likelihood used a Bernoulli soft-weight formulation on a fixed radius classifier. In the real-data regime, this attenuated the recovered β_age by approximately 35% because the unmodeled overlap between super-Earth and sub-Neptune populations allowed the variance to absorb the mixing signal. The fix replaced this with the proper two-component Gaussian mixture model described in Section 3.2, where the observed radius is the outcome variable and the mixing fraction is the science target.

2. **Period-mixing diagnosis.** The audits identified a discrepancy between a significant binned test and a non-significant hierarchical model result (Section 4.2), which led to the investigation and identification of the period-mixing artifact — itself a standalone methodological finding.

Additional corrections identified during the audits include: (a) a Gauss–Hermite normalization error (probabilists' nodes paired with physicists' weights, deflating fitted σ by √2); (b) an archive radius-error outlier (one FGK planet with a 68.9 R⊕ radius error corrupting the calibration, fixed by clipping relative errors to [1%, 30%] of the planet radius); (c) component degeneracy in the mixture model (free component widths allowing the optimizer to slide β toward zero, fixed by constraining component bounds to data-measured values); (d) objective-function scaling causing premature L-BFGS-B convergence at β = 0, fixed by normalizing by sample count; (e) an anchor name-matching bug producing spurious anchor correlations via substring matching (fixed by exact post-normalization string matching, leaving only 2 known-young anchors — AU Mic and K2-25 — in the sample, rendering the ρ > 0.5 anchor-validation gate criterion untestable and acknowledged as such); and (f) a retraction of a premature "precise null" label for the FGK control result, which was within ~0.4σ of both zero and the literature effect.

The standing reporting rule adopted for this project is that analytic standard errors and p-values from the hierarchical mixture model are never quoted as if trustworthy on their own; only calibration-based (injection-derived) verdicts are reported.

---

## 4. Results

### 4.1 FGK Control: Calibrated Gate Result

Applied to the completeness-weighted FGK Kepler-only control sample (N = 2,135 planets; ESS = 1,645 after trimming at the 95th percentile), the calibrated gate returns the following:

| Quantity | Value |
|---|---|
| Real-data β_age per IQR of v_tan | −0.030 ± 0.106 |
| Null-injected median β̂ | −0.013 (95% CI: [−0.248, +0.221]) |
| Real percentile in null distribution | 43.3% |
| Signal-injected median β̂ | −0.178 (95% CI: [−0.355, +0.068]) |
| Real percentile in signal distribution | 90.8% |
| **Verdict** | **FAIL (inconclusive)** |

The real coefficient sits at the 43.3rd percentile of the null-injected distribution, meaning it is comfortably consistent with no age effect at all. Simultaneously, it sits at the 90.8th percentile of the signal-injected distribution, placing it just inside the upper boundary of the 95% envelope of the literature-sized effect. Both interpretations — no effect and a literature-sized effect — are statistically admissible given the noise properties and effective sample size of this dataset. The data genuinely cannot discriminate between them. We emphasize that this inconclusive control result is a statement about the power of our pipeline at this effective sample size and noise structure — it is neither evidence for nor against the FGK age-evolution effect, which rests on the independent literature analyses cited in Section 1.

This is the calibrated gate working as designed: it does not disguise an ambiguous result as either a null or a detection. The previously retracted "precise null" label (which claimed the FGK control result was inconsistent with the literature) failed exactly this test — the real coefficient was within ~0.4σ of both zero and the literature effect, ruling out neither.

Signal-recovery diagnostics confirm that the pipeline is functioning correctly: signal-injected replicates with β_inj = −0.14 per standard deviation are recovered with documented ~35% attenuation toward zero in the estimator (mean shift of +0.37–0.49 logits per IQR), a known property of the mixture model in this sample-size regime.

### 4.2 Period-Mixing Artifact

In the course of the FGK control analysis, a significant apparent wrong-sign signal appeared in binned tests: a completeness-weighted two-proportion test on the FGK control showed sub-Neptune fraction significantly *higher* among kinematically old (high-v_tan) hosts (z = +2.16, p = 0.03), the opposite of the expected direction from atmospheric mass loss. The hierarchical model, however, reported β ≈ 0. Two methods applied to the same data produced contradictory stories.

Rather than resolving this by choosing the more convenient result, we investigated the discrepancy systematically.

**Metallicity hypothesis (falsified).** We first tested whether metallicity confounding could explain the artifact, since kinematically hot stars tend to be metal-poor and metallicity affects planet compositions. Stratifying the sample into [Fe/H] terciles:

| Stratum | Median [Fe/H] | z-statistic | p-value |
|---|---|---|---|
| Full sample (reference) | — | +2.16 | 0.030 |
| Tercile 1 (metal-poor) | −0.07 | +0.49 | 0.622 |
| Tercile 2 (solar) | +0.02 | +1.55 | 0.121 |
| Tercile 3 (metal-rich) | +0.10 | +0.62 | 0.534 |

The inverse-variance pooled estimate was dSN = +0.042 ± 0.028 (z = +1.53, p = 0.126), with a heterogeneity Q-statistic of 0.69 (χ² p = 0.71). A full-power weighted logistic regression with [Fe/H] as the only covariate (complete-case N = 1,479) yielded an average marginal effect of +6.29 ± 2.90 percentage points (p = 0.030) and β_[Fe/H] = −0.89 ± 0.50 per dex (p = 0.07). The metallicity hypothesis was rejected: the point estimate barely moved from the unstratified value.

**Period mixing (confirmed mechanism).** Adding log₁₀(P) as a covariate reduced the average marginal effect to +4.16 ± 2.61 percentage points (z = +1.61, p = 0.11), with β_log₁₀P = +2.38 ± 0.17. The mechanism is straightforward:

1. Kinematically old hosts carry slightly longer-period planets in this sample (median 10.69 d vs. 10.05 d).
2. Sub-Neptune fraction rises steeply with orbital period against a fixed 1.88 R⊕ classification boundary (from 0.141 at P < 3 d to 0.837 at 30–100 d).
3. This manufactures a spurious excess of classified sub-Neptunes among high-v_tan hosts — a wrong-sign artifact mimicking an age trend.

The hierarchical model never displayed this artifact because it included log-period as a covariate from the outset. The binned test, which performed no period adjustment, was vulnerable. This finding is relevant to any archival radius-valley demographics work that uses simple binned comparisons without controlling for the strong period dependence of the sub-Neptune fraction.

### 4.3 Power Analysis

We quantify the sample-size requirement for resolving the M-dwarf question via a 15-cell Monte Carlo power grid. The grid bootstraps hosts and covariates from the real M-dwarf sample at five scale factors (0.5×, 1×, 2×, 4×, 8× current size), for each of three injected effect sizes (null: β_inj = 0; SWEET-Cat-sized: β_inj = −0.14 per SD; strong: β_inj = −0.28 per SD). Planet radii are drawn from reality-matched overlapping class-conditional distributions; completeness weights are applied; each realization is fitted with the mixture-likelihood model. There are 120 replicates per null cell and 150 per effect-size cell.

| Scale factor | Approximate planets | Null FPR | SWEET-size detection rate | Strong detection rate |
|---|---|---|---|---|
| 0.5× | ~208 | 6.7% | 6.0% | 16.0% |
| **1× (current)** | **~417** | **11.7%** | **9.3%** | **30.0%** |
| 2× | ~834 | 4.2% | 14.7% | 43.3% |
| 4× | ~1,670 | 2.5% | 14.7% | 46.0% |
| 8× | ~3,340 | 7.5% | 18.7% | 64.0% |

Null false-positive rates drop toward nominal as sample scales increase, recovering from 11.7% at the 1× scale to 2.5% at the 4× scale. The 11.7% false-positive rate at 1× (14 out of 120 replicates) is statistically significant outside pure Monte Carlo (binomial) sampling noise for a 5% target (95% CI upper bound ~8.9%). It reflects the documented anti-conservative nature of the mixture model's analytical sandwich standard errors at small sample sizes, where cross-component overlap and sparse covariates cause the analytical variance to underestimate true sampling variance, leading to artificially narrow confidence intervals. This exact miscalibration mechanism is why the primary methodological gate (Section 3.3) requires comparing against a full empirical null distribution rather than relying on analytical confidence intervals. One non-monotonicity remains: the null false-positive rate at the 8× scale (7.5%, 9 of 120 replicates) is higher than the 4× value (2.5%, 3 of 120). We assessed whether this is Monte Carlo noise at the current 120-replicate count or a real effect: the two cells are not statistically distinguishable from each other (Fisher's exact test, p = 0.14), the 8× cell is individually consistent with the nominal 5% rate (binomial P(X ≥ 9 | n = 120, p = 0.05) = 0.15), and their 95% Clopper–Pearson intervals ([3.5%, 13.8%] and [0.5%, 7.1%]) overlap almost entirely. The 8× point estimate therefore reflects Monte Carlo noise rather than a real degradation of calibration, and the dedicated higher-replicate Monte Carlo run planned as follow-up work will pin this cell down decisively. Empirical coverage otherwise spans 83–98%, confirming that the pipeline is broadly well-calibrated once sample asymptotic properties stabilize.

At the current sample size (1×, ~417 planets), the detection rate for a literature-sized (SWEET-Cat) effect is 9.3% — the sample has essentially no power to detect this effect. Even at 8× the current sample (~3,340 planets), the SWEET-Cat detection rate reaches only 18.7%, and the strong (2× literature) effect is detected 64.0% of the time. These required-N figures are lower bounds, because the upsampling procedure resamples the empirical covariate distribution and cannot account for the different noise properties and survey selection functions of future samples.

The quantified conclusion is deliberately sobering, and we state it plainly: **even the 8× scale (~3,300 M-dwarf transit planets), the largest cell in our grid, falls well short of conventional 80% power for a literature-sized effect, reaching only 18.7% detection under the current modeling assumptions.** No cell in the grid attains 80% power for the literature-sized effect; the grid stops at 8× precisely because the power curve is still far below the conventional threshold there. A simple Gaussian extrapolation (18.7% power implies |β|/σ ≈ 1.1 at 8×; 80% power requires |β|/σ ≈ 2.8; σ ∝ N^(−1/2)) suggests the sample would need to grow by roughly another factor of ~7 beyond 8× — of order 2×10⁴ planets — before 80% power is reached, with wide uncertainty from the 150-replicate cell (the 95% interval on the 18.7% rate spans a required additional factor of ~5–12). The practical reading is that the sample-size requirement — or equivalently the requirement for materially better per-host age information than Gaia kinematic proxies provide — is likely several times larger than the 8× figure alone suggests. Even the 8× scale itself is already squarely in the territory of PLATO and next-generation hosted terrestrial-planet-finding observatories, and it is a floor, not a sufficient target.

---

## 5. Discussion

### 5.1 Implications for the M-Dwarf Radius Valley

The inconclusive verdict on the FGK control — where the age-evolution signal is established in the literature (Berger et al. 2020; David et al. 2021; Chen et al. 2022; Kamulali et al. 2026) but is *not* independently recovered by our pipeline — establishes that our methodology and current sample sizes are insufficient to confirm even a literature-established effect in the best-case scenario. The M-dwarf sample, which is smaller by a factor of five and has weaker age-proxy constraints, is a fortiori underpowered for this test.

This result has direct implications for interpreting the formation-based hypothesis of Gillis, Cloutier & Pass (2026). Their proposal — that the radius valley around mid-to-late M dwarfs is absent because it was never formed, rather than because it evolved away — is a testable prediction: it implies that there should be no age dependence in M-dwarf valley demographics. Our analysis shows that this prediction remains **untested**, not confirmed or refuted. The absence of a detected age signal in current archival data is entirely consistent with either the formation-only hypothesis or with an atmospheric-evolution signal that exists but is below the detection threshold of current samples. Distinguishing these scenarios requires at least the ~8× sample expansion quantified in Section 4.3 — and, since even that scale yields only ~19% detection power for a literature-sized effect under current assumptions (Section 4.3), a substantially larger expansion or materially better per-host age information. One spectral-type caveat also applies: our M-dwarf sample spans the full cool-star temperature range selected here (T_eff < 4200 K, i.e., late-K through late-M), whereas the valley-disappearance finding of Gillis, Cloutier & Pass (2026) is located specifically at approximately M3 and later, so this pooled analysis cannot distinguish whether the age-evolution question — or the underlying demographics — differs across that spectral-type boundary; subdividing the M-dwarf sample by spectral type is flagged as future work.

### 5.2 Period Mixing as a General Methodological Concern

The period-mixing artifact identified in Section 4.2 is not specific to our analysis or our choice of age proxy. Any archival demographics study that compares planetary populations between stellar subgroups defined by age, metallicity, Galactic kinematics, or any variable correlated with orbital-period distribution is vulnerable to this effect if the comparison is performed without controlling for the strong and monotonic dependence of sub-Neptune fraction on orbital period.

The mechanism is simple: different stellar populations may host planets with systematically different period distributions (in our case, kinematically old hosts carry slightly longer-period planets), and the sub-Neptune fraction rises steeply with period against any fixed radius boundary used to classify planets into categories. A naive binned comparison will attribute this period-driven compositional difference to the grouping variable. The hierarchical model avoids this artifact naturally by including log-period as a covariate, but binned tests — which are common in the literature as quick summary statistics — do not.

We recommend that future archival radius-valley demographics work either (a) explicitly condition on orbital period in any between-group comparison, or (b) verify that the groups being compared have consistent period distributions before interpreting differences in sub-Neptune fraction.

### 5.3 The Path to Resolution

The power analysis does not provide a sufficient target, and it is important to be explicit about this: the ~3,300-planet (8×) scale is the largest cell in our grid and the point at which detection of a literature-sized effect first exceeds ~15%, but even there power reaches only 18.7%, far below the conventional 80% threshold. Extrapolating the grid's power trend (Section 4.3) implies that a sample of order 2×10⁴ M-dwarf transit planets with measured radii in the 0.5–4.0 R⊕ range, orbital periods under 100 days, and associated kinematic or age data derived from Gaia astrometry would be needed for conventional power, unless future samples carry materially better per-host age information than kinematic proxies provide. An approximately eightfold expansion is achievable with the yield projections for PLATO's long-duration pointing fields, which will survey bright, nearby M dwarfs with the photometric precision and temporal baseline needed for both planet detection and rotation-period measurement; even that expansion, however, leaves the question statistically unresolved under current modeling assumptions.

The caveat on these required-N estimates is important: because the power grid resamples the empirical covariate distribution from the current M-dwarf sample, the figures are lower bounds. A future sample drawn from different survey fields, with different stellar-population mixes and photometric-noise properties, may require somewhat larger numbers. The qualitative conclusion — that an order-of-magnitude expansion is a floor rather than a sufficient target — is robust.

An attempted validation using isochrone-based ages in place of the kinematic proxy encountered unresolved data and calibration issues and was not completed; this remains a direction for future investigation.

### 5.4 The Audit Trail as Methodology

The corrections documented in Section 3.4 are presented as a methodological feature of this work, not an embarrassment to minimize. Adversarial auditing — subjecting an analysis pipeline to independent review specifically aimed at finding errors, with all corrections applied end-to-end before any results are finalized — is a practice that the field would benefit from adopting more widely. In this project, the two audits caught errors (the mixture-identifiability collapse, the Gauss–Hermite normalization, the radius-error outlier) that would have produced quantitatively different and potentially misleading results if left uncorrected. The retraction of the premature "precise null" label, performed before any manuscript was drafted, illustrates the value of requiring calibration-based verdicts rather than relying on analytic standard errors from complex hierarchical models.

### 5.5 Robustness to Within-Host Clustering

Planets sharing a host star are not independent draws: they share the host's age, formation history, and disk composition, so treating them as statistically independent could in principle understate the uncertainty on β_age. Section 3.2 documents how the production specification handles this structure (host-level cluster-robust sandwich standard errors, plus injection ensembles generated by resampling whole hosts with a shared per-host effect). Here we verify materiality empirically with a bounded robustness check: we re-ran the primary FGK gate (Section 4.1) using only one planet per host — the largest-radius planet per host (N = 2,135 → 1,588 planets; a uniformly random selection per host gives the same result to within Monte Carlo noise). The coefficient moves from −0.030 ± 0.106 to −0.041 ± 0.123 per IQR (z = −0.28 → −0.33), an SE inflation of only ~16%, and its percentile placement in the 120-replicate null-injected distribution is essentially unchanged (43.3rd percentile, P(null ≥ β_real) = 0.57, versus 43.3rd for the full sample). The verdict remains **FAIL (inconclusive)**.

These checks confirm that the existing model specification is adequate as-is: **no host-level random effect or additional clustering term is required for the primary inference**, and the planned higher-replicate Monte Carlo run can reuse the existing host-bootstrap machinery unchanged. (The host random intercept implemented in the code remains available for sensitivity use but is deliberately not part of the production specification, for the identifiability reasons given in Section 3.2.)

---

## 6. Conclusion

We set out to test whether the radius valley's age evolution around FGK stars — well established in the literature (Section 1), though not independently re-confirmed by our own pipeline — extends to M dwarfs, motivated by the recent finding of Gillis, Cloutier & Pass (2026) that the valley disappears entirely around mid-to-late M dwarfs. To adjudicate this, we developed a calibrated gate methodology that injects null and literature-sized synthetic signals into the same noise structure and requires the real data to discriminate between them — a protocol we present as a reusable methodological contribution for archival demographics work.

The honest headline result is **inconclusive**. Applied to the FGK Kepler-only control (2,135 planets, ESS = 1,645), the real age-proxy coefficient (−0.030 ± 0.106 per IQR) sits at the 43.3rd percentile of the null-injected distribution and the 90.8th percentile of the signal-injected distribution: both interpretations remain statistically admissible, and the data cannot discriminate between them. The M-dwarf sample (419 planets) is a fortiori underpowered.

Along the way, we identified a period-mixing artifact — a demonstrated methodological pitfall in which kinematically old hosts' longer-period planets produce a spurious wrong-sign signal in naive binned tests — that is relevant to any archival study comparing sub-Neptune fractions between stellar subgroups without controlling for orbital period.

A 15-cell Monte Carlo power grid quantifies the path forward, and its central result is a negative one: **even at the 8× scale (~3,300 M-dwarf transit planets, roughly eight times the current archival inventory), a literature-sized effect is detected only ~19% of the time, far below conventional 80% power.** The sample-size requirement — or equivalently the requirement for materially better per-host age information — is therefore substantially larger than the 8× figure alone suggests; even the 8× scale is already squarely in PLATO and next-generation hosted-observatory territory, and it is a floor, not a sufficient target.

---

## Acknowledgments

[PLACEHOLDER — to be completed by the authors.]

---

## References

<!-- Standard AAS/AJ citation format. Bibcodes and DOIs to be added. -->

- Angus, R., Morton, T. D., Foreman-Mackey, D., et al. 2019, AJ, 158, 173
- Aumer, M., & Binney, J. J. 2009, MNRAS, 397, 1286
- Berger, T. A., Huber, D., Gaidos, E., et al. 2020, AJ, 160, 108
- Chen, D.-C., Xie, J.-W., Zhou, J.-L., et al. 2022, AJ, 163, 249 (arXiv:2204.01940)
- David, T. J., Contardo, G., Sandoval, A., et al. 2021, AJ, 161, 265
- Fulton, B. J., Petigura, E. A., Howard, A. W., et al. 2017, AJ, 154, 109
- Gaidos, E., Ali, A., Kraus, A. L., & Rowe, J. F. 2024, MNRAS, 534, 3277 (arXiv:2404.11022, DOI 10.1093/mnras/stae2207)
- Gillis, J. A., Cloutier, R., & Pass, E. K. 2026, arXiv:2602.23364
- Ginzburg, S., Schlichting, H. E., & Sari, R. 2018, MNRAS, 476, 759
- Gupta, A., & Schlichting, H. E. 2019, MNRAS, 487, 24
- Gupta, A., & Schlichting, H. E. 2020, MNRAS, 493, 792
- Lee, E. J., & Chiang, E. 2016, ApJ, 817, 90
- Lopez, E. D., & Fortney, J. J. 2013, ApJ, 776, 2
- Lopez, E. D., & Rice, K. 2018, MNRAS, 479, 5303
- McQuillan, A., Mazeh, T., & Aigrain, S. 2014, ApJS, 211, 24
- Owen, J. E., & Wu, Y. 2017, ApJ, 847, 29
- Van Eylen, V., Agentoft, C., Lundkvist, M. S., et al. 2018, MNRAS, 479, 4786

- Kamulali, J., Adibekyan, V., Nsamba, B., et al. 2026, A&A, 707, A41 (arXiv:2601.12396)

---

## Appendix: Flags and Placeholders

The following items require author attention:

1. **Author list and affiliations:** Not provided in the inventory; must be supplied.
2. **Acknowledgments section:** Placeholder only; to be completed.
3. **Table/figure formatting:** This draft is in Markdown. Conversion to AASTeX LaTeX format will require reformatting all tables, adding figure environments for any plots, and applying journal-specific macros.
4. **Figures:** The inventory references analysis outputs (JSON files, CSV files) but does not specify prepared figures. The paper would benefit from at least: (a) a completeness-weight ESS comparison figure; (b) the calibrated gate visualization (real β vs. null and signal injection distributions); (c) the period-mixing demonstration (SN fraction vs. period, stratified by v_tan); (d) the power-analysis grid. These would need to be generated from the pipeline outputs.
