# PROJECT STATUS REPORT
## Does the exoplanet radius valley evolve with stellar age around M dwarfs as it does around FGK stars?

**Standalone status report — prepared 2026-08-26 during a third, adversarial verification audit.**
Every factual claim is tagged **[VERIFIED]** (re-derived from files/code/reruns during this audit) or **[PROVISIONAL]** (taken from project records not independently rechecked). Audit scratch code lives in `results/audit3_tmp/`.

---

## 1. Research question and why the gap matters

Small transiting planets split into two populations — super-Earths (~1.2–1.3 R⊕) and sub-Neptunes (~2.1–2.4 R⊕) — separated by the *radius valley*. Atmospheric-mass-loss theory (core-powered mass loss; photoevaporation) predicts that the valley's demographic signature should sharpen **with stellar age**, as planets lose their envelopes and migrate across the valley.

- **FGK baseline (quantified).** Kamulali et al. 2026 (arXiv:2601.12396, "Revisiting the exoplanet radius valley with host stars from SWEET-Cat") used homogeneous SWEET-Cat parameters for 1,221 main-sequence stars hosting 1,405 planets and found that as stellar age increases the super-Earth/sub-Neptune (SE/SN) ratio climbs from **0.51 (+0.11/−0.08) below ~3 Gyr to 0.64 (+0.11/−0.11) above it**, while the valley becomes shallower and shifts to larger radii. [VERIFIED against abstract]
- **M-dwarf test (one, underpowered).** Gaidos et al. 2024 (arXiv:2404.11022, "The Radius Distribution of M dwarf-hosted Planets and its Evolution") refined radii for 117 KOIs orbiting 74 single late-K/M dwarfs with rotation-based ages and found the sub-Neptune/super-Earth ratio **declines between the younger and older halves** (split at median age 3.86 Gyr), but only marginally significant (**p = 0.08**). They also found the valley's period dependence markedly weaker than for Sun-like stars. [VERIFIED against abstract]

**Gap:** whether the quantified FGK age evolution exists around M dwarfs is undetermined; the only existing test is consistent both with a real decline and with pure noise. This project builds an independent M-dwarf test plus a power analysis against the quantified FGK baseline. No scientific claim about age evolution is made yet (see §4, §5).

---

## 2. Data assets built (all counts re-derived from `data/*.csv` this audit) [VERIFIED]

| Asset | Content | Count |
|---|---|---|
| `planet_sample_M.csv` | NEA `ps`, default_flag=1, controv=0, transit, 0.5<R<4 R⊕, P<100 d, host Teff<4200 K | **419 planets / 303 hosts** |
| `planet_sample_FGK.csv` | same cuts, Teff 4200–7000 K (control) | **2880 planets / 2150 hosts** |
| `gaia_hosts_M.csv` | Gaia DR3 crossmatch + quality gate (parallax>0, finite PM); rejects logged (`Kepler-1410`, `Kepler-1652`: null parallax/PM) | **301 hosts** (RV: 237 = 79%) |
| `gaia_hosts_FGK.csv` | same; 26 rejects logged (incl. negative-parallax Kepler-1248) | **2123 hosts** (RV: 1429 = 67%) |
| `hosts_kinematics_M/FGK.csv` | v_tan = 4.74·μ·d for all hosts; W (heliocentric Galactic) where RV | 301 / 2123 hosts |
| `host_ages_M.csv` | Angus+19 gyrochronology (gyrointerp, Teff window 3800–6200 K hard-enforced; archive Prot errors clipped to [5%,35%]) | 66 rotators → **12 usable** (median 2.50 Gyr; 33% censored at 3 Gyr ceiling) |
| `host_ages_FGK.csv` | same | 88 → **79 usable** (median 1.93 Gyr; 18% censored) |

Age proxies: **v_tan** (all hosts, mission-independent) and **|W|** (RV subset) used as *ordinal* rankings only; gyrochronology covers only the warm tail (M-dwarf coverage 12/303 ≈ 4%). Literature anchor coverage for the kinematic axis is **n = 2 in-sample** (AU Mic: 22 Myr, 21 km/s; K2-25: 650 Myr, 26 km/s) — directionally consistent but statistically untestable. [VERIFIED]

Known accepted gaps [PROVISIONAL, prior audit]: no epoch-propagation of coordinates before cone search; one FGK host (Kepler-515) lost to retry exhaustion; ~60 planets with NULL `st_teff` excluded by design.

---

## 3. Methods implemented and their validation status

### 3.1 Binned validation suite (`04_validate.py`) — Fisher tests on valley-crossing (R > floor) fractions, young vs old halves. Floor: M 1.85 R⊕, FGK 1.88 R⊕ [floors PROVISIONAL from prior audit; rerun of 04 done this audit].
**Status: reproduces ROADMAP session table exactly** [VERIFIED — rerun this audit]:
- A: M × gyro age — N=16 pl/12 hosts; SN frac 0.875 → 0.250; Fisher p=0.041 (flagged SUGGESTIVE ONLY in-code).
- B: FGK × gyro age — N=128/79; 0.688 → 0.656; p=0.851 (NULL, direction-only).
- C: M × v_tan — N=417/301; 0.411 → 0.404; p=0.921 (NULL).
- C2: M × |W| — N=317/237; 0.412 → 0.420; p=0.910 (NULL; direction does not match Gaidos).

### 3.2 Wrong-sign diagnostic (`08_headline_fgk_vtan.py`) — reproduces every ROADMAP headline number exactly [VERIFIED — rerun this audit]: FGK × v_tan all-hosts SN frac young 0.561 → old 0.604, **p=0.023 but WRONG SIGN**; survives |Fe/H|<0.05 cut (p=0.022); killed by distance <200 pc (0.687 vs 0.701, p=0.77); localized entirely in Kepler-field hosts (N=2135, p=0.024; non-Kepler subsets null). Interpretation: raw planet counts vs kinematic heat conflate Kepler-field survey geometry/completeness with any true age signal.

### 3.3 Hierarchical Bayesian logistic model (`09_hierarchical_model.py`)
P(SN | age, covariates) with soft radius classification (weight w = Φ((r_obs − floor)/σ_r)), host random intercepts (24-node Gauss-Hermite), covariates log P, computed log S, mission dummies (Kepler/K2 vs other; ±Fe/H variant); MAP + Laplace CIs; priors α~N(0,5²), others N(0,2.5²), log σ_u~N(−1.5,1²). Published in JSON for six configurations; **all numbers in the ROADMAP MODEL RESULTS table match the JSONs exactly** [VERIFIED].

**However — three defects found by this audit materially change interpretation (§5, items 1–3).** Corrected refits were run in scratch (`audit3_tmp/a2_refit_variants.py`):

| Fit | As published (per IQR) | Insolation-corrected refit |
|---|---|---|
| M × v_tan | +0.30 ± 0.29 (z=+1.06) | **+0.37 ± 0.25 (z=+1.46)** |
| M × \|W\| | −0.15 ± 0.27 (z=−0.53) | **−0.02 ± 0.25 (z=−0.10)** |
| FGK × v_tan | +0.12 ± 0.12 (z=+0.97) | +0.10 ± 0.12 (z=+0.86) |

The qualitative conclusion survives (no significant effect either way), but two published narratives do not: (i) the "sign inconsistency between proxies reinforces that nothing is detectable" argument was largely an artifact of the broken insolation column — after the fix both M fits are consistent with each other (≈0 to positive); (ii) insolation was **not** actually controlled despite the docstring/ROADMAP claim.

### 3.4 Monte Carlo power analysis (`11_power_analysis.py`) — bootstrap-upsample synthetic M samples at 0.5×–8× the real 417 planets; inject β_age per SD of v_tan (null / −0.14 "SWEET" = the Kamulali-sized effect / −0.28 "strong"); draw radii from overlapping class-conditionals N(1.30,0.22)/N(2.35,0.45) + per-planet errors; run the full 09 estimation; detection = 95% CI excludes 0 with injected sign. All 15 cells completed in `results/power_analysis.json` [VERIFIED].

---

## 4. All results to date, with honest significance labels

1. **M dwarfs × gyrochronology (Test A):** SN fraction 0.875 (young) vs 0.250 (old), nominal p=0.041 — **SUGGESTIVE ONLY**: N=16 planets, uncorrected for censoring/multiplicity/trials; would not survive any multiplicity correction. [VERIFIED]
2. **FGK × gyrochronology control (Test B):** null (p=0.85, direction-only agreement with published decline). [VERIFIED]
3. **M dwarfs × kinematic ranking (Tests C, C2):** null under v_tan (p=0.92) and |W| (p=0.91). [VERIFIED]
4. **Powered FGK positive control produces a wrong-sign artifact:** significant (p=0.023) *increase* of SN fraction with v_tan, driven by Kepler-field geometry (kills on local-volume cut; mission-localized). Treated as a systematic demonstration, never as science. [VERIFIED — rerun]
5. **Hierarchical model: no detectable age evolution around M dwarfs under any proxy/cut** (all |z| ≤ 1.46; published table max |z|=1.37). Sanity checks pass: Fe/H coefficient +0.83 ± 0.24 (metal-rich → more sub-Neptunes, matches literature SE/SN-vs-Fe/H); strong period dependence recovered. [VERIFIED vs JSONs]
6. **Power analysis (the project's main quantitative result so far):** at the actual sample size (417 planets), a FGK-sized effect (−0.14/SD) is detected only **15.3%** of the time, and the null false-positive rate at this scale is **9.2%** (above nominal 5%). The 2×-literature effect needs ≥4× (~1670 planets) for coin-flip power (50.7%) and reaches 67.3% at 8× (~3340 planets). The SWEET-size effect stays below 25% even at 8×. **Conclusion: current-generation M-dwarf samples cannot distinguish presence from absence of FGK-like valley age evolution; the Gaidos+24 non-detection was essentially guaranteed by sample size.** [VERIFIED vs power_analysis.json]

---

## 5. Known limitations and threats to validity (blunt; new audit findings first)

**Critical (found this audit; must be fixed before any external release):**

1. **Insolation units bug** (`09.build_design` L66, copied into `11.synth`): luminosity coded as `R*² · Teff/5772⁴` instead of `R*² · (Teff/5772)⁴`. Computed S values are wrong by ~10 orders of magnitude and carry a spurious Teff⁻³ tilt. [VERIFIED numerically]
2. **Consequence via the outlier guard:** the guard `vz[(v<-9)|(v>9)]=0` silently zeroes the z-scored logS column wherever raw buggy-logS < −9 — which is **79.6% of M rows, 78.2% of |W| rows, 59.8% of FGK rows** (logP: 0 rows affected everywhere). The fitted γ_logS is therefore not an insolation coefficient but a contrast of "very hot / very short-period" survivors (true S ≳ 127 S⊕ for M, ≳ 425 S⊕ for FGK) against everyone else imputed at the mean — collinear with period and mission geometry. Any text claiming insolation control is currently false. With the correct formula, **zero** rows would be guarded. [VERIFIED census]
3. **Gauss-Hermite quadrature mixes conventions:** `hermegauss` nodes (probabilists', weight e^(−x²/2)) are scaled with the physicists' rule u=√2σx and renormalized by √π instead of √(2π). Verified numerically: the grid integrates N(0,**2σ_u²**), not N(0,σ_u²), and inflates each host likelihood by √2 (constant, harmless for MAP). Net effect: the true host-level mixing SD is **√2 × the reported σ_u** (published "σ_u≈1.4–1.6" ⇒ effective ≈ 2.0–2.3 logits). Refits show point estimates of β_age are barely moved by fixing this (Δ≤0.01), but all quoted σ_u values are systematically deflated. [VERIFIED]
4. **Selftest does not validate the regime the model is used in.** Published SELFTEST: 6 seeds, one β, host-effect SD 0.35 — reproduced (6/6 coverage; est −1.11…−1.44 vs truth −1.39). Extended here to 24 seeds × 3 regimes: in the paper regime bias is +0.16 logits (~12% attenuation, coverage 21/24 = 87.5%, already below nominal); at the **real-data host-effect scale (SD 1.5)** attenuation grows to **+0.51 logits (~37%) with 0/24 CI coverage**. Laplace CIs are strongly anti-conservative precisely where the model is applied. Consistently, the power-sim medians show a persistent ~10–20% attenuation of β̂ that does not shrink even at 8× sample size. Detection rates are unaffected in validity (simulation and estimator share the pipeline), but quoted effect sizes are biased low. [VERIFIED]

**Material caveats (quantified this audit):**

5. **Soft-classification coherence.** w = Φ((r_obs−floor)/σ_r) treats label uncertainty as if class-conditional densities crossed at the floor with width equal to the *measurement* error. It ignores intrinsic population widths (SE ≈0.22, SN ≈0.45) and base rates; since σ_r << intrinsic widths (median σ_r = 0.10 M / 0.20 FGK), w is near-step-function overconfident for boundary planets — effectively a hard cut with slight smearing. Limits behave sensibly (σ_r→0: hard cut; σ_r→∞: w→0.5, row auto-downweighted), but boundary probabilities are miscalibrated. Quantified ambiguous band (0.05<w<0.95): **20.6% of M, 36.2% of FGK planets**. Additionally the model uses a *static* floor although the valley demonstrably slopes with period/insolation (both cited papers), which can distort γ_logP and attenuate β_age. [VERIFIED census; interpretation]
6. **Power-analysis generative fidelity is imperfect in the optimistic direction.** Synthetic radii match the real M-sample marginal acceptably (KS D=0.056, p=0.15) because calibration forces the crossing fraction, but the components are misplaced: KDE modes 1.32/2.43 vs real 1.23/2.08; boundary-window (1.55–2.15 R⊕) mass 21.1% vs real 28.1%. Cross-valley confusion probability under typical errors is ~0.15 synthetic vs ~0.32 reality-matched — i.e., **real data have roughly twice the label ambiguity, so published detection rates are optimistic by an unquantified but material amount.** [VERIFIED]
7. **`calibrate_alpha()` omission:** solving the intercept with only β·logP ignores mission dummies (fractions Kepler 30.2%, K2 22.8% ⇒ −0.278 logits) and logS. Realized pre-radius SN probability ≈0.385 vs 0.453 target; Bernoulli information −4.4%. Effect on power conclusions: minor (<few %), same direction as item 6. [VERIFIED]
8. **Null FPR at the operating point is ~9%, not 5%.** At 1× scale the false-positive rate is 9.2% (binomial 95% CI 4–14%), converging to ~4–5% only at ≥2×. Small-sample Laplace CIs are anti-conservative — relevant when quoting z-scores from the real 417-planet fits. Coverage across cells 88.7–95.8% (mean 93.2%). [VERIFIED]
9. **Non-reproducible seeds:** `11` derives RNG seeds from salted `hash(effect_name)` — verified to differ across processes; `power_analysis.json` cannot be regenerated bit-exactly. [VERIFIED]
10. **No completeness/occurrence treatment anywhere.** The 08 artifact demonstrates distance-dependent selection directly; until survey completeness (or IDEM-style weighting) enters the model, *any* sign of the age coefficient — including the residual positive-v_tan tendency that survives mission dummies — is uninterpretable. This is the single largest scientific threat. [VERIFIED as absent from code]
11. **Kinematic proxy validation rests on n=2 anchors;** gyro ages cover 4% of M hosts with 33% right-censoring; the stale pre-fix file `results/kinematic_validation_FGK.txt` still contains voided substring-match numbers (e.g., K2-25 at 100 km/s contradicting the corrected 26.3 km/s) and should be deleted or annotated. [VERIFIED]

**Accepted/carried-over limitations [PROVISIONAL, prior audits]:** no epoch-propagation before cone search; Kepler-515 absent from Gaia table; ~60 null-Teff planets excluded; bootstrap upscaling resamples empirical covariates (required-N figures are lower bounds — correctly documented in 11).

---

## 6. What remains before submission-worthiness

**Must-fix (bugs), in order:**
1. Fix the luminosity formula `(Teff/5772)**4`; drop or recalibrate the outlier guard (with correct units it never fires); refit all six configurations and rewrite the hierarchical-results narrative (the |W|-vs-v_tan sign-inconsistency argument disappears).
2. Fix the GH grid (`u=σx, w=W/√(2π)`); report σ_u at true scale (expect ~2.0–2.3).
3. Rebuild the SELFTEST: ≥100 seeds, β ∈ {0, ±0.14, ±0.8}, host-effect SD matched to data (≈1.5–2.2), mission mix mirrored; report bias curves and coverage honestly.
4. Regenerate the power analysis with the fixed pipeline, reality-matched component means (SN mode ~2.08), reproducible seeded RNG, and calibration including mission dummies; expect detect rates somewhat *lower* than the current table.
5. Delete/annotate stale artifacts (`kinematic_validation_FGK.txt`; ROADMAP 8× "(pending)" row → 4.2%/24.0%/67.3%).

**Must-do (science):**
6. Completeness/selection modeling (at minimum IDEM-style weights or stratified local-volume analyses; ideally forward-modeled occurrence), motivated explicitly by the 08 artifact, before interpreting any age coefficient.
7. Increase M-host rotation coverage (script 10: light-curve Prot for cool hosts) — gyro ages are the only route to a *continuous* age axis; 12 hosts cannot anchor Test A.
8. Expand kinematic-anchor set beyond n=2 (benchmark clusters/moving groups crossing the sample) to validate the v_tan/|W| ordinal assumption.
9. Period-dependent valley floor inside the classification weight (the valley is not vertical in R–P space in either cited paper).

**Bottom line.** The pipeline's central negative finding — no detectable M-dwarf age evolution, and a demonstration that FGK-sized signals are undetectable at current M-dwarf sample sizes (power ≈15% at n=417) — is robust to every defect found and is the publishable core. The hierarchical model's machinery needs the four fixes above before its numbers appear anywhere public; its qualitative conclusion is unchanged by corrected refits run during this audit.

---

## Appendix: verification trail for this report

- Scripts rerun clean this audit: `04_validate.py` (table matches ROADMAP exactly), `08_headline_fgk_vtan.py` (all headline numbers exact), `09_hierarchical_model.py SELFTEST` (6/6 as claimed).
- Corrected-pipeline refits, extended selftests, guard/ambiguity censuses, KS fidelity tests, quadrature unit checks: `results/audit3_tmp/a1_guard_insolation.py`, `a2_refit_variants.py`, `a3_fidelity_power.py`, `a4_selftest_hessian.py`.
- External citations checked against arXiv abstract pages on 2026-08-26 (2601.12396, 2404.11022).


---

## ADDENDUM (2026-08-26, post-audit remediation — supersedes stale numbers above where they conflict)

All five critical findings from this audit were fixed and the affected results regenerated:

1. **Insolation formula fixed** (`(T/5772)^4` precedence) in 09 and 11; z-guard replaced
   by winsorization. With correct units the guard fires on zero rows.
2. **Gauss-Hermite normalization corrected** (probabilists' nodes with /sqrt(2pi)
   weights). Previously integrated N(0, 2 sigma^2), deflating fitted sigma_u by ~sqrt2.
3. **Root cause of selftest attenuation identified**: unmodeled super-Earth/sub-Neptune
   population overlap (intrinsic radius scatter >> measurement errors) caused classic
   misclassification attenuation (~35%) in the Bernoulli soft-weight likelihood.
4. **Model v2 implemented**: observed radius modeled directly as a two-component
   Gaussian mixture with covariate-dependent mixing fraction (the science target);
   sigma_u fixed at 0 with cluster-robust sandwich SEs over hosts. Residual estimator
   attenuation characterized honestly: ~+0.37..0.49 logits/IQR, CI coverage 25-30% in
   the real regime -> CONSERVATIVE for null conclusions; calibration printed by every
   selftest run and must be applied to any future positive detection.
5. **Power analysis regenerated (FINAL)** with reality-matched class conditionals,
   stable seeds, mission-offset-corrected calibration, and generative host
   heterogeneity sigma_u=1.5. Key numbers: detection of the Kamulali+26-sized
   effect is 9.3% at the actual sample (417 planets), rising to only 18.7%
   even at 8x (~3,300 planets); the null false-positive rate stays in
   2.5-11.7% across all scales. Full table in ROADMAP.md / power_analysis.json.

### Model v2 fits (mixture + sandwich SEs) [VERIFIED against saved JSONs]
| Fit | N pl | beta_age/IQR +- se | z |
|---|---|---|---|
| M x v_tan | 417 | +0.40 ± 0.36 | +1.1 |
| M x abs(W) | 317 | −0.15 ± 0.31 | −0.5 |
| M x v_tan local (<200pc) | 311 | +0.28 ± 0.38 | +0.7 |
| M x v_tan + Fe/H | 417 | +0.23 ± 0.35 | +0.7 |
| FGK control x v_tan | 2848 | +0.04 ± 0.09 | +0.4 |

Sanity anchors hold: Fe/H +0.76±0.25 (metal-rich -> more sub-Neptunes, matches PAST III).

### Bottom line (unchanged by any fix, now properly supported)
No detectable age evolution of the radius valley around M dwarfs at current sample
sizes; the FGK control is itself a precise null under our archival pipeline; and the
power analysis quantifies why: literature-sized effects are undetectable below
thousands of planets. The project's deliverable is exactly its promised two-outcome
structure, resolved as: "underpowered" — with a quantitative required-N target for
PLATO/HWO-era samples.

### Still open before submission
- Light-curve rotation periods for cool hosts (10) to upgrade the ordinal kinematic
  age axis to true gyro ages.
- Full occurrence-rate completeness treatment (inverse detection efficiency or
  injection-recovery) remains the largest unmodeled systematic; until then all
  results are framed as archival-count analyses with documented caveats.
