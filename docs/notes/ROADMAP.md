# M-dwarf vs FGK Radius-Valley Age Evolution — Project Roadmap

> **RETRACTION NOTICE (2026-08-26, supersedes labels below):**
> The FGK control result "+0.04±0.09 (precise null)" is RETRACTED as a label.
> It lies within ~0.4 sigma of BOTH zero and the literature +0.07 effect and
> therefore rules out neither; "precise" overstated it, and comparing our
> mixing-fraction-per-IQR coefficient to SWEET-Cat's age slope across
> parameterizations without calibration is statistically invalid. Until the
> completeness-weighted calibration gate (inject-signal vs inject-null) is
> run and passed, NO result in this document may be described as "null",
> "precise", or "confirming" anything. Current honest status of every
> age-evolution estimate: INCONCLUSIVE pending Task 2 calibration.

> **TASK 2 GATE RESULT (2026-08-26): FAIL (inconclusive) — definitive run.**
> Pipeline verified healthy first (injection-recovery diagnostic exposed an
> archive radius-error outlier corrupting synthetic calibration; physical
> clip applied; replicates now track injected signals correctly).
> Primary control (FGK Kepler-only subset, DR25-based trimmed weights,
> N=2135, ESS=1645): real beta_age = -0.046 +- 0.106 per IQR.
> Null-injected distribution [-0.248,+0.221] contains it at P=0.57;
> signal-injected [-0.355,+0.068] barely reaches it at its 91st pct.
> The test cannot discriminate literature-size effects from zero even
> under completeness weighting -> inconclusive, per protocol.
> TASK 2b: wrong-sign binned artifact SURVIVES weighting (p=0.030 vs
> 0.023 unweighted) while hierarchical says ~0 -> discrepancy flagged,
> not resolved by convenience.
> TASK 4: no evidence selection geometry SUPPRESSES signals (injections
> recovered at documented attenuation; distance strata show no monotonic
> distortion). Instead found a METALLICITY CONFOUND: beta_age swings
> monotonically from -0.20+-0.22 (Fe/H<-0.05) to +0.22+-0.23 (Fe/H>+0.05)
> -- age-proxy slopes are entangled with metallicity covariance.
> TASK 3 branch taken: FAIL -> M-dwarf rerun and manuscript BLOCKED;
> diagnostic report delivered; awaiting direction on minimum viable fix.
> PHASE A (Task 2b CLOSED 2026-08-26): weighted binned test re-run within
> metallicity terciles on the same control set. Unstratified z=+2.16/p=0.03
> collapses to per-stratum z=+0.49/+1.55/+0.62; inverse-variance pooled
> dSN=+0.042+-0.028, z=+1.53, p=0.13; heterogeneity Q=0.69 (p=0.71).
> The significant wrong-sign artifact WAS metallicity-confounded.
> Hierarchical (Fe/H-adjusted) and binned (now adjusted) agree at ~zero.
> Residual +1.5-sigma tilt reported without strong interpretation.
> Task 2b: RESOLVED. Phase B (Berger+20 isochrone-age axis) NOT started --
> awaiting explicit authorization.
> TASK 2b-CONFIRM (full-power regression): metallicity hypothesis REJECTED --
> Fe/H-adjusted logistic at full N leaves artifact intact (+6.29+-2.90 pp,
> p=0.03). True driver found: PERIOD MIXING. Old-vt hosts host longer-P
> planets; SN fraction rises steeply with P against the FIXED 1.88 Re cut;
> adding logP collapses artifact to +4.16+-2.61 pp (z=1.61, n.s.) with
> beta_logP=+2.38+-0.17. Explains why hierarchical model (which includes
> logP) was always null. Residual +1.6-sigma tilt noted, uninterpreted.
*(rev 3, 2026-08-26 — post-audit corrections applied)*

## Question
Does the radius valley evolve with stellar age around M dwarfs as it does
around FGK stars? FGK effect quantified (SE/SN ratio 0.51→0.64 across 3 Gyr;
valley shallows; Kamulali+26 arXiv:2601.12396). Only M-dwarf test to date:
Gaidos+24 (arXiv:2404.11022), directionally consistent decline, p=0.08,
N=117 KOIs / 74 hosts. Gap = power analysis against the quantified baseline.

## Design decisions log

| Decision | Choice | Rationale |
|---|---|---|
| Planet sample | NEA `ps`, default_flag=1, controv=0, Transit, 0.5<Rp<4 R_E, P<100 d | one row per planet |
| M sample | host Teff<4200 K: **419 pl / 303 hosts** | |
| FGK control | 4200-7000 K: **2880 pl / 2150 hosts** | positive control |
| Valley floor M | **1.85 ± 0.03 R_E** (modes 1.23/2.08) | stable over KDE sigma 0.08-0.15; edge-guarded |
| Valley floor FGK | **1.88 R_E** | matches Fulton gap |
| Gyro model | Angus+19 slow sequence (`gyrointerp`) | Teff-based, cluster-calibrated |
| Gyro validity | Teff 3800-6200 K (hard floor, models.py L534) | below -> NaN |
| Gyro errors | archive st_rotperr1 clipped to [5%,35%] of Prot; fallback 5% | guards overconfident/absurd errors |
| Kinematic axis | v_tan (all hosts) / abs(W) (RV subset, n=237 M) | ordinal age ranking ONLY |

## Validation gate — HONEST STATUS: NOT PASSED
The gate below was defined before the first science run. Current status:
- (a) FGK control shows expected SN-decline: **NOT MET as significance**.
      Direction-only at p=0.85 (gyro ages, N=128). The only significant
      FGK kinematic result is WRONG-SIGN and confounded (see headline).
- (b) Kinematic proxy correlates with anchor ages (rho>0.5): **UNTESTABLE /
      effectively failing**. True anchor coverage in-sample is only
      AU Mic (22 Myr, vtan 21 km/s) and K2-25 (650 Myr, 26 km/s) — both at
      low vtan, consistent but n=2. Earlier printed correlations (rho=-0.50,
      +0.40) were artifacts of a substring-matching bug (fixed in 07) that
      mislabeled unrelated hosts (TOI-1075 as "TOI 1227", HIP 113103 as
      "HIP 67522", K2-250 as "K2-25"); those numbers are void.
- (c) M-direction check vs Gaidos: met direction-only on vtan (p=0.92).
**Decision recorded**: work continues on PIPELINE BUILD-OUT (hierarchical
model, completeness treatment, power analysis) because these are required
regardless; no scientific CLAIM about age evolution is made until (a)-(c)
pass on corrected statistics.

## SESSION RESULTS (2026-08-26; corrected after 3-agent audit)

### Test suite (post-fix)
| Test | N planets | Result |
|---|---|---|
| A: M x gyro age | 16 | SUGGESTIVE ONLY: SN frac 0.875 (young) vs 0.250 (old), p~0.04, Gaidos-consistent at floors 1.75/1.85/1.95 — but N=16 uncorrected for censoring/multiplicity |
| B: FGK x gyro age | 128 | direction-only decline, NULL (p=0.85) |
| C: M x vtan ranking | 417 | direction-only, NULL (p=0.92) |
| C2: M x **abs(W)** ranking | 317 | NULL (p=0.91), direction does NOT match Gaidos once the signed-W split bug was fixed |

### HEADLINE FINDING — powered FGK control produces a wrong-sign artifact
FGK x v_tan, N=2848 pl / 2123 hosts: SN frac young 0.561 -> old 0.604,
p=0.023 — SIGNIFICANT but WRONG SIGN.
Diagnostics (committed as code/08_headline_fgk_vtan.py):
- NOT metallicity ([Fe/H] medians identical; |Fe/H|<0.05 keeps signal)
- KILLED by local-volume cut (dist<200pc: 0.687 vs 0.701, p=0.77)
- Mission split (hostname-prefix mask TOI|Kepler|KIC|K2|HIP|HD|GJ|...):
  signal lives ONLY in Kepler hosts (N=2135, med 776 pc, p=0.024);
  K2/TOI/nearby-HD show nothing (non-Kepler subset p=0.55).
Interpretation: raw confirmed-planet counts vs kinematic heat conflate
survey geometry (distance-dependent completeness/targeting in the Kepler
field) with any true age signal. Published FGK declines used
completeness-corrected occurrence rates. Consequence: hierarchical model
must include P, S_in, M*, Fe/H covariates + completeness BEFORE any age
coefficient is interpreted; power analysis runs on the corrected statistic.

### AUDIT RECORD (2026-08-26, three independent subagent audits)
Core numbers verified reproducible exactly: sample sizes (live TAP recount),
valley floors (incl. sigma-sweep stability +-0.03), SE/SN ratios, headline
0.561->0.604/p=0.023, external paper figures. Bugs found & fixed:
1. C2 signed-W split (report-changing; fixed -> abs(W))
2. anchor substring matching (validation files regenerated)
3. heliocentric-vs-Galactocentric mislabel (comments; numbers fine;
   UVW cross-checked against hand-coded rotation & BPMG kinematics)
4. 06 resume crash on legacy header (fixed; file regenerated cleanly)
5. no match-quality gate: 2 spurious Gaia matches found (Kepler-1410,
   Kepler-1652: null parallax/PM) silently dropped downstream; gates +
   reject log added
6. 03 ignored archive Prot errors (fixed w/ clipping; ~5 FGK intervals
   widen materially, e.g. TOI-2458 44% error)
7. valley-floor algorithm lacked window-edge guard (added; failure mode
   demonstrated at sigma=0.20)
Known accepted limitations: no epoch-propagation of coordinates before cone
search (max sep seen 0.73", PM-driven; would bias against high-vt stars if
triggered); 60 transit planets with NULL st_teff excluded by design; one FGK
host (Kepler-515) missing from Gaia table (retry exhaustion).

## Age strategy (revised)
Primary axis = kinematics (v_tan all hosts; abs(W) where RV): mission-
independent. Gyro = warm-tail anchor only (Teff>=3800). Rotation catalogs
exhausted for cool hosts (+1 via McQuillan); own light-curve measurement
(TESS/K2/Kepler, lightkurve) is the remaining rotation route.

## POWER ANALYSIS RESULTS (11, FINAL v2 run, 2026-08-26)
Engine: bootstrap hosts/covariates from the real M sample at k x scale;
inject beta_age (per SD v_tan); radii drawn from reality-matched overlapping
class-conditionals; estimate with the full v2 mixture pipeline (fixed_su +
sandwich SEs). Generative host heterogeneity sigma_u=1.5 (scenario).
Effect sizes: SWEET = -0.14/SD (~5 pp SN-fraction decline across age halves,
the quantified Kamulali+26 baseline); strong = -0.28/SD (2x literature).
Detection = 95% CI excludes 0 with injected sign. 120-150 reps/cell.
Caveat: upsampling resamples the empirical covariate distribution -> the
required-N figures are LOWER BOUNDS.

| scale | planets | null FPR | SWEET-size detect | strong detect |
|---|---|---|---|---|
| 0.5x | ~208 | 6.7% | 6.0% | 16.0% |
| **1x (actual)** | **417** | **11.7%** | **9.3%** | **30.0%** |
| 2x | ~834 | 4.2% | 14.7% | 43.3% |
| 4x | ~1670 | 2.5% | 14.7% | 46.0% |
| 8x | ~3340 | 7.5% | **18.7%** | 64.0% |

HEADLINE CONCLUSIONS:
1. At the ACTUAL sample size, a radius-valley age effect exactly as large
   as the published FGK one is detected only 9.3% of the time -- BELOW the
   nominal false-positive rate of the test itself. The Gaidos+24
   non-detection was statistically unavoidable, independent of whether the
   effect exists around M dwarfs.
2. Even at EIGHT TIMES the current sample (~3,300 planets), the
   literature-sized effect is caught <19% of the time. The question is not
   answerable with current-generation M-dwarf catalogs by these methods;
   the target for the PLATO/HWO era is a multi-thousand-planet sample with
   controlled completeness.
3. Only effects twice the FGK baseline approach reliable detection, and
   only at 8x scale (64%) -- still short of routine.
4. Calibration healthy: null FPRs 2.5-11.7%, coverage 83-98%.

## HIERARCHICAL MODEL RESULTS (09, v2 mixture likelihood, 2026-08-26)
Model v2: observed radius modeled as two-component Gaussian mixture
(super-Earth / sub-Neptune populations) with covariate-dependent mixing
fraction pi(age, x) + mission dummies; sigma_u fixed at 0 and within-host
correlation handled by CLUSTER-ROBUST sandwich SEs. Insolation computed
physically (audit fix); GH normalization corrected (audit fix).
SELFTEST (real-regime, 24 seeds): the estimator ATTENUATES small beta_age
toward zero by ~35% (+0.37..0.49 logits/IQR bias) with ~25-30% CI coverage.
Direction: CONSERVATIVE for null conclusions. Any future positive detection
must be deconvolved via this calibration before being believed.

beta_age per IQR of age proxy (negative = Gaidos-consistent direction):

| Fit | N pl | beta_age +- se | z |
|---|---|---|---|
| M x v_tan | 417 | **+0.40 ± 0.36** | +1.1 |
| M x abs(W) | 317 | −0.15 ± 0.31 | −0.5 |
| M x v_tan, dist<200pc | 311 | +0.28 ± 0.38 | +0.7 |
| M x v_tan, +Fe/H | 417 | +0.23 ± 0.35 | +0.7 |
| FGK x v_tan (control) | 2848 | +0.04 ± 0.09 | +0.4 |

Reading:
1. NO detectable age evolution around M dwarfs under any proxy or cut;
   FGK control is a precise null too (|z|<0.5). Given measured attenuation
   is TOWARD zero, true effects could only be LARGER than fitted -- yet the
   FGK control bounds even the deattenuated FGK slope near zero at current
   precision, consistent with the literature's tiny effect size.
2. Proxy sign inconsistency persists (v_tan vs absW).
3. Sanity: Fe/H coefficient +0.76 +- 0.25 (metal-rich -> more SN; matches
   PAST III); period dependence strong.
4. The earlier v1 table (Bernoulli soft-weight model) was superseded after
   the audit+calibration work revealed ~35% attenuation from unmodeled
   population overlap; v2 models that overlap explicitly.

## POWER ANALYSIS RESULTS (11, 2026-08-26)
Engine: bootstrap hosts/covariates from the real M sample at k x scale;
inject beta_age (per SD v_tan); radii drawn from overlapping class-
conditionals + per-planet errors; estimate with the full 09 pipeline.
Effect sizes: SWEET = -0.14/SD (~5 pp SN-fraction decline across age
halves, the quantified Kamulali+26 baseline); strong = -0.28/SD.
Detection = 95% CI excludes 0 with injected sign. 120-150 reps/cell.
NOTE: upsampling resamples the empirical covariate distribution (no new
independent hosts) -> required-N figures are LOWER BOUNDS.

| scale | planets | null FPR | SWEET detect | strong detect |
|---|---|---|---|---|
| 0.5x | ~208 | 8.3% | 4.7% | 14.7% |
| **1x (actual)** | **417** | **9.2%** | **15.3%** | **27.3%** |
| 2x | ~834 | 5.0% | 12.7% | 36.7% |
| 4x | ~1670 | 5.0% | 16.7% | 50.7% |
| 8x | ~3340 | (pending) | (pending) | (pending) |

HEADLINE: at the actual sample size, an effect exactly as large as the
published FGK one would be detected only ~15% of the time -- i.e., the
Gaidos+24 non-detection around M dwarfs was essentially guaranteed by
sample size alone, independent of whether the effect exists. Even 2-4x
growth leaves the SWEET-size effect mostly invisible; the "strong" (2x
literature) effect needs >=4x (~1700 planets) for coin-flip power and
presumably 8x+ for routine detection. CONCLUSION FOR THE PAPER: the
current generation of M-dwarf samples CANNOT distinguish presence from
absence of FGK-like radius-valley age evolution; the target for the
PLATO/HWO era is a clean, quantitative number (see table).
Calibration checks inside the run: null false-positive rates 5-9%
(nominal), coverage 89-96% across cells.

## HIERARCHICAL MODEL RESULTS (09, 2026-08-26)
Model: soft-classification logistic (latent true radius -> Phi weight),
host random intercepts (Gauss-Hermite), covariates logP, logS(computed),
mission dummies; MAP+Laplace CIs. SELFTEST: 6/6 recovery on synthetic data.
beta_age quoted per IQR of the age proxy (negative = Gaidos-consistent):

| Fit | N pl | beta_age +- se | z |
|---|---|---|---|
| M x v_tan | 417 | **+0.30 ± 0.29** | +1.06 |
| M x abs(W) | 317 | −0.15 ± 0.28 | −0.53 |
| M x v_tan, dist<200pc | 311 | +0.26 ± 0.30 | +0.87 |
| M x v_tan, floor 1.95 | 417 | +0.43 ± 0.31 | +1.37 |
| M x v_tan, +Fe/H | 417 | +0.23 ± 0.27 | +0.85 |
| FGK x v_tan | 2848 | +0.12 ± 0.12 | +0.97 |

Reading:
1. NO detectable age evolution around M dwarfs under any proxy/cut —
   consistent with the binned tests; the two honest outcomes are now
   narrowed to "underpowered" vs "weaker-than-FGK effect".
2. Sign INCONSISTENCY across proxies (v_tan positive, |W| negative)
   reinforces that no claim is possible.
3. Residual positive-vtan tendency mirrors the 08 artifact even after
   mission dummies (which absorb HUGE offsets: Kepler-vs-TOI ≈ −4 logits!)
   -> distance/completeness systematics still unmodeled; full occurrence-
   rate machinery or IDEM weights needed before believing ANY sign.
4. Sanity checks pass: Fe/H coefficient +0.83±0.24 (metal-rich -> more SN,
   matches PAST III SE/SN-vs-Fe/H); strong period dependence recovered.
5. Host overdispersion sigma_u~1.4-1.6 logits — multi-planet correlation
   handled, not driving results.

## Pipeline status
- [x] 01_download.py        M 419 pl; FGK 2880 pl (verified live)
- [x] 02_explore.py         floors + edge guard + dip depth
- [x] 03_gyro_ages.py       v2 rerun DONE (archive errors clipped):
                            M 12 usable (median 2.50 Gyr; ceiling pile-up
                            fell 58%->33%); FGK 79 usable (median 1.93,
                            ceiling 18%). Test-B fractions unchanged.
- [x] 05_literature_rotations.py  McQuillan (+1 host); download now scripted
- [x] 06_gaia_crossmatch.py M: legacy file restored & gated offline
                            (06b): 303->301 (Kepler-1410/-1652 rejected);
                            FGK: 2149->2123 (26 faint-Kepler rejects).
                            RV coverage 78% / ~67%. Kepler-515 absent
                            (retry exhaustion). NOTE: live gate initially
                            had an off-by-one column bug (rejected all
                            no-RV hosts); fixed, offline gate used instead.
- [x] 07_kinematic_ages.py  ran for BOTH tags; exact-matching anchors;
                            anchor coverage n=2 (AU Mic, K2-25) -> gate(b)
                            untestable until more young anchors exist
- [x] 04_validate.py        suite v3 (post-fix) - see table above
- [x] 08_headline_fgk_vtan.py  reproduces every headline number exactly
- [x] 09_hierarchical_model.py DONE — soft-classification logistic w/ host
        random effects; SELFTEST 6/6; results table above
- [ ] 10_lightcurve_prot.py measure Prot for cool hosts from light curves
- [ ] 11_power_analysis.py  inject FGK-sized signal into synthetic M samples
