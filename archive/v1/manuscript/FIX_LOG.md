# FIX LOG — critic-driven revision (2026-09-27)

Base: `C:\Windows\Temp\opencode\manuscript_unzip\sample701.tex` (finalised draft, 444 lines).
Fixed manuscript: `sample701.tex` (edited in place) + copies
`sample701_FIXED.tex` (same folder) and
`radius-valley-age-evolution/manuscript/sample701_FIXED.tex`.

## Code fixes (authoritative tree `radius-valley-age-evolution/code/`)

1. Injection host-effect bug — opt-in corrected form, legacy default preserved:
   - `13_calibration_fgk.py`: `synth(..., independent_hosts=False)`; per-occurrence
     `np.repeat(uh, sizes)` when True; added `B_INJ_ALT = -0.127`.
   - `11_power_analysis.py`: same pattern (`SIGMA_U`, block sizes via count).
   - `14_one_per_host_check.py`: `synth_null(..., independent_hosts=False)`.
   - `15_stratified_gate.py`: `synth(..., independent_hosts=False)`,
     `--independent-hosts` CLI flag wired through `_init_worker/_one_rep/run_arm`
     legacy + spawned paths; RNG calls unchanged (stream-comparable).
   - Verified: `py_compile` exit 0; import OK; `synth` signatures show flag;
     legacy vs fixed on M base (seed 7): same N=410, 8.5% `pl_rade` draws differ.
   - Quantified waste (computed, seed 0, n=301, 2000 bootstraps):
     mean distinct hosts 190.5/301 (63.3%), ~110.5 draws (36.7%) unused in legacy.
     Manuscript Sec audits now states ≈190 distinct, ≈111/301 (≈37%) unused.

2. `--verify` still pins legacy stream (intended); corrected form is separate flag.

## Manuscript fixes (`sample701.tex`)

- A1–A4: "no estimator ... can distinguish" / "unanswerable / in principle
  regardless estimator" → "under the tested mixture model and calibrated gate
  (arms 0.19 SD)". Gate-method paragraph softened accordingly.
- B2–B3: Teff<4200K scope guard strengthened (heterogeneous <4200K, cannot test
  M3+ directly). Title unchanged.
- C1–C2: `-0.127` alternative leaves FAILs unchanged; distribution-injection
  stated as future work.
- Extrapolation: kept 1.7e4/4.7e4 bracket + "neither demonstrated" guard;
  abstract still quotes no N80.
- Sample footnote (Sec 2.1): 419/303 → 417/301 (Kepler-1410/1652, no Gaia);
  2880/2150 → 2848/2123 (32 pl/27 hosts) → Kepler-only 2135/1588. From direct
  CSV row counts + join logic in `15:135-142`, `13:113-115`.
- Gaia Sec: added parallax-S/N + uncertainty-propagation sentence; cone-radius/
  epoch recorded as not-in-manifest limitation (was TODO comment).
- Audits Sec: bug paragraph now quantified + points to `independent_hosts` opt-in.
- Availability: Zenodo DOI placeholder line (was TODO); AI tool named
  "OpenCode with Muse Spark 1.3 Free" (was generic + TODO).
- Checks: `TODO` 0 occurrences; "no estimator"/"unanswerable" 0; py_compile 0.

## Remaining before submission (not run here)

- Full 1000/arm + 500/cell rerun with `--independent-hosts` (needs ~6h Colab/CPU
  per manifest); report arm delta appendix.
- Beta-distribution injection (Kamulali ±0.11) sensitivity.
- Mint Zenodo DOI; attach prompt log; compile on Overleaf (no local TeX here).
- RUWE ≤1.4 + parallax-S/N sensitivity (code currently naive `d=1000/plx`).
