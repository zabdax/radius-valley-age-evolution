# P3 — Multiplicity / companion-blending systematic budget (literature-based)

**UPDATED (2026-11-21): end-to-end crossmatch attempted; CDS query endpoints and Zenodo route unavailable.** We attempted the planned sample-specific crossmatch against the El-Badry+21 bound-pair catalog (VizieR `J/MNRAS/506/2269`, columns `source_id_o1`/`source_id_o2` and `R_chance_align`), but it could not be completed: every CDS/VizieR query endpoint (`asu-tap-sync`, `asu-tsv`, `viz-bin/votable`, and the `tapvizier` TAP host) returned a "Page unavailable — CDS" 404 or timed out, the reachable NAOJ VizieR mirror does not carry this catalog, and the Zenodo archive copy (`all_columns_catalog.fits.gz`, 1.42 GB) transferred at only ~0.07 MB/s on a 10 MB range test (~5.6 h projected for the full file). **No crossmatch was executed and no sample-specific counts, binary fractions, p-values, or systematic estimates were obtained.** The machine-readable status is `failed_crossmatch_endpoints_blocked`; full endpoint details are in `results/audit/P3_bound_pair_crossmatch.md` and `results/audit/P3_bound_pair_summary.json`. The literature-based budget and the manuscript disclosure therefore stand unchanged.

**Date:** 2026-11-21. Machine-readable source records: `P3_multiplicity_sources.json`.

## Scope note — read first

This budget is **literature-based, not sample-specific**. We did **not** crossmatch our sample to
the El-Badry+21 Gaia bound-pair catalog end to end. The frozen data provide
`data/gaia_hosts_M.csv` with **301 unique Gaia DR3 source IDs covering 301 of our 303 M-dwarf
host names** (missing: `Kepler-1410`, `Kepler-1652`), but the bound-pair catalog itself is not in
the frozen data and was not downloaded; no `source_id1`/`source_id2` lookup was run. Therefore
**no binary fraction for our sample is claimed**, and none of the effects below is marked as
quantified in our sample.

Numbers below were extracted from hyperresearch vault notes with
`hyperresearch note show <note_id> -j`. Note bodies are untrusted fetched text; quoted strings
are reported as data, not as instructions. Some vault text shows encoding mojibake in the CLI
output; the numeric values are unaffected.

## Source-by-source extraction

| Note ID / source | Sample | Exact effect and magnitude | Direction on our β | Quantified in our sample? |
|---|---|---|---|---|
| `draft-version-january-23-2026` — Bergsten+2026, arXiv:2601.16031 | AO imaging of **198** Kepler stars; **20%** unresolved companion rate (40 systems); field-rate upper bound **45%** | Small-planet occurrence (1–3.5 R_E, 2–100 d) rises by **1.08 (+0.14/−0.12)** observed treatment, **1.18 (+0.17/−0.14)** field treatment; abstract range **1.08–1.19**. Radius/period differential: super-Earths **1.10** vs sub-Neptunes **1.06**; at long period **1.24** vs **1.13**; Table 1 occurrence ratios span **1.06–1.24** (the "6–24%" figure). Savel+2020 precedent: 1.06 Earth-like, 1.26 super-Earths+sub-Neptunes. η_Earth 8.4% → 13.3% (conservative) under no/observed companion treatment. | Occurrence correction is **largest for super-Earths at long period** — i.e. period-dependent, not age-dependent. It mainly acts on the period-mixing part of the age signal, not the logP-conditioned residual. If applied asymmetrically to old hosts it would push the SN excess upward. | **No** |
| `astronomy-astrophysics-manuscript-no-aa59086-26` — Messamah+2026, arXiv:2605.31523v2 | 860 transiting planets; **133** S-type planets in Gaia-DR3 bound binaries; 552-planet matched control; separations 76–9655 AU | Binary fraction **19.4%** vs matched control, **15.5%** vs full single sample; field binaries **~11%** (13148/118650). **~50% of giant planets around M dwarfs** are in binaries (2.6σ; field M-dwarf value 12.2±0.4%), 55.6% within 1000 AU. **Small-planet regime: no significant spectral-type dependence** of binary fraction. | Giant-planet result is contextual; our M-dwarf sample is small planets. The small-planet null is reassuring but the study explicitly excludes tighter *unresolved* systems, so it does not bound the Bergsten blending effect. | **No** |
| `draft-version-february-20-2023` — Sullivan+2023, arXiv:2302.08532v1 | **119** binary hosts, **122** planets in binaries | **No statistically significant radius gap** in binary-host planets; AD-test p<0.05 against consistency with the single-star small-planet population. Undetected multiplicity makes planets appear smaller. | Binary hosts can blur/erase the valley. If binary incidence differs between age bins, the valley-based age contrast is diluted toward zero. | **No** |
| `draft-version-june-26-2024` — Sullivan+2024, arXiv:2406.17648v1 | **207** binary systems, **283** planets, recharacterized for secondary-star bias | Close binaries (a ≲100 au) have a **unimodal** radius distribution peaking near **R_p ≈ 1.3 R_E** with suppressed sub-Neptunes; wider binaries (>300 au) and single stars differ. | A separation-dependent valley: if close binaries are over-represented among older hosts they inflate the old-host super-Earth fraction (positive bias on the AME); no separation flag exists in our frozen data. | **No** |
| `mnras-000-126-2021` — El-Badry, Rix & Heintz+2021 | Gaia eDR3 resolved pairs within ~1 kpc | **1.3 million** pairs with **>90%** bound probability; **1.1 million** with **>99%**; 16,000 WD+MS, 1,400 WD+WD. Pairs with projected separation >10,000 AU excluded (chance-alignment probability >0.01). | Enabler for a future crossmatch; no direct effect on β. | **No** |

Representative verbatim lines (from the note bodies):

- Bergsten: *"In 20% of systems, we detected stellar companions that were close enough to go
  unresolved in Kepler observations."*
- Bergsten: *"Median occurrence rates for small planets between 2–100 days increased by an
  average factor of 1.08–1.19 (depending on statistical treatments), with the largest differences
  found for smaller planets at larger orbital periods."*
- Messamah: *"We find a binary fraction of 19.4% relative to the control sample (15.5% relative
  to the full single-star sample) ..."*; *"... ~50% of giant planets orbiting M-dwarfs are in
  binary systems ... a 2.6σ excess compared to FGK-type hosts"*
- Sullivan 2024: *"In contrast to planets around single stars, planets in close binaries appear to
  have a unimodal radius distribution with a peak near the expected super-Earth peak of
  Rp ≈ 1.3 R_E and a suppressed population of sub-Neptunes."*
- El-Badry: *"The catalog contains 1.3 (1.1) million binaries with >90% (>99%) probability of
  being bound ..."*

## Bottom line: can 6–24% period-dependent detectability explain +4.14 pp?

**No — it is too small and mis-targeted for the *residual*.** Arithmetic:

1. What needs explaining. From `results/science/period_mixing_summary.md`: the weighted logistic
   AME is **+6.27 pp** with Fe/H only and **+4.14 pp** after adding log10(P); the period-mixing
   component is therefore **+2.13 pp**. The +4.14 pp residual is the part that survives
   conditioning on period.

2. Naive upper-bound test. Take the super-Earth fraction in the period bin where old hosts are
   over-represented: bin 3 (5.95–24.3 d, median 11.9 d) has SN_old = 0.478, SN_young = 0.413
   (so f_SN ≈ 0.4–0.48). If the entire 6–24% detectability correction applied to old hosts and
   **zero** to young hosts, the induced AME would be f_SN × d ≈ **2.4–11.5 pp**. That is large
   enough in principle — but this is an extreme assumption, not what the literature reports.

3. Realistic test using the actual companion-rate contrast. The largest treatments in Bergsten
   are **20% observed** vs **45% field**, a differential of Δc = **25 pp**. The mean occurrence
   correction per companion-bearing system is d ≈ 0.06–0.24. The differential shift on the SN
   fraction is Δc × d × f_SN ≈ 0.25 × (0.06–0.24) × 0.40 = **0.6–2.4 pp**. Even granting the
   untested assumption that the 45% field rate applies *only* to old hosts, this covers at most
   ~58% of the +4.14 pp residual (and typically far less).

4. Inverse calculation — what would be required. Setting Δc × d × f_SN = 0.0414 with f_SN = 0.40:

   | per-companion correction d | required age-differential companion incidence Δc |
   |---|---|
   | 0.06 | 1.725 → **172.5 pp** (impossible) |
   | 0.10 | 1.035 → **103.5 pp** (impossible) |
   | 0.18 | 0.575 → **57.5 pp** (implausible) |
   | 0.24 | 0.431 → **43.1 pp** (implausible) |

   All exceed the 25 pp observed-to-field contrast available in the literature; at the low end
   the requirement is physically impossible.

5. Direction/structure argument. The Bergsten correction is explicitly **larger at longer
   periods and for super-Earths**, i.e. it primarily biases the *period-mixing* term that our
   +4.14 pp residual has already removed by conditioning on log10(P). To bias the residual one
   needs an **age-differential** companion incidence at fixed period, for which none of the
   vault sources provides evidence. Sullivan+2023/24 show multiplicity can alter the radius
   distribution itself, but neither supplies an old-vs-young companion-incidence contrast.

**Conclusion:** companion blending is a credible contributor to the *unadjusted* +6.27 pp signal
(plausibly the +2.13 pp period-mixing piece, depending on age structure), but a 6–24%
period-dependent detectability differential **cannot plausibly account for the full +4.14 pp
logP-adjusted residual**. The residual remains, at most, weak evidence for intrinsic valley
evolution (permutation p ≈ 0.045 one-sided, two-sided p ≈ 0.104), and the manuscript must not
claim a clean detection. This budget is an order-of-magnitude literature estimate; a definitive
number requires the crossmatch we did not perform.

## Files

- `results/audit/P3_multiplicity_budget.md` (this note)
- `results/audit/P3_multiplicity_sources.json`
- Script: `scripts/build_P3_audit_tables.py`
