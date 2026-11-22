# Zenodo upload checklist for Q1 paper v2 (2026-11-21)

Preparation only. No DOI has been minted and no Zenodo record has been created.

## Suggested record title

Replication package for "Can Current Archival Data Detect Age Evolution of the Radius Valley Around M Dwarfs? A Calibrated Null Result"

## Files to upload (46 files, 1532762 bytes)

Exact list (paths relative to manuscripts/q1-paper/v2-2026-11/; checksums in ZENODO_FILE_MANIFEST.csv and ZENODO_RELEASE_MANIFEST.json):

- `aasjournalv7.bst`
- `aastex701.cls`
- `BASELINE_SHA256.txt`
- `figures/fig1_radius_valley.png`
- `figures/fig2_gate.png`
- `figures/fig3_period.png`
- `figures/fig4_power.png`
- `main.tex`
- `NOVELTY_AND_RIGOR.md`
- `README.md`
- `references.bib`
- `results/audit/P1_period_mixing_stats.json`
- `results/audit/P1_statistical_conventions.md`
- `results/audit/P3_host_demographics.csv`
- `results/audit/P3_host_demographics.md`
- `results/audit/P3_multiplicity_budget.md`
- `results/audit/P3_multiplicity_sources.json`
- `results/audit/P3_parallax_orchestrator_check.md`
- `results/audit/STOP_GATE_1_LEDGER.md`
- `results/audit/STOP_GATE_2_LEDGER.md`
- `results/science/checkpoints/curve_null.json`
- `results/science/checkpoints/curve_sig_0.json`
- `results/science/checkpoints/curve_sig_1.json`
- `results/science/checkpoints/curve_sig_2.json`
- `results/science/effect_size_derivation.json`
- `results/science/effect_size_derivation_orchestrator_check.md`
- `results/science/P2_period_results_verification.md`
- `results/science/parse_time_table.py`
- `results/science/period_baseline_orchestrator_check.md`
- `results/science/period_classifier_results.csv`
- `results/science/period_mixing_summary.md`
- `results/science/permutation_control.json`
- `results/science/power_curve_analytic_crosscheck.md`
- `results/science/power_vs_effect.png`
- `results/science/power_vs_effect_figure_data.csv`
- `results/science/power_vs_effect_raw.csv`
- `results/science/power_vs_effect_summary.json`
- `results/science/time_table_crossmatch.json`
- `results/science/time_table_limits.md`
- `results/science/time_table_mrt_source.txt`
- `results/science/time_table_parsed.csv`
- `results/science/within_period_test.csv`
- `scripts/parse_time_table.py`
- `scripts/period_mixing_analysis.py`
- `scripts/power_vs_effect.py`
- `scripts/write_partial_power.py`

Do not upload (30 files excluded):

- `results/audit/_note_astronomy-astrophysics-manuscript-no-aa59086-26.json` -- scratch note-cache artifact; filename starts with _; contains an absolute local vault path
- `results/audit/_note_draft-version-february-20-2023.json` -- scratch note-cache artifact; filename starts with _; contains an absolute local vault path
- `results/audit/_note_draft-version-january-23-2026.json` -- scratch note-cache artifact; filename starts with _; contains an absolute local vault path
- `results/audit/_note_draft-version-june-26-2024.json` -- scratch note-cache artifact; filename starts with _; contains an absolute local vault path
- `results/audit/_note_excerpts.txt` -- scratch note excerpt dump; filename starts with _; contains untrusted third-party text
- `results/audit/_note_mnras-000-126-2021.json` -- scratch note-cache artifact; filename starts with _; contains an absolute local vault path
- `results/audit/AUDIT_VERIFICATION_v2.md` -- contains absolute local Windows paths (C:\Users\..., F:\...) and a user-specific Python install path
- `results/audit/gaia_cut_summary.json` -- records an absolute local frozen-project root (F:\...)
- `results/audit/injected_effect_derivation.json` -- contains a user-specific local Python install path
- `results/audit/P1_calibration_provenance.md` -- opens with an absolute local workspace path (F:\Research_Radius-Valley-Age-Evolution)
- `results/audit/P1_counts_conventions.md` -- states paths relative to an absolute local workspace root (F:\...)
- `results/audit/P3_gaia_production_cuts.md` -- records an absolute local working root (F:\...)
- `results/audit/P4_REWRITE_PLAN.md` -- internal working/planning note, not a scientific result product
- `results/audit/PROVENANCE_LOCK.md` -- contains a user-specific local Python install path (C:\Users\...)
- `results/audit/recompute_rho_fgk_v2.json` -- records a user-specific local Python install path
- `results/audit/recompute_rho_fgk_v2.py` -- hard-codes an absolute local workspace root (F:\...)
- `results/audit/verify_archived_runs_v2.json` -- contains 16 absolute local F:\ paths to archived files and scripts
- `results/audit/verify_archived_runs_v2.py` -- hard-codes an absolute local workspace root (F:\...)
- `results/science/_tmp_raw_0.csv` -- scratch/temporary diagnostic artifact, not a result product
- `results/science/_tmp_raw_1.csv` -- scratch/temporary diagnostic artifact, not a result product
- `scripts/__pycache__/period_mixing_analysis.cpython-312.pyc` -- compiled bytecode; regenerable and interpreter-version specific
- `scripts/_explore_demographics.py` -- scratch exploration script hard-coded to an absolute local F:\ path
- `scripts/_explore_gaia.py` -- scratch exploration script hard-coded to an absolute local F:\ frozen-project path
- `scripts/_extract_notes.py` -- scratch utility hard-coded to an absolute local F:\ audit directory
- `scripts/build_P3_audit_tables.py` -- hard-codes an absolute local working root (F:\...)
- `scripts/find_table_end.py` -- one-off script hard-coded to an absolute local F:\ path
- `scripts/gaia_production_cuts.py` -- one-off verification script bound to the local frozen-project path
- `scripts/make_power_figure.py` -- hard-codes an absolute local F:\ results directory
- `scripts/quick_large_effect_check.py` -- scratch diagnostic hard-coded to an absolute local F:\ path
- `scripts/summarize_power_checkpoints.py` -- hard-codes an absolute local F:\ results directory

## Suggested licences

- Code (scripts/, analysis/recompute .py): MIT.
- Paper (main.tex, .bib, figures, compiled PDF if added): CC BY 4.0.
- Result data and audit notes: CC BY 4.0.
- astex701.cls and asjournalv7.bst are official AAS build files; confirm their redistribution terms before upload rather than assuming the MIT/CC-BY grant covers them.

## Prerequisites

- ORCID iDs (already present in main.tex): Zubayer Hasan Shaad 0009-0002-2322-8553; Md Mahdin Al Abrar 0009-0009-9645-152X. Confirm both are linked to the submitting account.
- Affiliations: Government Tolaram College; Chandpur Government College.
- Description: use the description field in ZENODO_RELEASE_MANIFEST.json (3 sentences, grounded in the manuscript).
- Keywords: Exoplanets; Exoplanet evolution; M dwarf stars; Astrostatistics; Exoplanet catalogs (UAT 498, 491, 982, 1882, 488).
- Related identifier: 10.5281/zenodo.23007310, relation IsVersionOf (switch to IsIdenticalTo only if the deposit is byte-identical reuse).
- Confirm the AASTeX build files may be redistributed, and confirm whether author email addresses in main.tex should remain visible.

## Verification after upload

Re-check every checksum in the deposited CSV from the project root:

``bash
python - <<'PY'
import csv, hashlib, pathlib
root = pathlib.Path("manuscripts/q1-paper/v2-2026-11")
bad = []
with open(root / "results/audit/ZENODO_FILE_MANIFEST.csv", newline="", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        p = root / row["relative_path"]
        if hashlib.sha256(p.read_bytes()).hexdigest() != row["sha256"]:
            bad.append(row["relative_path"])
print("MISMATCHES:", bad if bad else "none")
PY
``

Frozen data checksums (independent of the upload set), from 
adius-valley-age-evolution/data/:

``bash
cd radius-valley-age-evolution/data && sha256sum -c checksums.sha256
``

## Frozen data checksum verification (re-run for this manifest)


adius-valley-age-evolution/data/checksums.sha256: **31/31 verified, 0 mismatches, 0 missing files** (SHA-256, re-checked 2026-11-21 before preparation of this manifest). The frozen project was read-only throughout.

## Scope

The three ZENODO_* files in 
esults/audit/ are release metadata created by this curation pass and are intentionally not part of the inventoried file list or the upload set.
