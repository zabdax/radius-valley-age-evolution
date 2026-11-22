# Repository status: which results are authoritative

This repository contains two generations of material. Read this file first so
that historical numbers are not mistaken for the current ones.

## Current (v2 / Q1) results

**`repro/v2/` is authoritative for the v2 manuscript.** Use the JSON/CSV files
there, plus the top-level `data/` (frozen inputs) and `code/` (shared analysis
machinery), for any current claim. Verify the bundle with:

```bash
cd repro/v2 && sha256sum -c provenance/BUNDLE_SHA256.txt
```

The release metadata is `v2.0.3` (`.zenodo.json`, `CITATION.cff`).

## Historical material (do not use for v2 claims)

The v1 notes, pre-audit draft, and earlier revision records were removed from
the active repository. They remain available in the immutable Zenodo v2.0.1
snapshot (doi:10.5281/zenodo.23007310), so the audit trail is preserved
without cluttering the current tree.

Still in place but explicitly historical:

- `results/colab_audit/AUDIT_REPORT.md` — the v1 audit report (banner added at top)
- top-level `figures/` — earlier figures plus the two active animations; the v2
  publication figures are in `repro/v2/figures/`

## Superseded -> current values

| Quantity | Older docs (historical) | v2 value (`repro/v2/`) |
|---|---|---|
| Null false-positive rate | 2.5--11.7% (120 replicates) | 5.2--8.2% (500 replicates) |
| Coverage | 83--98% | 85.2--95.4% (pooled 91.3%) |
| 8x literature-sized detection | 18.7% (v1), then 23.2% | 23.2% [19.6, 27.2] |
| 80%-power sample scale | ~1.7e4 planets [1.4, 2.2]e4 | **bracket 1.7e4--4.7e4** planets |
| Kamulali 3 Gyr split | 136 young / 687 old = 19.8th percentile | 687 is the **old** subsample; 136/823 = 16.5th percentile, 1.80 SD, −0.126/SD |
| Cool-host age proxy | "unvalidated" | pipeline-domain only (rho = 0.25, n = 12); external TIME Table ordinal validation rho = 0.372, p = 0.0023, n = 65 |
| Gate verdict | FAIL (inconclusive) | unchanged: FAIL |

The 80%-power number in particular is a **model-dependent extrapolation**, not
an empirically demonstrated requirement. v2 reports both the Gaussian
independent-host estimate (1.7e4) and the scaling measured in the grid (4.7e4),
and recommends quoting the bracket.

## Two easily confused Gaia counts

- **2,424** host source IDs were re-queried for the RUWE / parallax-S/N
  sensitivity gate (301 cool + 2,123 FGK).
- **2,420** unique `source_id` rows are present in the frozen quality table
  `data/gaia_quality_DR3_v2.csv`, which is why the audit files under
  `repro/v2/results/audit/` quote 2,420.

Both are correct: the second is the deduplicated row count of the downloaded
quality table, not a correction to the first.

## Manuscript source

LaTeX source (`.tex`, `.bib`) and the compiled manuscript PDF are deliberately
not distributed. Figures, data, code, results, and provenance records are.
