# Manuscript

> **Manuscript source (`.tex`/`.bib`) is withheld from this public repository
> until publication**, to protect authorship. The figures, revision log, and
> numerical results below are shared; the paper text itself will appear via the
> journal / arXiv record, which carries the authoritative author list and date.

## Results the manuscript reports (high-replicate Monte Carlo)

- gate: 1,000 replicates per arm (real percentile 45.6% null / 81.2% signal,
  verdict FAIL inconclusive);
- power grid: 500 replicates per cell (8x SWEET detection 23.2%, strong 60.6%);
- extrapolated 80%-power scale ~1.7e4 planets, explicitly model-dependent.

## Files

- `fig1_ess.png`, `fig2_gate.png`, `fig3_period.png`, `fig4_power.png` — figures
  (Fig 2 and Fig 4 are the high-replicate versions; identical files are in
  `../figures/`)
- `aastex701.cls`, `aasjournalv7.bst` — AASTeX class and bibliography style
  (kept so the build environment is documented; the `.tex`/`.bib` are private)
- `FIX_LOG.md`, `REVISION_NOTES.md` — revision logs with the final
  numerical-consistency verification results

## Citing

Cite the paper via its journal/arXiv record (author list + date are fixed
there) and this repository via `CITATION.cff`.
