# Manuscript

`sample701.tex` is the authoritative submission-ready manuscript (AASTeX
v7.01, two-column), fully updated to the high-replicate Monte Carlo results:

- gate: 1,000 replicates per arm (real percentile 45.6% null / 81.2% signal,
  verdict FAIL inconclusive);
- power grid: 500 replicates per cell (8x SWEET detection 23.2%, strong 60.6%);
- extrapolated 80%-power scale ~1.7e4 planets, explicitly model-dependent.

## Files

- `sample701.tex`, `sample701.bib` — manuscript source and references
- `aastex701.cls`, `aasjournalv7.bst` — AASTeX class and bibliography style
- `fig1_ess.png`, `fig2_gate.png`, `fig3_period.png`, `fig4_power.png` — figures
  (Fig 2 and Fig 4 are the high-replicate versions; identical files are in
  `../figures/`)
- `REVISION_NOTES.md` — complete revision log across all passes, with the
  final numerical-consistency verification results

## Compiling

Upload the whole folder to Overleaf (compiler: pdfLaTeX), or locally:

    pdflatex sample701.tex
    bibtex sample701
    pdflatex sample701.tex
    pdflatex sample701.tex

`orcid-ID.png` is required by the class for the ORCID marker. All four
figures are referenced with `\includegraphics` at explicit `\textwidth`
fractions and placed as `figure*[!t]` two-column floats; float placement
parameters are tuned in the preamble so no figure occupies a float page, and
the references start on their own page (`\clearpage` before the bibliography).