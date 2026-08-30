# REVISION NOTES — sample701.tex (post-audit publication hardening)

Basis: high-replicate Monte Carlo audit (results/colab_audit/AUDIT_REPORT.md).
Raw replicate CSVs = ground truth; JSONs = cross-checks. Figures regenerated
from high-replicate data by results/colab_audit/regen_manuscript_figs.py.

| # | Section | Change | Reason |
|---|---|---|---|
| 1 | Abstract | Full tail rewrite: 1,000/arm gate values (45.6% / 81.2%, beta -0.0296+-0.1061), 500/cell grid, 23.2% [19.6,27.2] at 8x, extrapolation order 1e4-2e4 (point ~1.7e4) labeled model-dependent, host-diversity + age-proxy caveats, mild CI anti-conservatism note | Superseded values; Sections XVIII, IX, X, XII, VII |
| 2 | Methods 3.3 (gate) | 120 -> 1,000 replicates per arm; pointer to availability statement | High-replicate run is authoritative |
| 3 | Table 2 (tab:gate) | New values: -0.0296+-0.1061; null -0.0186 [-0.2476,+0.2246]; 45.6% [42.5,48.7]; P(null>=real)=0.544 row added; signal -0.1365 [-0.3898,+0.0836]; 81.2% [78.6,83.6]; comment rewritten (empirical MC rankings, binomial intervals) | Section XIX |
| 4 | Results 4.1 text | Percentile updates + explicit "comfortably compatible with null and signal"; 81.2nd percentile is placement, not detection | Sections III, XIX |
| 5 | Fig 2 | Regenerated from raw 1,000-replicate draws (empirical histograms replace Gaussian approximations); caption updated; figure* two-column span | Ground-truth data; user request for larger figures |
| 6 | Results 4.2 | Subsection retitled "Period Mixing: An Identified Mechanism"; caption already dataset-specific | Section XIV |
| 7 | Results 4.3 intro | 120/150 -> 500 replicates per cell; raw = ground truth statement | Section V |
| 8 | Table 3 (tab:power) | All 15 cell values updated (7.2/9.0/16.2; 8.2/9.4/27.0; 5.2/13.6/38.8; 5.6/18.8/50.2; 7.2/23.2/60.6); N 3,340 -> 3,337 | Section V |
| 9 | Results 4.3 FPR paragraph | Full rewrite: old non-monotonicity resolved as MC noise (no 8x anomaly, Fisher p>=0.24); pooled FPR 167/2500=6.68%, p~1e-4; explicit gate-vs-analytical-interval distinction | Sections VI, XX |
| 10 | Results 4.3 calibration paragraph (NEW) | Coverage 6850/7500=91.3%, range 85.2-95.4%, worst 8x strong; anti-conservative sandwich CIs; does not disappear at 8x; gate unaffected | Sections VII, XX |
| 11 | Results 4.3 power paragraph | 9.4% current; 8x SWEET 23.2% [19.6,27.2], strong 60.6% [56.2,64.9]; monotone curves | Section V |
| 12 | Results 4.3 extrapolation paragraph | N80 ~= 1.7e4 [(1.4-2.2)e4] with the x=Phi^-1(p)+1.96=1.23 arithmetic shown; host-diversity saturation caveat (SE ratios 2.37/1.73/1.36/1.12/1.00 vs 4.00/2.83/2.00/1.41/1.00); two-direction uncertainty | Sections IX, X |
| 13 | Fig 4 | Regenerated from power_analysis_rep500.json (same plotting style); caption: 23.2%, 60.6%, no scale reaches 80%, extrapolation model-dependent; figure* two-column span | Sections XXIII, V |
| 14 | Discussion 5.3 | Three-layer rewrite: empirical benchmark / model-dependent extrapolation / caveats; "lower bound planning benchmarks" removed; planet count != independent age information added | Sections XXI, X, XI, XII |
| 15 | Discussion 5.5 (clustering) | Percentile clarification: check's own 120-rep ensemble (43.3rd) vs full-sample 1,000-rep gate (45.6%) | Consistency after gate update |
| 16 | Conclusion | All numbers updated (45.6/81.2, 23.2 [19.6,27.2], ~3,337, N80 ~1.7e4 model-dependent); 8-point conclusion structure | Sections XXIV, XVII |
| 17 | Data & Reproducibility Availability (NEW section*) | Data/code/raw-CSV ground-truth statement; Colab execution config; honest disclosure that colab_highrep_gate_fast_v3.py was NOT archived (no fabricated hash); fingerprint evidence cited | Sections XVI, XXVI |
| 18 | Figures 1-4 sizing | Figs 2, 3, 4 converted to figure* (full page width); Fig 1 kept single-column | User request: larger figures, layout preserved |

## Pass 2 — final publication hardening (2026-08-30)

| # | Section | Change | Reason |
|---|---|---|---|
| 19 | Abstract | Condensed from ~380 to ~230 words; retains question, calibrated gate, 1,000/arm FGK result (45.6%/81.2%, FAIL), period mixing, 8x detection 23.2% [19.6,27.2], extrapolation order 1e4-2e4 (~1.7e4), age-precision/host-diversity phrase; removed CDPP counts, audit taxonomy, CI diagnostics | Directive 2 |
| 20 | Abstract + all | "recovered in only 23.2%" -> "detected in only 23.2%" | Directive 3 (detection rate, not recovery) |
| 21 | Figure 4 | y-axis 0-125 -> 0-100 (detection rate cannot exceed 100%); regenerated | Directive 4 |
| 22 | Figure 4 | legend moved lower-right -> center-left (it occluded the SWEET and null curves near 4x-8x); visually verified | Directive 18 |
| 23 | Discussion 5.5 | Reworded to make explicit that 43.3% (robustness check) and 45.6% (primary gate) come from different calibration ensembles | Directive 5 |
| 24 | Acknowledgments | "adversarial review passes" -> "adversarial audit passes" | Directive 6 |
| 25 | Conclusion | "The honest headline result is inconclusive." -> "The primary statistical result is inconclusive." | Directives 7, 20 |
| 26 | Methods 3.3 + Results 4.3 | "deliberately conservative" -> "intentionally strict" / "conservative" (style) | Directive 20 |
| 27 | Introduction | "honest statistical adjudication" -> "rigorous statistical adjudication" (style) | Directive 20 |
| 28 | Availability (tex, README, colab/README) | Executed Colab script `colab_highrep_gate_fast_v3.py` obtained from the author and archived in `colab/` (SHA-256 fed4824d...a69e59); the previous "not preserved" disclosure replaced with the archived-script statement + hash; fingerprint evidence retained as supporting verification | The reproducibility limitation recorded in the audit is now closed |

## Final verification results (pass 2)

- Stale values: 18.7 / 64.0 / 90.8 / "43.3rd percentile of the null" / "83-98" / "lower bound" / "honest headline" / "recovered in only" / "deliberately conservative" -> **0 occurrences each**.
- Intentional old-run references retained (correctly): 11.7% and 120-replicate mentions only in the FPR history sentence (4.3), the clustering robustness check (5.5), and the availability statement's fingerprint evidence; 43.3% only in those two contexts.
- High-replicate values present: 23.2% (7x), 60.6% (4x), 45.6% (5x), 81.2% (4x), P=0.544 (2x), 6.68% (1x), 91.3% (1x), N80 ~1.7e4 (4x), 1,000 replicates/arm (4x), 500 replicates/cell (4x).
- Figures 2 and 4 visually inspected: Fig 2 shows empirical 1,000-replicate histograms with correct 45.6/81.2 annotation, no detection implication; Fig 4 has y-axis 0-100%, 80% reference line, unoccluded curves (legend moved), 23.2%/60.6% at 8x, correct annotations, title notes 500 replicates/cell.
## Pass 3 — figure layout & references page (2026-08-30)

| # | Section | Change | Reason |
|---|---|---|---|
| 28 | Preamble | Float placement tuning (topfraction/dbltopfraction 0.95, textfraction 0.05, floatpagefraction/dblfloatpagefraction 0.80) | Prevents lone figures from occupying full float pages |
| 29 | Figures 1-4 | \plotone (hardcoded 0.85*linewidth) replaced with explicit \includegraphics widths: fig1 0.72\textwidth, fig2 0.80\textwidth, fig3 0.98\textwidth, fig4 0.85\textwidth; all figures now two-column figure*[!t] with \centering | \plotone inside figure* left images one-column small and pushed floats onto near-empty float pages (Figure 3 took a full page) |
| 30 | References | \clearpage before \bibliography — references start on a standalone page and all floats are flushed before them | User request |

Compiled-page expectation: every figure spans the full text width at 2.5-3.3 in height and is placed at the top of a text page, never alone on a float page; references begin on their own page.

PDF compilation status: NOT performed locally (no LaTeX distribution installed: pdflatex/xelatex/lualatex/tectonic/latexmk/miktex all absent). Source-level QA passed (balanced environments, labels/refs consistent, document closes correctly). Compile on Overleaf (aastex701.cls included) before submission.

Deliberately unchanged: sample construction, completeness weighting, mixture
model specification, gate rule, period-mixing analysis numbers (not superseded),
host-clustering robustness check numbers (separate experiment), audit trail
content, citation structure, spectral-type caveats (already present).
