# P3 — End-to-end Gaia DR3 bound-pair crossmatch: ATTEMPTED, NOT COMPLETED

**Date:** 2026-11-21.
**Status:** `failed_crossmatch_endpoints_blocked` — **no crossmatch was performed; no counts, fractions, p-values, or systematic estimates from sample data are reported. Nothing here is fabricated.**

## Catalog identified (this part succeeded)

From the vault note `mnras-000-126-2021` (El-Badry, Rix & Heintz 2021, *A million binaries from Gaia eDR3*, MNRAS 506, 2269; DOI 10.1093/mnras/stab323; arXiv:2101.05282), section 7 "Catalog Description":

- **VizieR catalog: `J/MNRAS/506/2269`** — "Validation of Gaia EDR3 parallaxes".
- Landing page (reachable, HTTP 200): `https://cdsarc.cds.unistra.fr/viz-bin/cat/J/MNRAS/506/2269`
- Tables: **`catalog.dat`** (the binary catalog, 1,817,594 pairs) and `shift.dat` (shifted chance-alignment catalog, 517,993 pairs).
- Column names (note: these differ from the `source_id1`/`source_id2`/`probability` names assumed in the task):
  - Component identifiers: **`source_id_o1`, `source_id_o2`** (suffixes `_o1`/`_o2` = primary/secondary; primary is the brighter `G` component). All `gaiaedr3.gaia_source` columns are copied over for both components with these suffixes.
  - Separation: **`pairdistance`** (angular separation, degrees) and **`sep_AU`** (projected separation, au).
  - Bound/chance flag: **`R_chance_align`** (chance-alignment probability; high-confidence bound pairs selected with a low threshold, e.g. `< 0.1`; the paper's own high-purity selections use `< 0.01`).
  - Other: `binary_type` (e.g. MSMS, WDMS), plus DR2 cross-matched columns prefixed `dr2_`.
- eDR3 source IDs are the same as DR3 source IDs, so an ID-based match against our host IDs is scientifically valid — we simply could not execute it.
- Archive copy: Zenodo record **4435257**, DOI `10.5281/zenodo.4435257`, "Wide binaries from Gaia eDR3".

## Host list prepared (not used — match never ran)

- `radius-valley-age-evolution/data/gaia_hosts_M.csv`: **301** hosts with Gaia DR3 `source_id`.
- `radius-valley-age-evolution/data/gaia_hosts_FGK.csv`: **2,123** hosts with Gaia DR3 `source_id`.
- Planned union: **2,424** source IDs; planned filter: `source_id_o1 IN (...) OR source_id_o2 IN (...)` with `R_chance_align < 0.1`; planned chunk size 200 IDs.
- **These counts describe the input lists only. No matched-host count exists.**

## Exact endpoints attempted and their failures

All attempts were made on 2026-11-21 from `F:\Research_Radius-Valley-Age-Evolution` using Python 3.12.10 (`requests`, `astropy`, `astroquery` all installed).

### CDS VizieR query endpoints — blocked
| Endpoint | Request | Result |
|---|---|---|
| `https://vizier.cds.unistra.fr/viz-bin/asu-tap-sync` | GET, and POST with `REQUEST=doQuery&LANG=ADQL&FORMAT=json\|csv`, query against `tap_schema.tables` filtered on `506/2269`, `El-Badry`, `million binaries` | **HTTP 404**, identical 4,245-byte "Page unavailable — CDS" HTML block page. |
| `https://vizier.u-strasbg.fr/viz-bin/asu-tap-sync` | same POST/GET | **HTTP 404**, same 4,245-byte "Page unavailable — CDS" block page. |
| `https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=J/MNRAS/506/2269&-out.max=1` | GET; also with browser `User-Agent`; also for URL-encoded source and for control catalog `II/246/out` | **HTTP 404**, same block page for every variant — the block is endpoint-wide, not catalog-specific. |
| `https://cdsarc.cds.unistra.fr/viz-bin/ReadMe/J/MNRAS/506/2269?format=html&tex=true` | GET | **HTTP 404**, same block page. |
| `https://cdsarc.cds.unistra.fr/ftp/J/MNRAS/506/2269/ReadMe` | GET | **HTTP 404**, same block page. |
| `https://cdsarc.cds.unistra.fr/viz-bin/cat/J/MNRAS/506/2269` | GET | **HTTP 200** — landing-page metadata only; it describes the catalog but offers no query result. |
| `https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync` and `/tap/tables` | GET, HTTPS and HTTP | **`ConnectTimeout`** (Max retries exceeded, port 443 / port 80). |
| `http://tapvizier.u-strasbg.fr/TAPVizieR/tap/sync` | GET | **`ConnectTimeout`** (port 80). |
| `https://vizier.cds.unistra.fr/viz-bin/votable` (via `astroquery.vizier.Vizier(catalog='J/MNRAS/506/2269', row_limit=1).query_constraints()`) | astroquery POST | **`HTTPError: 404 Client Error: Not Found`** — same endpoint-wide block. |

### VizieR mirrors — catalog absent / unreachable
| Endpoint | Result |
|---|---|
| `http://vizier.nao.ac.jp/viz-bin/asu-tsv?-source=J/MNRAS/506/2269&-out.max=1` | **HTTP 200** (VizieR/7.6.1) but body reports `#INFO Error=Table or Catalog not found: J/MNRAS/506/2269`. Same for `.../catalog` and `.../shift`. The mirror is reachable but does not host this catalog. |
| `https://vizier.nao.ac.jp/viz-bin/VizieR?-source=J/MNRAS/506/2269` | **`ConnectTimeout`** (port 443). |

### Zenodo route — reachable, but transfer rate prohibitive
| Endpoint | Result |
|---|---|
| `https://zenodo.org/api/records/4435257` | **HTTP 200**. Files: `all_columns_catalog.fits.gz` **1,419,826,533 bytes (~1.42 GB)**; `all_columns_catalog_shift.fits.gz` 354,073,286 bytes; `num_neighbors_edr3.py` 3,616 bytes; `find_binaries_edr3.py` 14,017 bytes. |
| `https://zenodo.org/api/records/4435257/files/all_columns_catalog.fits.gz/content` with `Range: bytes=0-9999999` | **HTTP 206 Partial Content**, but 10 MB took **135.4 s ≈ 0.07 MB/s**. At this rate the 1.42 GB catalog would need ~5.6 hours — not feasible. |
| Range tests on `https://zenodo.org/records/4435257/files/all_columns_catalog.fits.gz?download=1` and the API content URL (`Range: bytes=0-19999999`) | **Aborted by instruction** before producing a result. No throughput number from these tests is claimed. |

Per the orchestrator's final instruction, all network attempts were stopped and no further fetching was done.

## Results

| Quantity | Value |
|---|---|
| Bound-pair catalog actually queried | **none** |
| Crossmatch executed | **no** |
| M hosts matched to ≥1 bound companion | **not available** |
| FGK hosts matched to ≥1 bound companion | **not available** |
| Young/old binary fractions and Fisher/chi-square p-value | **not available** |
| Measured age-differential binary incidence | **not available** |
| induced change in SN fraction / AME residual | **not available from sample data** |

Because the crossmatch never ran, **there is no sample-specific measurement of the age-differential companion incidence**. The existing literature-based estimate in `P3_multiplicity_budget.md` (Δc ≈ 25 pp observed-to-field contrast; induced shift ≈ 0.6–2.4 pp on the SN fraction, short of the +4.14 pp residual) remains the only statement available, and it is explicitly *not* a measurement in our sample.

## Verdict

- **The multiplicity disclosure in the manuscript must stay.** Nothing in this audit removes it: the end-to-end crossmatch that was supposed to quantify the effect in our sample could not be executed because CDS/VizieR query endpoints are blocked from this environment, the available mirror does not carry the catalog, and the Zenodo copy is too large to transfer at the observed rate.
- Next action if this is to be resolved later: obtain the catalog through an unblocked route (e.g. an institutional VizieR mirror that hosts `J/MNRAS/506/2269`, a local cached copy of `all_columns_catalog.fits.gz`, or CDS TAP access from a different network), then run the chunked ID query on `source_id_o1`/`source_id_o2` with `R_chance_align < 0.1`.

## Related files

- Machine-readable status: `results/audit/P3_bound_pair_summary.json`
- Literature budget (top note updated): `results/audit/P3_multiplicity_budget.md`
- No crossmatch script was written: the orchestrator's final instruction was to stop all network attempts and record the failure only, so the intended chunked-ID query was never implemented or run.
