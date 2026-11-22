# Zenodo release guide for v2.0.3

This repository is connected to Zenodo through GitHub. A new Zenodo version is
created automatically when you publish a new GitHub release, so the correct way
to mint the v2.0.3 DOI is to publish a release/tag here rather than uploading
files by hand.

## Status

- Current archived release: **v2.0.3**, DOI 10.5281/zenodo.24444755
- Earlier release: **v2.0.1**, DOI 10.5281/zenodo.23007310 (preserves the
  superseded v1 material)
- The v2.0.3 GitHub release is published and Zenodo has minted the DOI above.

## Record metadata (use these verbatim)

Title:
Radius-valley age evolution around cool stars: calibrated gate, power grid, data, and code

Creators:
- Shaad, Zubayer Hasan (ORCID 0009-0002-2322-8553)
- Al Abrar, Md Mahdin (ORCID 0009-0009-9645-152X)

License: MIT. Upload type: Software. Access: Open.

The full description, keywords, and author list are already in the repository's
`.zenodo.json` and `CITATION.cff`; the GitHub-to-Zenodo integration reads these,
so keep them as the single source of truth.

## Files to release

Publish the whole repository at tag `v2.0.3`. It contains 292 tracked files, of
which the reproducibility bundle is `repro/v2/` (85 files verified by
`provenance/BUNDLE_SHA256.txt`).

Deliberately not distributed: LaTeX source (`main.tex`, `references.bib`),
AASTeX styles, and the compiled manuscript PDF. Data, code, results, figures,
and provenance records are included.

## Step-by-step

1. Confirm CI is green on `main` (checks frozen inputs, v2 bundle, all scripts,
   and the archived high-replicate audit).
2. In GitHub, go to Releases -> Draft a new release.
3. Choose **Create a new tag**: `v2.0.3`; target `main`.
4. Release title: `v2.0.3 - reproducible v2 bundle`.
5. In the description, state what is new: the v2 reproducibility bundle
   (`repro/v2/`), TIME Table age-ordering validation, the full power-versus-effect
   curve, Gaia production cuts, regenerated figures, and the new animation.
6. Publish the release. GitHub notifies Zenodo, which archives the tag and mints
   a new DOI (a new version of the same concept DOI).
7. The Zenodo version DOI is **10.5281/zenodo.24444755**. It is already in the
   README release paragraph and in the manuscript's Data Availability section.
8. If the manuscript or README changes again, commit and push normally. The
   Zenodo archive of v2.0.3 is immutable, so a later edit needs no new DOI.

## If you must upload manually instead

Create a clean zip of the repository at the `v2.0.3` tag (no `.git` directory),
then upload it to a new Zenodo record with the metadata above. Verify the upload
with the checksums in `repro/v2/provenance/BUNDLE_SHA256.txt` and
`data/checksums.sha256`.
