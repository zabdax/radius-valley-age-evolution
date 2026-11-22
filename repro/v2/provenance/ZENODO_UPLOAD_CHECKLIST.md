# Zenodo release guide for v2.1.0

This repository is connected to Zenodo through GitHub. A new Zenodo version is
created automatically when you publish a new GitHub release, so the correct way
to mint the v2.1.0 DOI is to publish a release/tag here rather than uploading
files by hand.

## Status

- Current archived release: **v2.0.1**, DOI 10.5281/zenodo.23007310
- Next release to publish: **v2.1.0**
- No v2.1.0 DOI exists yet. Do not insert a guessed DOI into the manuscript.

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

Publish the whole repository at tag `v2.1.0`. It contains 292 tracked files, of
which the reproducibility bundle is `repro/v2/` (85 files verified by
`provenance/BUNDLE_SHA256.txt`).

Deliberately not distributed: LaTeX source (`main.tex`, `references.bib`),
AASTeX styles, and the compiled manuscript PDF. Data, code, results, figures,
and provenance records are included.

## Step-by-step

1. Confirm CI is green on `main` (checks frozen inputs, v2 bundle, all scripts,
   and the archived high-replicate audit).
2. In GitHub, go to Releases -> Draft a new release.
3. Choose **Create a new tag**: `v2.1.0`; target `main`.
4. Release title: `v2.1.0 - reproducible v2 bundle`.
5. In the description, state what is new: the v2 reproducibility bundle
   (`repro/v2/`), TIME Table age-ordering validation, the full power-versus-effect
   curve, Gaia production cuts, regenerated figures, and the new animation.
6. Publish the release. GitHub notifies Zenodo, which archives the tag and mints
   a new DOI (a new version of the same concept DOI).
7. Wait until the Zenodo record shows version v2.1.0 and copy its **version DOI**
   (not the concept DOI).
8. Update the README release paragraph with that DOI, and optionally paste it
   into the Data Availability section of the manuscript before submission.
9. Commit and push the README change. The Zenodo archive of v2.1.0 is immutable;
   a README fix afterwards is fine and does not require a new DOI.

## If you must upload manually instead

Create a clean zip of the repository at the `v2.1.0` tag (no `.git` directory),
then upload it to a new Zenodo record with the metadata above. Verify the upload
with the checksums in `repro/v2/provenance/BUNDLE_SHA256.txt` and
`data/checksums.sha256`.
