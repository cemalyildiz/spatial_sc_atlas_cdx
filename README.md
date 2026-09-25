# Human Spatial Atlas

A human-only spatial transcriptomics and single-cell / single-nucleus RNA-seq discovery catalogue.

**Website:** https://cemalyildiz.github.io/spatial_sc_atlas_cdx/

## Explore

- Full-text search, disease and subtype aliases, source, assay, measurement approach, tissue, access and release-date filters.
- Spatial, scRNA-seq and snRNA-seq remain separate modalities.
- Evidence-aware pairing: matched donor/sample subset, verified same study, same collection, or both methods reported but unverified.
- Source links, related accessions, publication links, shareable filters and CSV export.
- English interface, responsive layout, keyboard controls and accessible native dialogs.

## Public metadata and coverage

`collect.py` uses NCBI GEO E-utilities, the BioStudies ArrayExpress collection, CELLxGENE Discover and Europe PMC (including PubMed records). It downloads metadata only. No raw reads, expression matrices, or patient-level records are hosted here.

The catalogue is **not an exhaustive census**. Each broad GEO search has a 1,500-result cap; 37 disease/subtype query groups independently retrieve up to 35 hits each; each ArrayExpress query retrieves up to 200 hits. Europe PMC supplies up to 300 recent publication leads. Full query terms, observed totals, retrieved counts and timestamps are available in the site’s Sources & updates panel and `catalog.json`. CELLxGENE processes the collections endpoint, then retains human disease datasets. Broad GEO searches include historical spatial studies and single-cell studies released from 2024 onward; targeted queries also discover older studies. Source query caps can be changed in the workflow and collector arguments.

Metadata-derived disease, subtype and method labels are **mentions, not verified sample annotations**. Missing values remain unspecified. Single-cell-resolution imaging alone does not imply independent scRNA-seq. GEO sample entries may represent sections or regions rather than unique donors. A collection can include healthy controls. Mixed-species GEO and ArrayExpress records are excluded; only human datasets are selected from CELLxGENE.

`curation.json` stores reviewed overrides with evidence URLs. The breast cancer matched example applies to four overlapping tumors in a six-sample spatial cohort, not every tumor. Some accessions contain only one component of a multimodal study; additional data links and pairing evidence explain this. Identical accessions are deduplicated. Shared publication IDs create related-record links rather than merging distinct experiments. Source-specific records may still refer to the same biological cohort.

## Run and verify

Requires Python 3.9+ and Node.js 18+. No third-party runtime packages.

```sh
python3 -B collect.py --limit 1500 --ae-limit 200
python3 -B validate.py
python3 -B test_collect.py
node --test test_search.mjs
python3 -m http.server 8891 --bind 127.0.0.1
```

Open http://127.0.0.1:8891/. Serve over HTTP; opening `index.html` directly with `file://` will not load the JSON catalogue.

The disposable `cache_cdx` folder caches public API responses for 24 hours, is excluded from Git, and is never published. The website is packaged into `site_cdx`. GitHub requires its conventional `.github/workflows` path for automation.

## Updates and publishing

GitHub Pages uses the GitHub Actions source. The `Update human atlas and publish` workflow runs daily at 04:23 UTC (07:23 Europe/Istanbul), on pushes, and on manual dispatch. Pushes publish committed data without re-querying sources. Scheduled/manual refreshes collect, validate, commit `catalog.json`, and publish the same checkout. Source failures retain previous records and show a failure status; an empty aggregate never overwrites the catalogue. GitHub may delay scheduled jobs or disable inactive public-repository schedules; check the source timestamps and Actions history rather than assuming freshness.

Repository metadata is public. Original datasets retain their own licensing, controlled-access rules and citation requirements. Consult the [NCBI disclaimer](https://www.ncbi.nlm.nih.gov/About/disclaimer.html) and the linked source record before reuse.
