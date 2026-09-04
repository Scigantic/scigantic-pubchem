# Changelog

## 0.9.0

Fixes two client-side gaps found running a real batch workload (a few
hundred independent similarity searches, each followed by an assay-data
join), plus adds the batch primitive that workload needed:

- **A resolved search result is now cached.** `request_search()`
  (`similar_compounds()`/`substructure_search()`'s underlying async-search
  protocol) previously never cached anything, not even the final resolved
  body once polling ended, only the intermediate "still waiting" polls were
  meant to be excluded. So re-running the exact same search, e.g. resuming
  a crashed batch job partway through, always repeated the full live
  search (including any async poll) even though the `resolve_many()` call
  that follows it was already cached. The final body is now written to the
  same on-disk cache as any other lookup, under the search's own
  `(path, params)`; the polling calls themselves stay uncached, unchanged.
- **A malformed-but-real PubChem response no longer crashes the caller.**
  Verified live 2026-09-04, mid a large batch: a `concise`/`assaysummary`
  response can carry a depositor free-text field (an assay comment or
  description) with a raw, unescaped ASCII control character, which the
  standard JSON decoder (`strict=True`, the default for both `json.loads`
  and `requests.Response.json()`) rejects outright as "Invalid control
  character", regardless of retries, since the identical bytes come back
  every time. `request()` now falls back to `json.loads(text, strict=False)`,
  the stdlib's own documented escape hatch for exactly this, before giving
  up; a response that's corrupted for a real reason (truncation, garbage)
  still raises a clear `PubChemError` instead of a bare
  `json.JSONDecodeError`.
- **New `similar_compounds_many(smiles_list, ...)`**: `similar_compounds()`
  for many query structures at once, dispatched through a small thread
  pool (`max_workers`, default 5, deliberately lower than
  `resolve_many()`'s 8; see the function's docstring for why) rather than
  one query at a time. Returns a `dict` keyed by input SMILES so every
  query gets an entry, empty list included. Pairs with the existing
  `compound_assay_results_many()` for "find every close neighbor of a set
  of compounds, then pull assay data for all of them" in two calls instead
  of a hand-rolled loop. No batched form of `substructure_search()` yet;
  add one if a real workload needs it the way this one needed
  `similar_compounds_many()`.

## 0.8.0

Adds `dose_response()`/`download_dose_response()`: the raw, un-curated
per-concentration qHTS Data Table underneath `assay_results()`'s `concise`
summary. `concise`/`assaysummary` reduce a qHTS screen to outcome plus a
single fitted potency; the raw layer -- every tested concentration's raw
percent response, `Max_Response`, and curve-fit parameters, present even
for a compound with no fitted curve at all -- is a structurally different
PUG REST operation (the plain, non-`concise` `/CSV` Data Table), verified
live 2026-08-31 against AID 1851 (a 5-isoform CYP inhibition qHTS panel,
17,143 compounds): `concise` carries none of it.

- `dose_response(aid, sids=..., cids=...)` returns typed `RawAssayResult`
  rows (dataclass, with a `dose_response: tuple[DoseResponsePoint, ...]`
  field for the per-concentration readout) for a bounded set of up to 200
  identifiers in one call. Capped there rather than auto-chunked: this
  operation does not scale linearly the way `concise` does (measured live:
  250 SIDs under a second, 2000 SIDs took 94s against the same server),
  and PUG REST itself hard-caps it at 10,000 identifiers per request
  regardless.
- `download_dose_response(aid, dest, ...)` pulls a whole assay, chunked
  (200 SIDs by default) and resumable by default: progress is tracked in a
  `{dest}.progress.json` sidecar (assay identity, chunk size, and SID
  count, so a resume against different arguments is detected and
  discarded with a warning instead of silently misapplied) plus the exact
  byte offset `dest` was last known-good at, so an interruption mid-chunk
  truncates back to that offset and retries cleanly rather than risking a
  duplicated or corrupt file. Sequential, not parallelized the way
  `compound_assay_results_many()` is, given the demonstrated per-request
  cost on a live public NIH service.
- Only PubChem's reserved `PUBCHEM_*` columns, the standard NCGC/NCATS
  qHTS curve-fit vocabulary (`Potency`, `Curve_Description`, `Fit_*`,
  `Max_Response`), and the `Panel_*` columns (present on multi-target
  panel assays) are modelled as named `RawAssayResult` fields; everything
  else is depositor-specific and lands in `extra` rather than being
  hardcoded from one assay's shape. Note `Panel Name` carries the
  mnemonic (e.g. `"p450-cyp1a2"`) and `Panel Target` the protein
  accession (e.g. `"NP_000752.2"`) -- easy to mix up, verified live.
- New CLI commands `assay-dose-response`/`assay-dose-response-download`,
  matching the existing `assay-results`/`assay-download` pattern.

## 0.7.1

Documentation and CLI parity only. No change to any existing function's
behavior, arguments, or return values.

- Seven functions were exported in `__init__.py`/`__all__` and fully tested,
  but never reached the README or (for two of them) the CLI: `assay_cids`,
  `assay_sids`, `compound_assay_results_many` (`bioassay.py`); `pdb_structures_many`
  (`xrefs.py`); `inchi_keys`, `parent`, `substance_ids` (`mirror.py`). Added a
  README section/example for each.
- New CLI subcommands `assay-cids` and `assay-sids`, matching the existing
  `aids-for-compound`/`aids-for-target` pattern (a single identifier in, a
  flat list of ids out, one per line). The batch `_many` functions
  (`compound_assay_results_many`, `pdb_structures_many`) do not get CLI
  mirrors, matching the package's own existing convention: `resolve_many()`,
  `xrefs_many()`, and `chembl_ids_many()` were already CLI-less before this
  release, since a dict-of-lists keyed by CID has no natural one-line CLI
  shape the way a single-identifier lookup does. `mirror.py`'s `inchi_keys`,
  `parent`, and `substance_ids` are documented but also stay CLI-less, for
  the same reason their sibling `synonyms()` already was: only
  `mirror-identifiers`/`mirror-download` have ever had CLI mirrors in that
  module.
- `requests`, `duckdb`, and `pandas` now carry an explicit upper bound
  (`<3`, `<2`, `<4` respectively, checked live against each package's
  current PyPI major version), so a future breaking major release can't
  silently break an install the way an unbounded floor already broke this
  package's own `mcp` sibling once. Still real dependencies, no version
  pinned tighter than that.
- Repository topics set on GitHub (`pubchem`, `cheminformatics`, `bioassay`,
  `drug-discovery`, `chemistry`, `bioinformatics`); previously unset.

## 0.7.0 (2026-08-30, #8)

Adds `mirror.py`: reads `s3://scigantic-pubchem`, a new weekly-refreshed
parquet mirror of PubChem's own compound identifier/name/mass registry
(SMILES, titles, IUPAC names, InChI/InChIKeys, formula/mass, parent CIDs,
CID-to-SID provenance, synonyms), for CID-keyed bulk/offline work that
doesn't need a live PUG REST round trip.

- `identifiers()`, `inchi_keys()`, `synonyms()`, `parent()`, `substance_ids()`:
  live DuckDB `httpfs` reads against the mirror's per-table parquet files,
  reusing `bridge.py`'s shared connection.
- `download_mirror()` pulls the parquet files down once for real local/offline
  use (query them afterward with your own DuckDB/pandas session, no
  per-query network round trip). Skips a file already present at the
  destination with a matching size.
- Deliberately does not add name/InChIKey/SMILES -> CID resolution against
  the mirror: every table is CID-sorted, so a CID filter gets real
  row-group pruning from a remote parquet read, but any other filter would
  mean a full remote column scan. Use the existing live `resolve()`/`xrefs()`
  for that instead.
- Corrected the package's own top-level docstring and README opening line,
  which had claimed "No mirror, no download" since the 0.1.0 release; no
  longer true as of this version.
- New CLI: `mirror-identifiers`, `mirror-download`.

## 0.6.2 (2026-08-29)

Metadata only. Dropped the "...PubChemPy doesn't have" phrasing from the
package description on PyPI; read as a dig at that project rather than a
plain statement of what this one does.

## 0.6.1 (2026-08-29, #7)

Adds `pdb_structures()`/`pdb_structures_many()`: cross-references into
NCBI's Molecular Modeling Database (MMDB), PUG REST's closest reachable
proxy for "structures with this compound bound as a ligand."

- Verified live that PUG REST's `xrefs` endpoint has no type that returns a
  4-character PDB ID directly: `MMDBID` is valid, `PDBID` 400s with "Invalid
  xrefs type". Mapping an MMDB ID to its PDB ID needs a separate NCBI
  service outside PUG REST, out of scope for a package that only ever calls
  `pubchem.ncbi.nlm.nih.gov/rest/pug`; each MMDB ID is still a real, directly
  usable pointer to a deposited structure (imatinib/CID 5291 resolves to 27
  of them).
- This surfaced and fixed a real documentation bug in the already-shipped
  `xrefs_many()`: its docstring claimed a CID with zero xrefs is "simply
  absent" from the result, but that's only true when every CID in a chunk
  misses at once (the whole chunk 404s). A chunk with at least one match
  returns an explicit empty-list entry for every CID in it that missed, so
  presence in the result dict actually depends on chunk-mates, not just on
  the CID itself. `pdb_structures_many()` sidesteps this by looping over the
  input CIDs explicitly, the same pattern `chembl_ids_many()` already used.
- New CLI: `pdb-structures`.

## 0.6.0 (2026-08-29, #6)

Adds `tox21.py`: the 12-endpoint Tox21 Data Challenge nuclear-receptor and
stress-response panel (NR-AR, NR-AR-LBD, NR-AhR, NR-Aromatase, NR-ER,
NR-ER-LBD, NR-PPAR-gamma, SR-ARE, SR-ATAD5, SR-HSE, SR-MMP, SR-p53), the
same 12 columns MoleculeNet/DeepChem's `tox21.csv` exposes for ML
benchmarking.

- `TOX21_ENDPOINTS` maps each endpoint name to its PubChem AID. Each AID
  was verified live three independent ways: cross-checked against Huang et
  al. 2016 (Frontiers in Environmental Science 3:85, Table 1), confirmed
  each AID's own PubChem `summary` name matches the claimed target/pathway,
  and confirmed all 12 share the same 8,099- or 7,329-CID Tox21 10K library
  count.
- `tox21_results()` is a thin wrapper over the existing `assay_results()`
  (one batched request for the whole panel). `tox21_matrix()` does the
  assembly step raw PubChem access doesn't: a wide CID x endpoint matrix,
  consolidating each compound's replicate rows into a single active(1) /
  inactive(0) / untested(`None`) label per endpoint.
- Documented as a live reconstruction of the panel from PubChem's own qHTS
  data, not a byte-for-bit copy of NCATS's original challenge SDF (which
  applied its own compound-level SMILES canonicalization and train/test
  split).
- New CLI: `tox21-results`, `tox21-matrix`.

## 0.5.0 (2026-08-29, #5)

The BioAssay and gene/protein domains, PubChemPy has no equivalent for
either (verified by reading its source directly: `concise`, `assaysummary`,
`summary` never appear in it, and it has no `Gene`/`Protein` class at all).
This version also carries three real fixes merged as PRs #1-#3 the same
day, whose own version bumps (0.3.2 through 0.3.4) were never independently
tagged or published to PyPI, so everything below reached users together as
0.5.0.

**BioAssay** (`bioassay.py`): `assay_summary()` (flat overview plus
active/inactive/total counts), `assay_results()`/`compound_assay_results()`
(row-level bioactivity table, both directions, batched across AIDs/CIDs in
one PUG REST request), `assay_cids()`/`assay_sids()`, `aids_for_compound()`/
`aids_for_target()`, and `download_assay_results()` (streams a `concise`
table straight to disk instead of buffering it in memory; a modern qHTS
screen can run to hundreds of thousands of rows).

**Gene and protein** (`gene_protein.py`): `gene_info()`/`protein_info()`
give the flat overview PubChem's own `summary` operations return.
`gene_assay_results()`/`protein_assay_results()` are a third direction of
the bioactivity table, keyed by the target itself rather than requiring
`aids_for_target()` plus a per-AID fetch first.

Shared table-parsing logic (column-name-keyed, since assay/compound/gene/
protein `concise` tables return overlapping but not identical column sets,
verified live) was factored into `_tables.py`, used by both modules.

New CLI subcommands: `assay-summary`, `assay-results`, `compound-assay-results`,
`aids-for-compound`, `aids-for-target`, `assay-download`, `gene-info`,
`protein-info`, `gene-assay-results`, `protein-assay-results`.

Also bundled in this release, from the same day's earlier PRs:

- (#1) A token-bucket rate limiter paced to PubChem's documented 5
  requests/second, acquired before every request attempt including
  retries, so a burst is paced proactively rather than relying on
  PubChem's `X-Throttling-Control` header to say "slow down" after the
  fact. `resolve_many()`'s chunk loop was also parallelized through a
  bounded thread pool (measured: 30.5s for 50 CIDs via one batched request
  path versus 90.7s doing them one at a time, 3.0x).
- (#2) Async search polling (for slow similarity/substructure searches)
  now backs off adaptively, starting at 0.5s and doubling to a 5s cap,
  instead of a flat 2s interval that wasted round trips on fast jobs and
  underserved slow ones (a real measured case took 30-60s).
- (#3) `xrefs_many()`/`chembl_ids_many()`: cross-referencing N compounds to
  ChEMBL used to mean N round trips. PUG REST's `xrefs` endpoint accepts a
  comma-separated CID list the same way `resolve_many()`'s property
  endpoint does, verified live at 200 CIDs in one request. `bridge.py`'s
  DuckDB connection is now created once and reused via `cursor()` per call,
  instead of reopened (re-running `INSTALL`/`LOAD httpfs` and
  re-registering the S3 secret) on every `chembl_context()`/
  `bindingdb_measurements()` call.

## 0.3.1 (2026-08-28)

Fixed a real concurrency bug: `cache.put()`'s temp filename was derived
only from the cache key, so two threads racing to fill the same key shared
one temp path, and whichever thread's `os.replace()` ran second raised
`FileNotFoundError` because the first had already consumed it. Found via a
32-thread stress run against the published package. Fixed by making the
temp filename unique per write.

## 0.3.0 (2026-08-27)

Stress/concurrency test coverage: multi-chunk `resolve_many()` (250 CIDs
forcing a real 200+50 split), concurrent resolution from 8 threads (which
surfaced and fixed a real lazy-init race in `_get_session()`, two threads
could each see the shared session as unset and construct their own; fixed
with double-checked locking), and the async search-polling protocol
verified deterministically with a scripted mock response sequence rather
than checked only against documentation (every live query tried during
development had resolved synchronously).

Not independently tagged or published; superseded within a day by 0.3.1.

## Earlier, unpublished development versions

`0.1.0` (first commit: live PUG REST queries with retry/backoff and an
always-on cache, no mirror) and `0.2.0` (live similarity/substructure
search over PubChem's full ~120M-compound corpus; an expiring cache,
default 30-day TTL; and a real bug fix, where async search polling could
cache an in-progress `{"Waiting": ...}` response under the search's own
key) were real version bumps in this repository's history but were never
tagged or published to PyPI. The earliest version on PyPI is `0.3.0`.
