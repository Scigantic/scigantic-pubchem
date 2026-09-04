<h1 align="center">scigantic-pubchem</h1>

<p align="center">
    <a href="https://github.com/Scigantic/scigantic-pubchem/actions/workflows/ci.yml">
        <img alt="CI" src="https://github.com/Scigantic/scigantic-pubchem/actions/workflows/ci.yml/badge.svg" /></a>
    <a href="https://pypi.org/project/scigantic-pubchem/">
        <img alt="PyPI" src="https://img.shields.io/pypi/v/scigantic-pubchem" /></a>
    <a href="https://pypi.org/project/scigantic-pubchem/">
        <img alt="PyPI - Python Version" src="https://img.shields.io/pypi/pyversions/scigantic-pubchem" /></a>
    <a href="https://github.com/Scigantic/scigantic-pubchem/blob/main/LICENSE">
        <img alt="License" src="https://img.shields.io/github/license/Scigantic/scigantic-pubchem" /></a>
</p>

Query PubChem live via PUG REST. No mirror, no download, no local database.

```python
import scigantic_pubchem as pubchem

aspirin = pubchem.resolve("aspirin")
print(aspirin.cid, aspirin.smiles, aspirin.inchi_key)
```

## Installation

```console
$ pip install scigantic-pubchem
```

## Why this exists

[PubChemPy](https://github.com/mcs07/PubChemPy) has been the standard way to script against PubChem in Python for years, and it covers a lot this package doesn't try to replace: 3D conformers, substances, atoms and bonds. This package is narrower, focused on identifier resolution, cross-referencing, and bioassay/gene/protein data, and adds a few things that matter specifically for that: resilience under PubChem's own rate limiting, a cache that keeps a notebook fast without going stale, live search over PubChem's full corpus rather than a local index, and (see BioAssay and Gene and protein below) whole domains PubChemPy doesn't reach at all.

Measured, not asserted, on 2026-08-27:

| | measured |
|---|---|
| Cached lookup vs a live PUG REST round trip | under 1ms vs 250ms-2.4s in testing, roughly 1,000-10,000x depending on the identifier |
| Batch resolution (50 CIDs) | 30.5s in one request via `resolve_many()`, versus 90.7s across 50 separate requests done one at a time (3.0x) |
| Live 503 during an 8-thread concurrent test | recovered automatically via retry, no failed lookups |
| Similarity search corpus size | PubChem's full ~120M compounds, live, versus a precomputed local index bounded to whatever was indexed ahead of time |

The rate-limit handling works at two levels. Every request first goes through a token bucket paced to PubChem's documented 5 requests/second, so a burst (parallel batch chunks, a caller's own thread pool) is paced up front rather than relying on PubChem to say "slow down" after the fact. On top of that, every response is checked against PUG REST's own `X-Throttling-Control` header, which reports live status across three dimensions (request count, request time, service load), and backs off further if it's elevated. PubChemPy's source doesn't read this header at all, so it has no way to back off before hitting a hard limit, and raises immediately on an error with no retry. That's a reasonable design for a general-purpose client; this package leans the other way on purpose, since backing off proactively and retrying transient failures matters more when a notebook is doing dozens of lookups in a loop.

### Caching

On by default, and it expires. A lookup already made is never re-fetched, cached to `~/.cache/scigantic-pubchem` (override with `enable_cache(cache_dir=...)` or `SCIGANTIC_PUBCHEM_CACHE`). This is a deliberate difference from `scigantic-chembl` and `scigantic-bindingdb`, whose caching defaults off: those read a public S3 mirror with no meaningful rate limit, so caching there is pure convenience. This package calls a rate-limited live API for every lookup, so re-fetching the same identifier in a loop is both slow and the exact thing the throttling-aware client otherwise works to avoid. Entries expire after 30 days by default, so the cache can't quietly turn into a stale snapshot:

```python
pubchem.disable_cache()               # every call hits the network fresh
pubchem.enable_cache(ttl_days=7)      # shorter freshness window
pubchem.enable_cache(ttl_days=None)   # never expire
```

### Cross-references

```python
pubchem.chembl_id(2244)   # 'CHEMBL25', read live from PubChem's own xrefs, not a static table

pubchem.xrefs_many([2244, 3672, 2519])       # {2244: [...], 3672: [...], 2519: [...]}, chunked and parallelized
pubchem.chembl_ids_many([2244, 3672, 2519])  # {2244: 'CHEMBL25', 3672: 'CHEMBL521', 2519: 'CHEMBL113'}

pubchem.pdb_structures(5291)   # [110242, 131625, ...], 27 deposited structures with imatinib bound as a ligand
pubchem.pdb_structures_many([5291, 2244])   # {5291: [110242, ...], 2244: [...]}, every requested CID gets an entry
```

PubChem's own `xrefs/RegistryID` endpoint already carries the ChEMBL ID for a compound when one exists (verified live against aspirin, CID 2244 resolves to CHEMBL25). PubChemPy can reach the same endpoint through its low-level `request()`/`get()` functions; this package wraps it as a named, documented function. The batch versions accept a comma-separated CID list in one PUG REST request, the same way `resolve_many()` does (see below). `pdb_structures_many()` is the batch form of `pdb_structures()`, looping over the input CIDs explicitly so every one gets an entry (an empty list if it has no structures on file) regardless of `xrefs_many()`'s own chunk-dependent presence (see its docstring).

`pdb_structures()` reads a different xrefs type, `MMDBID`: NCBI's Molecular Modeling Database, its own mirror of PDB structure data. Not the 4-character PDB ID itself -- verified live 2026-08-29 that PUG REST has no xrefs type that returns one directly (`MMDBID` is valid, `PDBID` 400s), and mapping an MMDB ID to its PDB ID needs a separate NCBI service outside PUG REST, out of scope for a package that only ever calls `pubchem.ncbi.nlm.nih.gov/rest/pug`. Each ID is still a real, usable pointer to a deposited structure, viewable at `ncbi.nlm.nih.gov/Structure/mmdb/mmdbsrv.cgi?uid={id}`.

### Similarity and substructure search

```python
hits = pubchem.similar_compounds("CC(=O)OC1=CC=CC=C1C(=O)O", threshold=95, max_records=5)
# [Compound(cid=2244, title='Aspirin', ...), Compound(cid=4133, title='Methyl Salicylate', ...), ...]

pubchem.substructure_search("c1ccccc1", query_type="smiles")     # every compound containing a benzene ring
pubchem.substructure_search("[#6]1[#6][#6][#6][#6][#6]1", query_type="smarts")  # the SMARTS equivalent

pubchem.similar_compounds_many([smiles_1, smiles_2, ...])  # {smiles_1: [...], smiles_2: [...], ...}
```

`scigantic-chembl`'s `similar_compounds()`/`substructure_search()` precompute fingerprints once and search them locally: fast, but bounded to the roughly 1.68M ChEMBL compounds that carry a comparable measurement. This runs the search on PubChem's own servers, live, over the full ~120M-compound corpus, verified sub-second for a typical query and with no local fingerprint database to build or hold in memory. PubChemPy exposes the same PUG REST capability as a raw `searchtype="similarity"`/`"substructure"` parameter to its generic `get_compounds()`; this package gives it a named function, and keeps `query_type="smiles"` and `"smarts"` as separate, explicit paths rather than guessing between them, since they're genuinely different endpoints with different matching semantics (verified live: the same ring given as SMILES versus SMARTS returns overlapping but not identical results).

An expensive search can respond asynchronously, with PubChem handing back a job to poll rather than blocking the connection; handled transparently, using the same underlying protocol PubChemPy implements. Polling starts at 0.5s and doubles up to a 5s cap rather than a flat interval, so a fast job gets checked sooner and a slow one (a real, measured case took 30-60s) stops paying for a tight interval it never needed. Every live query tried during development resolved synchronously, even a maximally broad single-carbon substructure search, so the polling loop itself is verified with a scripted mock response sequence rather than left checked only against documentation. Once a search resolves, the result is cached under its own query the same as any other lookup, so re-running the same search (resuming a batch job after a crash, say) doesn't re-pay the live search, poll included.

`similar_compounds_many()` runs several independent similarity searches concurrently (there's no PubChem-side batched form of this endpoint, unlike `resolve_many()`'s comma-separated CID list), returning a dict keyed by the input SMILES so nothing is silently dropped. It's the composable half of "find every close neighbor of a set of compounds, then pull assay data for all of them": pair it with `compound_assay_results_many()` on the returned CIDs.

### BioAssay

```python
summary = pubchem.assay_summary(1)
# AssaySummary(aid=1, name='NCI human tumor cell line growth inhibition assay...',
#              cid_active=3370, cid_inactive=52324, cid_total=55532, ...)

results = pubchem.assay_results(1)          # every (SID, CID, outcome) row PubChem has for this assay
pubchem.compound_assay_results(2244)        # every assay result recorded for aspirin, the reverse direction
pubchem.compound_assay_results_many([2244, 3672])  # {2244: [...], 3672: [...]}, chunked and parallelized

pubchem.assay_cids(1)                       # every CID tested in this assay
pubchem.assay_sids(1)                       # every SID (pre-standardization substance ID) tested in this assay

pubchem.aids_for_compound(2244)             # every AID that tested aspirin
pubchem.aids_for_target("EGFR")             # every AID run against a gene target
```

`compound_assay_results_many()` is the batch form of `compound_assay_results()`, chunked at 200 CIDs per request and parallelized the same way `resolve_many()`/`xrefs_many()` are; a CID with no assay results at all is simply absent from the returned dict. `assay_cids()`/`assay_sids()` read PUG REST's own `cids`/`sids` operations directly: the same CID/SID columns `assay_results()`'s `concise` table already carries per row, as a plain list when all you need is "which compounds/substances did this assay touch," with no per-row outcome data to parse.

PubChemPy's `Assay`/`get_assays()` only reach PUG REST's `description` operation, the raw, deeply nested record built to round-trip a depositor's original submission (protocol text, full result-column schema, revision history), not to be read programmatically. Verified 2026-08-29 by reading PubChemPy's source directly: the strings `concise`, `assaysummary`, and `summary` never appear in it, under the assay domain or otherwise, so it has no path to the tabular bioactivity data (AID/SID/CID/Activity Outcome/...) most callers actually want. This package wraps those operations instead. `assay_summary()` gives the flat overview (name, description, target, active/inactive/total counts) `description` buries in that nested record; `assay_results()`/`compound_assay_results()` give the row-level bioactivity table, in both directions, using PUG REST's `concise` and `assaysummary` operations, which are PubChem's own purpose-built compact formats for exactly this.

A large assay's result table does not fit comfortably in memory as a Python list. Verified 2026-08-29: AID 3, a DTP/NCI screen from the 1990s, is 54,003 rows and about 9MB as `concise` CSV; a modern qHTS screen can run to hundreds of thousands of rows. `download_assay_results()` streams a `concise` table straight to a file, chunk by chunk, rather than buffering the whole response in memory first the way `assay_results()` does:

```python
pubchem.download_assay_results(1259416, "aid1259416.csv")             # CSV, PubChem's own bulk format for this
pubchem.download_assay_results([1, 3], "combined.csv")                 # multiple AIDs in one request, one file
pubchem.download_assay_results(1259416, "aid1259416.json", fmt="json")
```

`concise`/`assaysummary` reduce a qHTS screen down to outcome plus a single fitted potency. The raw layer underneath that -- every tested concentration's raw percent response, including for a compound PubChem marked Inactive with no fitted curve at all -- is a different PUG REST operation entirely (the plain, non-`concise` `/CSV` Data Table), verified live 2026-08-31 against AID 1851 (a 5-isoform CYP inhibition qHTS panel, 17,143 compounds):

```python
results = pubchem.dose_response(1851, sids=[842238])   # or cids=[...]
# [RawAssayResult(aid=1851, sid=842238, panel_name='p450-cyp1a2', activity_outcome='Inactive',
#                  max_response=-11.9, potency_um=None,
#                  dose_response=(DoseResponsePoint(concentration_um=0.0007, response_percent=-11.9), ...), ...), ...]

pubchem.download_dose_response(1259416, "aid1259416_raw.csv")   # a whole assay, chunked and resumable
```

`dose_response()` is capped at 200 identifiers per call (raise past that with a pointer to `download_dose_response()`): this operation does not scale linearly the way `concise` does -- measured live against AID 1851, 250 SIDs returned in under a second, 2000 SIDs took 94s against the same server -- and PUG REST itself hard-caps it at 10,000 identifiers per request regardless ("Assay record retrieval is limited to 10000 SIDs"). `download_dose_response()` chunks conservatively (200 SIDs by default) and writes to disk as it goes, resumable by default: an interruption partway through a large panel assay (realistically many minutes even chunked) leaves a `{dest}.progress.json` sidecar, and re-running the same call picks up after the last confirmed chunk rather than starting over or risking a duplicated/corrupt file.

### Gene and protein

```python
pubchem.gene_info("EGFR")
# GeneInfo(gene_id=1956, symbol='EGFR', name='epidermal growth factor receptor',
#          taxonomy='Homo sapiens (human)', synonyms=['ERBB', 'HER1', ...], ...)

pubchem.protein_info("P00533")
# ProteinInfo(accession='P00533', name='Epidermal growth factor receptor', ...)

pubchem.gene_assay_results("EGFR")           # every bioactivity row recorded against this gene, across every assay
pubchem.protein_assay_results("P00533")      # same, keyed by protein accession instead
```

A third direction alongside `assay_results()`/`compound_assay_results()` above, this time keyed by the target itself. `AssayResult` already carries `target_accession`/`target_gene_id` from the bioassay tables; these give that field somewhere to resolve to, and a bioactivity table read directly by target rather than requiring `aids_for_target()` plus a per-AID fetch first. PubChemPy has nothing here at all: verified 2026-08-29 by reading its source directly, it has no `Gene`/`Protein` class and none of `genesymbol`, `geneid`, or `ProteinAccession` appear in it anywhere. `gene_info()` takes `namespace="genesymbol"` (default) or `namespace="geneid"` for PubChem's numeric Entrez ID; `protein_info()` takes a protein accession (e.g. UniProt's `P00533`).

### Tox21

```python
pubchem.tox21_results(["NR-AhR"])          # raw bioactivity rows for one endpoint
pubchem.tox21_results()                    # all 12 endpoints, one batched request

matrix = pubchem.tox21_matrix()
# {2244: {'NR-AR': 0, 'NR-AR-LBD': 0, 'NR-AhR': 0, ..., 'SR-p53': 0}, ...}
```

The 12 nuclear-receptor and stress-response qHTS assays from the Tox21 Data Challenge (NR-AR, NR-AR-LBD, NR-AhR, NR-Aromatase, NR-ER, NR-ER-LBD, NR-PPAR-gamma, SR-ARE, SR-ATAD5, SR-HSE, SR-MMP, SR-p53) -- the same 12 columns MoleculeNet/DeepChem's `tox21.csv` exposes. `tox21_results()` is a thin wrapper over `assay_results()` with the 12 endpoint names mapped to their PubChem AIDs (each AID verified live against its own PubChem `summary` name and cross-checked against Huang et al. 2016, Table 1). `tox21_matrix()` does the assembly step raw PubChem access doesn't: one row per CID, one column per endpoint, consolidating a compound's replicate rows into a single 1 (active) / 0 (inactive) / `None` (untested or unresolved) label -- Tox21 is a canonical multi-task benchmark with real missing labels, and this is a live reconstruction of the panel, not a byte-for-bit copy of NCATS's original challenge SDF (which applied its own compound-level SMILES canonicalization and train/test split).

## Thread safety

Safe to call from multiple threads, a plausible real pattern: resolving a list of names via a `ThreadPoolExecutor`, say. The shared HTTP session is created once behind a lock rather than raced into existence by whichever thread gets there first. `enable_cache()`/`disable_cache()` are not synchronized against concurrent reads, the same way mutating `os.environ` isn't: call them once at the start of a script, not from multiple threads at once.

## Live bridge into scigantic-chembl and scigantic-bindingdb

```python
pubchem.chembl_context(2244)
# {'molregno': ..., 'chembl_id': 'CHEMBL25', 'pref_name': 'ASPIRIN', 'max_phase': 4}

pubchem.bindingdb_measurements(2244)
# DataFrame of every BindingDB measurement recorded against this CID
```

Both are on-demand DuckDB queries against the public [scigantic-chembl](https://github.com/Scigantic/scigantic-chembl) and [scigantic-bindingdb](https://github.com/Scigantic/scigantic-bindingdb) mirrors. `chembl_context` resolves the ChEMBL ID live (see above), then looks up its assay context; `bindingdb_measurements` filters BindingDB's own `pubchem_cid` column directly, since BindingDB already carries that mapping natively. Neither needs a precomputed bridge table, and neither goes stale the way a static one would. Needs `duckdb`:

```console
$ pip install "scigantic-pubchem[bridge]"
```

## Local mirror

Everything above is a live PUG REST call. For CID-keyed bulk or offline work, this package also reads [s3://scigantic-pubchem](https://github.com/Scigantic/scigantic-pubchem), a weekly-refreshed local parquet mirror of PubChem's compound identifier/name/mass registry (~15.5 GB: SMILES, titles, IUPAC names, InChI/InChIKeys, formula/mass, parent CIDs, CID-to-SID provenance, and synonyms):

```python
pubchem.identifiers(2244)
# {'cid': 2244, 'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O', 'title': 'Aspirin', ...}

pubchem.synonyms(2244)
# ['Aspirin', 'Acetylsalicylic acid', '2-Acetoxybenzoic acid', ...]

pubchem.inchi_keys(2244)
# DataFrame: one (inchi, inchikey) row per protonation state PubChem has on file for this CID

pubchem.parent(2244)          # 2244; a CID may be its own parent, or have none (returns None)

pubchem.substance_ids(2244)
# DataFrame of (sid, link_type) rows: every SID this CID derives from, per PubChem's own README-Extras
```

`identifiers()` already returns the first InChIKey it finds; `inchi_keys()` is for when a CID has more than one (one per protonation state) and you need all of them, not just the first. `substance_ids()`'s `link_type` is 1 for the standardized form of the deposited substance, 2 for a component of it.

Every mirror table is exported in ascending CID order, so a CID-filtered query gets real row-group pruning from a remote parquet read; a name/InChIKey/SMILES query would not, so use `resolve()`/`xrefs()` above for those instead of this mirror.

For real bulk/offline work, `download_mirror()` pulls the parquet files down once so you can query them locally with your own DuckDB or pandas session, no network round trip per query:

```python
pubchem.download_mirror("./pubchem_mirror")  # all 8 tables, ~15.5 GB
pubchem.download_mirror("./pubchem_mirror", tables=["smiles", "titles"])  # just these, ~3.4 GB
```

Needs `duckdb` for the live lookup functions (`identifiers`, `inchi_keys`, `synonyms`, `parent`, `substance_ids`), same extra as the ChEMBL/BindingDB bridge above; `download_mirror` only needs `requests`, already a base dependency.

## Batch resolution

```python
compounds = pubchem.resolve_many([2244, 2519, 1983])  # aspirin, caffeine, acetaminophen
```

PUG REST accepts a comma-separated CID list in a single request (verified live); this chunks at 200 CIDs per call rather than sending one unbounded URL for a long list. More than one chunk runs concurrently through a small thread pool, paced by the same rate limiter every request goes through, so a large list doesn't pay for each chunk's round trip in sequence.

## Command line

```console
$ scigantic-pubchem resolve aspirin
$ scigantic-pubchem chembl-id 2244
$ scigantic-pubchem pdb-structures 5291
$ scigantic-pubchem xrefs 2244 --type RegistryID
$ scigantic-pubchem similar "CC(=O)OC1=CC=CC=C1C(=O)O" --threshold 95
$ scigantic-pubchem substructure c1ccccc1
$ scigantic-pubchem assay-summary 1
$ scigantic-pubchem assay-results 1 3
$ scigantic-pubchem compound-assay-results 2244
$ scigantic-pubchem assay-cids 1
$ scigantic-pubchem assay-sids 1
$ scigantic-pubchem aids-for-compound 2244
$ scigantic-pubchem aids-for-target EGFR
$ scigantic-pubchem assay-download 1259416 aid1259416.csv
$ scigantic-pubchem assay-dose-response 1851 --sid 842238
$ scigantic-pubchem assay-dose-response-download 1259416 aid1259416_raw.csv
$ scigantic-pubchem gene-info EGFR
$ scigantic-pubchem protein-info P00533
$ scigantic-pubchem gene-assay-results EGFR
$ scigantic-pubchem protein-assay-results P00533
$ scigantic-pubchem tox21-results NR-AhR
$ scigantic-pubchem tox21-matrix
$ scigantic-pubchem mirror-identifiers 2244
$ scigantic-pubchem mirror-download ./pubchem_mirror --table smiles --table titles
```

## License

MIT-0. See [LICENSE](LICENSE).
