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

[PubChemPy](https://github.com/mcs07/PubChemPy) has been the standard way to script against PubChem in Python for years, and it covers a lot this package doesn't try to replace: 3D conformers, assays, substances, atoms and bonds. This package is narrower, focused on identifier resolution and cross-referencing, and adds a few things that matter specifically for that: resilience under PubChem's own rate limiting, a cache that keeps a notebook fast without going stale, and live search over PubChem's full corpus rather than a local index.

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
```

PubChem's own `xrefs/RegistryID` endpoint already carries the ChEMBL ID for a compound when one exists (verified live against aspirin, CID 2244 resolves to CHEMBL25). PubChemPy can reach the same endpoint through its low-level `request()`/`get()` functions; this package wraps it as a named, documented function. The batch versions accept a comma-separated CID list in one PUG REST request, the same way `resolve_many()` does (see below).

### Similarity and substructure search

```python
hits = pubchem.similar_compounds("CC(=O)OC1=CC=CC=C1C(=O)O", threshold=95, max_records=5)
# [Compound(cid=2244, title='Aspirin', ...), Compound(cid=4133, title='Methyl Salicylate', ...), ...]

pubchem.substructure_search("c1ccccc1", query_type="smiles")     # every compound containing a benzene ring
pubchem.substructure_search("[#6]1[#6][#6][#6][#6][#6]1", query_type="smarts")  # the SMARTS equivalent
```

`scigantic-chembl`'s `similar_compounds()`/`substructure_search()` precompute fingerprints once and search them locally: fast, but bounded to the roughly 1.68M ChEMBL compounds that carry a comparable measurement. This runs the search on PubChem's own servers, live, over the full ~120M-compound corpus, verified sub-second for a typical query and with no local fingerprint database to build or hold in memory. PubChemPy exposes the same PUG REST capability as a raw `searchtype="similarity"`/`"substructure"` parameter to its generic `get_compounds()`; this package gives it a named function, and keeps `query_type="smiles"` and `"smarts"` as separate, explicit paths rather than guessing between them, since they're genuinely different endpoints with different matching semantics (verified live: the same ring given as SMILES versus SMARTS returns overlapping but not identical results).

An expensive search can respond asynchronously, with PubChem handing back a job to poll rather than blocking the connection; handled transparently, using the same underlying protocol PubChemPy implements. Polling starts at 0.5s and doubles up to a 5s cap rather than a flat interval, so a fast job gets checked sooner and a slow one (a real, measured case took 30-60s) stops paying for a tight interval it never needed. Every live query tried during development resolved synchronously, even a maximally broad single-carbon substructure search, so the polling loop itself is verified with a scripted mock response sequence rather than left checked only against documentation.

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

## Batch resolution

```python
compounds = pubchem.resolve_many([2244, 2519, 1983])  # aspirin, caffeine, acetaminophen
```

PUG REST accepts a comma-separated CID list in a single request (verified live); this chunks at 200 CIDs per call rather than sending one unbounded URL for a long list. More than one chunk runs concurrently through a small thread pool, paced by the same rate limiter every request goes through, so a large list doesn't pay for each chunk's round trip in sequence.

## Command line

```console
$ scigantic-pubchem resolve aspirin
$ scigantic-pubchem chembl-id 2244
$ scigantic-pubchem xrefs 2244 --type RegistryID
$ scigantic-pubchem similar "CC(=O)OC1=CC=CC=C1C(=O)O" --threshold 95
$ scigantic-pubchem substructure c1ccccc1
```

## License

MIT-0. See [LICENSE](LICENSE).
