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

Query PubChem live via PUG REST -- no mirror, no download, no local database.

```python
import scigantic_pubchem as pubchem

aspirin = pubchem.resolve("aspirin")
print(aspirin.cid, aspirin.smiles, aspirin.inchi_key)
```

## Installation

```console
$ pip install scigantic-pubchem
```

## Compared to PubChemPy

[PubChemPy](https://github.com/mcs07/PubChemPy) is the standard way to script against PubChem in Python, and this package doesn't try to replace its full surface (3D conformers, assays, substances, atoms/bonds). What it does differently, verified against PubChemPy 1.0.5's actual source rather than assumed:

**Rate-limit awareness.** Every PUG REST response carries an `X-Throttling-Control` header reporting live status across three dimensions (request count, request time, service load). PubChemPy's source never reads this header -- on an error it raises immediately, no retry. This package reads it on every call, backs off proactively before hitting a hard limit, and retries 429/5xx with backoff instead of failing a notebook cell over a transient blip.

**Caching, on by default.** A lookup you've already made is never re-fetched -- cached to `~/.cache/scigantic-pubchem` (override with `enable_cache(cache_dir=...)` or `SCIGANTIC_PUBCHEM_CACHE`). This is a deliberate difference from `scigantic-chembl`/`scigantic-bindingdb`, whose caching defaults off: those read a public S3 mirror with no meaningful rate limit, so caching there is pure convenience. This package calls a rate-limited live API for every lookup, so re-fetching the same identifier in a loop is both slow and the exact thing the throttling-aware client otherwise works to avoid.

```python
pubchem.disable_cache()  # if you want every call to hit the network fresh
```

**Live cross-references, not just structures.**

```python
pubchem.chembl_id(2244)   # 'CHEMBL25' -- read live from PubChem's own xrefs, not a static table
```

PubChem's own `xrefs/RegistryID` endpoint already carries the ChEMBL ID for a compound when one exists (verified live against aspirin, CID 2244 -> CHEMBL25) -- reachable through PubChemPy's low-level `request()`/`get()` functions, but not wrapped as a named, documented convenience there.

## Live bridge into scigantic-chembl and scigantic-bindingdb

```python
pubchem.chembl_context(2244)
# {'molregno': ..., 'chembl_id': 'CHEMBL25', 'pref_name': 'ASPIRIN', 'max_phase': 4}

pubchem.bindingdb_measurements(2244)
# DataFrame of every BindingDB measurement recorded against this CID
```

Both are on-demand DuckDB queries against the public [scigantic-chembl](https://github.com/Scigantic/scigantic-chembl) and [scigantic-bindingdb](https://github.com/Scigantic/scigantic-bindingdb) mirrors -- `chembl_context` resolves the ChEMBL ID live (see above) then looks up its assay context; `bindingdb_measurements` filters BindingDB's own `pubchem_cid` column directly, since BindingDB already carries that mapping natively. Neither needs a precomputed bridge table, and neither goes stale the way a static one would. Needs `duckdb`:

```console
$ pip install "scigantic-pubchem[bridge]"
```

## Batch resolution

```python
compounds = pubchem.resolve_many([2244, 2519, 1983])  # aspirin, caffeine, acetaminophen
```

PUG REST accepts a comma-separated CID list in a single request (verified live); this chunks at 200 CIDs per call rather than sending one unbounded URL for a long list.

## Command line

```console
$ scigantic-pubchem resolve aspirin
$ scigantic-pubchem chembl-id 2244
$ scigantic-pubchem xrefs 2244 --type RegistryID
```

## License

MIT-0. See [LICENSE](LICENSE).
