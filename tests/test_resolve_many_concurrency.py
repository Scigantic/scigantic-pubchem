"""resolve_many()'s multi-chunk dispatch: mocked, so the parallelism claim
is verified deterministically by wall-clock elapsed time, rather than
inferred from a live call whose real latency can't distinguish sequential
from concurrent execution either way.
"""

import sys
import time
from unittest import mock

import scigantic_pubchem.resolve  # noqa: F401 -- populates sys.modules below

# Not "from scigantic_pubchem import resolve", and not
# "import scigantic_pubchem.resolve as _resolve_module" either:
# scigantic_pubchem/__init__.py does "from .resolve import resolve", which
# rebinds the *package's* `resolve` attribute to that function. Both of
# those forms resolve through attribute access on the package and hit that
# same shadowed name. sys.modules is keyed by the dotted path and is
# unaffected by it.
_resolve_module = sys.modules["scigantic_pubchem.resolve"]

_PER_CHUNK_LATENCY = 0.2


def _slow_fetch_chunk(chunk):
    time.sleep(_PER_CHUNK_LATENCY)
    return [
        _resolve_module.Compound(
            cid=int(c),
            title=None,
            smiles=None,
            inchi=None,
            inchi_key=None,
            iupac_name=None,
            molecular_formula=None,
            molecular_weight=None,
        )
        for c in chunk
    ]


def test_multiple_chunks_run_concurrently_not_sequentially():
    cids = list(range(1, 801))  # 800 CIDs -> 4 chunks of 200
    with mock.patch.object(_resolve_module, "_fetch_chunk", side_effect=_slow_fetch_chunk):
        start = time.monotonic()
        compounds = _resolve_module.resolve_many(cids)
        elapsed = time.monotonic() - start

    assert len(compounds) == 800
    # Sequential would be ~4 * 0.2s = 0.8s. 4 chunks with up to 8 free
    # workers should land close to one chunk's latency. Generous margin
    # against CI scheduling noise, while still well under sequential.
    assert elapsed < 3 * _PER_CHUNK_LATENCY


def test_single_chunk_skips_the_pool_entirely():
    cids = list(range(1, 51))  # 50 CIDs -> 1 chunk, no pool needed
    with mock.patch.object(_resolve_module, "_fetch_chunk", side_effect=_slow_fetch_chunk) as fetch:
        compounds = _resolve_module.resolve_many(cids)
    fetch.assert_called_once()
    assert len(compounds) == 50


def test_empty_input_returns_empty_without_calling_fetch():
    with mock.patch.object(_resolve_module, "_fetch_chunk", side_effect=_slow_fetch_chunk) as fetch:
        compounds = _resolve_module.resolve_many([])
    fetch.assert_not_called()
    assert compounds == []
