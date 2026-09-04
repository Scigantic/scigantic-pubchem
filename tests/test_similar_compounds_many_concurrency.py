"""similar_compounds_many()'s parallel dispatch: mocked, so the parallelism
claim is verified deterministically by wall-clock elapsed time, the same
approach test_resolve_many_concurrency.py uses for resolve_many() and for
the same reason -- a live call's real latency can't distinguish sequential
from concurrent execution either way.
"""

import time
from unittest import mock

import scigantic_pubchem as pubchem
from scigantic_pubchem import similarity

_PER_QUERY_LATENCY = 0.2


def _slow_similar_compounds(smiles, threshold=90, max_records=100, resolve=True):
    time.sleep(_PER_QUERY_LATENCY)
    return [len(smiles)]  # cheap, deterministic stand-in for a real CID list


def test_multiple_queries_run_concurrently_not_sequentially():
    queries = [f"C{i}" for i in range(10)]  # 10 distinct one-carbon-ish queries
    with mock.patch.object(similarity, "similar_compounds", side_effect=_slow_similar_compounds):
        start = time.monotonic()
        results = pubchem.similar_compounds_many(queries, resolve=False)
        elapsed = time.monotonic() - start

    assert set(results) == set(queries)
    # Sequential would be ~10 * 0.2s = 2.0s. 10 queries over the default
    # 5-worker pool should land close to two queries' latency in sequence
    # per worker (~0.4s). Generous margin against CI scheduling noise,
    # while still well under sequential.
    assert elapsed < 5 * _PER_QUERY_LATENCY


def test_each_query_is_dispatched_individually():
    queries = ["C", "CC", "CCC"]
    with mock.patch.object(similarity, "similar_compounds", side_effect=_slow_similar_compounds) as fetch:
        results = pubchem.similar_compounds_many(queries, resolve=False, max_workers=5)
    assert fetch.call_count == 3
    assert results == {"C": [1], "CC": [2], "CCC": [3]}


def test_repeated_query_collapses_to_one_key():
    with mock.patch.object(similarity, "similar_compounds", side_effect=_slow_similar_compounds) as fetch:
        results = pubchem.similar_compounds_many(["C", "C", "CC"], resolve=False)
    assert fetch.call_count == 3  # still called once per input, not deduped up front
    assert results == {"C": [1], "CC": [2]}
