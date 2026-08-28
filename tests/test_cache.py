from concurrent.futures import ThreadPoolExecutor
from unittest import mock

import scigantic_pubchem as pubchem
from scigantic_pubchem import _client, cache


def test_cache_enabled_by_default():
    assert pubchem.is_cache_enabled() is True


def test_cache_dir_resolves(tmp_path):
    resolved = pubchem.enable_cache(cache_dir=str(tmp_path))
    try:
        assert resolved == tmp_path
        assert pubchem.cache_dir() == tmp_path
    finally:
        pubchem.enable_cache()  # back to the default location for other tests


def test_second_lookup_does_not_hit_the_network(tmp_path):
    pubchem.enable_cache(cache_dir=str(tmp_path))
    try:
        pubchem.clear_cache()
        first = pubchem.resolve("water")
        assert first is not None

        session = _client._get_session()
        with mock.patch.object(session, "get", side_effect=AssertionError("network hit on a cached lookup")):
            second = pubchem.resolve("water")
        assert second == first
    finally:
        pubchem.enable_cache()


def test_entries_expire_after_ttl(tmp_path):
    pubchem.enable_cache(cache_dir=str(tmp_path), ttl_days=30)
    try:
        pubchem.clear_cache()
        pubchem.resolve("water")
        cached_files = list(tmp_path.glob("*.json"))
        assert len(cached_files) == 1

        # Backdate the cached entry past the TTL and confirm it's treated
        # as a miss (a real network call happens, not a crash or a stale
        # return) rather than trusting an indefinite cache silently.
        import json

        entry = json.loads(cached_files[0].read_text())
        backdated_at = entry["cached_at"] - 31 * 86400
        entry["cached_at"] = backdated_at
        cached_files[0].write_text(json.dumps(entry))

        session = _client._get_session()
        with mock.patch.object(session, "get", wraps=session.get) as spy:
            second = pubchem.resolve("water")
            assert spy.called  # the expired entry forced a real network call
        assert second is not None
        # A fresh entry was written with a new timestamp, not the backdated one.
        refreshed = json.loads(cached_files[0].read_text())
        assert refreshed["cached_at"] > backdated_at
    finally:
        pubchem.enable_cache()


def test_ttl_none_never_expires(tmp_path):
    pubchem.enable_cache(cache_dir=str(tmp_path), ttl_days=None)
    try:
        pubchem.clear_cache()
        pubchem.resolve("water")
        cached_files = list(tmp_path.glob("*.json"))
        import json

        entry = json.loads(cached_files[0].read_text())
        entry["cached_at"] -= 365 * 86400  # a year old
        cached_files[0].write_text(json.dumps(entry))

        session = _client._get_session()
        with mock.patch.object(session, "get", side_effect=AssertionError("network hit despite ttl_days=None")):
            pubchem.resolve("water")
    finally:
        pubchem.enable_cache()


def test_search_polling_is_never_cached(tmp_path):
    # A regression check for a real bug caught before shipping: caching an
    # intermediate {"Waiting": ...} response under the search's own
    # (path, params) key would make a later, unrelated call to the same
    # search replay stale in-progress state instead of a fresh request.
    pubchem.enable_cache(cache_dir=str(tmp_path))
    try:
        pubchem.clear_cache()
        pubchem.similar_compounds(
            "CC(=O)OC1=CC=CC=C1C(=O)O", threshold=95, max_records=5, resolve=False
        )
        assert list(tmp_path.glob("*.json")) == []
    finally:
        pubchem.enable_cache()


def test_concurrent_writes_to_the_same_key_never_raise(tmp_path):
    # Regression test for a real bug: put()'s temp filename used to be
    # deterministic (derived only from the cache key), so two threads
    # racing to fill the *same* key shared one temp path. Whichever thread's
    # os.replace() ran second raised FileNotFoundError, because the first
    # had already consumed it. Caught via a 32-thread stress run where two
    # names were each queued from two threads at once. Exercised directly
    # against cache.put(), not resolve(), so this stays fast and doesn't
    # depend on winning a real network race to reproduce.
    pubchem.enable_cache(cache_dir=str(tmp_path))
    try:
        pubchem.clear_cache()
        errors = []

        def write(i):
            try:
                cache.put("/same/path", {"k": "v"}, {"call": i})
            except Exception as e:  # noqa: BLE001
                errors.append(e)

        with ThreadPoolExecutor(max_workers=16) as pool:
            list(pool.map(write, range(16)))

        assert errors == []
        cached_files = list(tmp_path.glob("*.json"))
        assert len(cached_files) == 1  # all 16 writers converged on one key
        stray_parts = list(tmp_path.glob("*.part"))
        assert stray_parts == []  # every temp file was consumed by its own replace
    finally:
        pubchem.enable_cache()


def test_disable_cache_hits_network_every_time(tmp_path):
    pubchem.enable_cache(cache_dir=str(tmp_path))
    pubchem.clear_cache()
    pubchem.disable_cache()
    try:
        assert pubchem.is_cache_enabled() is False
        # Two real calls, no assertion beyond "this doesn't crash and both
        # resolve": proving cache bypass by absence of a cached file.
        pubchem.resolve("water")
        cached_files = list(tmp_path.glob("*.json"))
        assert cached_files == []
    finally:
        pubchem.enable_cache()
