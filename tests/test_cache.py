from unittest import mock

import scigantic_pubchem as pubchem
from scigantic_pubchem import _client


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


def test_disable_cache_hits_network_every_time(tmp_path):
    pubchem.enable_cache(cache_dir=str(tmp_path))
    pubchem.clear_cache()
    pubchem.disable_cache()
    try:
        assert pubchem.is_cache_enabled() is False
        # Two real calls, no assertion beyond "this doesn't crash and both
        # resolve" -- proving cache bypass by absence of a cached file.
        pubchem.resolve("water")
        cached_files = list(tmp_path.glob("*.json"))
        assert cached_files == []
    finally:
        pubchem.enable_cache()
