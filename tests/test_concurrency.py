"""Concurrent use from multiple threads: a real, plausible pattern (e.g.
a ThreadPoolExecutor resolving a list of names) that the correctness tests
elsewhere don't exercise. Confirms the lazily-created shared Session (see
_client._get_session's lock) survives concurrent first-use without
crashing or duplicating, and that a burst of concurrent lookups all
complete correctly.
"""

from concurrent.futures import ThreadPoolExecutor

import scigantic_pubchem as pubchem
from scigantic_pubchem import _client


def test_concurrent_first_use_creates_one_session():
    _client._session = None  # force the lazy-init race on every thread
    names = ["aspirin", "caffeine", "acetaminophen", "water", "ethanol", "glucose", "ibuprofen", "morphine"]

    with ThreadPoolExecutor(max_workers=8) as pool:
        results = list(pool.map(pubchem.resolve, names))

    assert all(r is not None for r in results)
    assert len({r.cid for r in results}) == len(names)  # all distinct, all resolved
    assert _client._session is not None  # exactly one session object exists now
