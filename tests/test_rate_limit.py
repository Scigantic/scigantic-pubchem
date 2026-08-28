"""_RateLimiter's token bucket: pure logic, no real network and no real
waiting. time.monotonic and time.sleep are both mocked with a shared fake
clock, since the bucket's correctness depends on real elapsed time between
acquire() calls -- something test_async_search.py's mocked-sleep-only
approach didn't need, because that test never depends on how much time a
mocked sleep represents, only on how many times it's called.
"""

import pytest

from scigantic_pubchem._client import _RateLimiter


def _fake_clock(monkeypatch):
    now = [0.0]
    monkeypatch.setattr("time.monotonic", lambda: now[0])
    monkeypatch.setattr("time.sleep", lambda seconds: now.__setitem__(0, now[0] + seconds))
    return now


def test_capacity_tokens_available_with_no_wait(monkeypatch):
    _fake_clock(monkeypatch)
    limiter = _RateLimiter(rate=5.0, capacity=5.0)
    for _ in range(5):
        limiter.acquire()  # every one of the starting 5 tokens, no sleep needed
    assert limiter._tokens == pytest.approx(0.0, abs=1e-9)


def test_exhausted_bucket_waits_exactly_one_token_worth(monkeypatch):
    now = _fake_clock(monkeypatch)
    limiter = _RateLimiter(rate=5.0, capacity=1.0)
    limiter.acquire()  # drains the single starting token immediately
    limiter.acquire()  # bucket empty: must wait 1/rate seconds for the next one
    assert now[0] == pytest.approx(1 / 5.0)


def test_refills_at_the_configured_rate_not_faster(monkeypatch):
    now = _fake_clock(monkeypatch)
    limiter = _RateLimiter(rate=2.0, capacity=1.0)
    limiter.acquire()
    now[0] += 0.1  # half a token's worth of time has passed, not a full one
    limiter.acquire()
    # Should have waited the remaining 0.4s (0.5s needed - 0.1s already
    # elapsed), landing at 0.5s total, not 0.1s.
    assert now[0] == pytest.approx(0.5)


def test_acquire_is_safe_from_multiple_threads(monkeypatch):
    # Frozen clock: no real time passes, so no refill happens either, and
    # capacity exactly matches the number of calls. The only way this lands
    # on exactly zero after 200 concurrent acquires from 8 real threads is
    # if the lock correctly serializes every read-check-decrement with no
    # lost updates; a real (unfrozen) clock can't give an exact assertion
    # here, since wall-clock jitter across threads refills a real, variable
    # number of tokens between calls.
    from concurrent.futures import ThreadPoolExecutor

    monkeypatch.setattr("time.monotonic", lambda: 0.0)
    limiter = _RateLimiter(rate=1.0, capacity=200.0)
    with ThreadPoolExecutor(max_workers=8) as pool:
        list(pool.map(lambda _: limiter.acquire(), range(200)))
    assert limiter._tokens == pytest.approx(0.0, abs=1e-9)
