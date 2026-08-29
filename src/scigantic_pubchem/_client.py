"""HTTP client for PubChem's PUG REST API, with the three things PubChemPy
does not do: respect PubChem's own live rate-limit signal, pace requests
against PubChem's documented limit before they're ever sent, and retry
transient failures instead of raising immediately.

Verified 2026-08-27: every PUG REST response carries an X-Throttling-Control
header reporting live status across three dimensions, e.g.
    Request Count status: Green (0%), Request Time status: Green (0%), Service status: Green (27%)
This is PubChem's own documented, real-time signal for how close a caller is
to being throttled. PubChemPy 1.0.5's source (inspected directly, not
assumed) never reads this header at all: on an HTTPError it raises
immediately with no retry, and it has no rate-awareness beyond that. Every
call here reads the header and backs off proactively when status is not
"Green" on any dimension, rather than waiting to be rejected.

That header-based backoff is still reactive: it only slows a caller down
after PubChem has already reported elevated load. A burst of concurrent
requests (resolve_many()'s parallel chunks, or a caller's own thread pool)
can fire well past PubChem's documented ceiling of 5 requests/second before
any of them has seen a throttle signal at all. _RateLimiter below is a
token bucket that every request acquires from first, so a burst is paced
to the documented limit up front instead of recovering from it after the
fact.
"""

from __future__ import annotations

import re
import threading
import time
import warnings
from typing import Any

import requests

from . import cache

_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
_USER_AGENT = "scigantic-pubchem/0.1.0 (+https://scigantic.com; mailto:support@scigantic.com)"

_THROTTLE_RE = re.compile(r"(\w[\w ]*?) status: (Green|Yellow|Red) \((\d+)%\)")

# Backoff seconds by worst observed status across the three dimensions.
_BACKOFF_SECONDS = {"Green": 0.0, "Yellow": 2.0, "Red": 10.0}

_MAX_RETRIES = 5
_RETRY_STATUS_CODES = {429, 500, 502, 503, 504}


class PubChemError(Exception):
    """Raised for a PUG REST error response (4xx/5xx after retries exhausted)."""


class CompoundNotFoundError(PubChemError):
    """Raised for a 404 (PUGREST.NotFound): a real, expected outcome for
    an unmatched identifier, not a transient failure. Verified live
    2026-08-27: PUG REST returns HTTP 404 with a JSON Fault body for a
    miss, not an empty 200. resolve() catches this and returns None
    rather than letting every caller handle the exception itself."""


def _parse_throttle_header(value: str | None) -> str:
    """Worst status ('Green'/'Yellow'/'Red') across all reported dimensions."""
    if not value:
        return "Green"
    statuses = [m.group(2) for m in _THROTTLE_RE.finditer(value)]
    if "Red" in statuses:
        return "Red"
    if "Yellow" in statuses:
        return "Yellow"
    return "Green"


class _RateLimiter:
    """A token bucket, so a burst of requests is paced up front instead of
    relying on X-Throttling-Control to catch it after the fact.

    capacity tokens available immediately (a small burst is fine, PubChem's
    limit is itself a rate, not a hard cap on simultaneous requests), then
    refilling at rate tokens/second. acquire() blocks the calling thread
    until a token is available; safe to call from multiple threads at once,
    the way resolve_many()'s parallel chunks and a caller's own
    ThreadPoolExecutor both do.
    """

    def __init__(self, rate: float, capacity: float) -> None:
        self._rate = rate
        self._capacity = capacity
        self._tokens = capacity
        self._last = time.monotonic()
        self._lock = threading.Lock()

    def acquire(self) -> None:
        while True:
            with self._lock:
                now = time.monotonic()
                self._tokens = min(self._capacity, self._tokens + (now - self._last) * self._rate)
                self._last = now
                if self._tokens >= 1:
                    self._tokens -= 1
                    return
                wait = (1 - self._tokens) / self._rate
            time.sleep(wait)


# PubChem's documented limit: no more than 5 requests/second per user.
_limiter = _RateLimiter(rate=5.0, capacity=5.0)

_session: requests.Session | None = None
_session_lock = threading.Lock()


def _get_session() -> requests.Session:
    """Lazily create the shared Session under a lock, so two threads
    resolving identifiers concurrently (a real, plausible pattern: a
    ThreadPoolExecutor over a list of names) can't both see _session as
    None and each construct one, leaking a connection pool. requests.Session
    itself is documented thread-safe for issuing requests once constructed;
    only the lazy-init race needed guarding."""
    global _session
    if _session is None:
        with _session_lock:
            if _session is None:  # re-check: another thread may have won the race
                _session = requests.Session()
                _session.headers["User-Agent"] = _USER_AGENT
    return _session


def request(
    path: str,
    params: dict[str, Any] | None = None,
    method: str = "GET",
    cacheable: bool = True,
) -> dict[str, Any]:
    """GET (or POST, for long identifier lists) a PUG REST path and return
    the parsed JSON body.

    path is appended to https://pubchem.ncbi.nlm.nih.gov/rest/pug, e.g.
    "/compound/name/aspirin/property/MolecularWeight/JSON".

    Paces every attempt through the shared token bucket (PubChem's
    documented 5 req/s) before it's sent, backs off further based on
    X-Throttling-Control before it becomes a 429, and retries 429/5xx with
    exponential backoff (PubChemPy does none of this) rather than raising
    on the first transient failure.

    Cached (see cache.py) keyed by (path, params) when caching is enabled
    (the default) and cacheable=True (the default). cacheable=False is for
    request_search()'s polling calls specifically: caching an intermediate
    "still waiting" response under the search's own (path, params) key
    would mean a later, unrelated call to the same search replays that
    stale in-progress state instead of making a fresh request.
    """
    if method == "GET" and cacheable:
        cached = cache.get(path, params)
        if cached is not None:
            return cached

    session = _get_session()
    url = f"{_BASE}{path}"
    last_exc: Exception | None = None

    for attempt in range(_MAX_RETRIES):
        _limiter.acquire()
        try:
            if method == "POST":
                response = session.post(url, data=params, timeout=30)
            else:
                response = session.get(url, params=params, timeout=30)
        except requests.RequestException as exc:
            last_exc = exc
            time.sleep(2 ** attempt)
            continue

        status = _parse_throttle_header(response.headers.get("X-Throttling-Control"))
        if response.status_code in _RETRY_STATUS_CODES and attempt < _MAX_RETRIES - 1:
            wait = max(_BACKOFF_SECONDS[status], 2 ** attempt)
            warnings.warn(
                f"PUG REST returned {response.status_code} (throttle status {status}); "
                f"retrying in {wait:.1f}s (attempt {attempt + 1}/{_MAX_RETRIES})",
                stacklevel=2,
            )
            time.sleep(wait)
            continue
        if response.status_code == 404:
            raise CompoundNotFoundError(f"no match for {url}")
        if response.status_code >= 400:
            raise PubChemError(f"PUG REST {response.status_code} for {url}: {response.text[:500]}")

        if status != "Green":
            # Proactive backoff even on success, so a run of calls doesn't
            # walk itself into a 429 a few requests later.
            time.sleep(_BACKOFF_SECONDS[status])
        body: dict[str, Any] = response.json()
        if method == "GET" and cacheable:
            cache.put(path, params, body)
        return body

    raise PubChemError(f"PUG REST request failed after {_MAX_RETRIES} attempts: {last_exc}")


_POLL_INITIAL_SECONDS = 0.5
_POLL_MAX_SECONDS = 5.0
_POLL_BACKOFF_FACTOR = 2.0


def request_search(path: str, params: dict[str, Any] | None = None) -> dict[str, Any]:
    """Like request(), for PUG REST's search-type endpoints (fastsimilarity_2d,
    fastsubstructure, fastformula, ...), which can respond either
    immediately or asynchronously.

    A slow search returns {"Waiting": {"ListKey": "..."}} instead of a
    result, and this polls /compound/listkey/{key}/cids/JSON until the real
    result is ready. The wait between polls starts at
    _POLL_INITIAL_SECONDS and doubles each time up to _POLL_MAX_SECONDS,
    rather than a flat interval: a fixed 2s wait (what PubChemPy uses, its
    source was read directly to confirm the underlying protocol, not
    assumed from docs) checks sooner than necessary for a job that's about
    to finish and, symmetrically, checks far more often than necessary for
    one that's going to take 30-60s. A real, measured case burns 15-30
    wasted round trips at a flat interval. Starting below 2s actually
    catches a fast job sooner than the old fixed wait did; the doubling is
    what keeps a slow job from paying for that on every subsequent poll.

    Neither the initial call nor the polling calls are cached
    (cacheable=False): a ListKey is single-use, and caching an
    intermediate "still waiting" response under the search's own
    (path, params) key would mean a later, unrelated call to the same
    search replays that stale in-progress state instead of making a fresh
    request. Callers that want to avoid repeat searches should cache the
    resolved CIDs themselves.
    """
    body = request(path, params, method="GET", cacheable=False)
    poll_wait = _POLL_INITIAL_SECONDS
    while "Waiting" in body and "ListKey" in body.get("Waiting", {}):
        time.sleep(poll_wait)
        poll_wait = min(poll_wait * _POLL_BACKOFF_FACTOR, _POLL_MAX_SECONDS)
        listkey = body["Waiting"]["ListKey"]
        body = request(f"/compound/listkey/{listkey}/cids/JSON", method="GET", cacheable=False)
    return body
