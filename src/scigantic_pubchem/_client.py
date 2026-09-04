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

import json
import os
import re
import threading
import time
import uuid
import warnings
from pathlib import Path
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


def _send_with_retry(
    method: str,
    url: str,
    params: dict[str, Any] | None,
    stream: bool = False,
    timeout: float = 30.0,
) -> requests.Response:
    """The retry/backoff loop request() and stream_to_file() both need:
    paces every attempt through the token bucket, retries 429/5xx with
    backoff informed by X-Throttling-Control, and raises
    CompoundNotFoundError/PubChemError the same way for either caller.

    Returns the still-open Response on success. stream=True (used by
    stream_to_file(), which writes a CSV that can run to tens of MB for a
    large assay directly to disk) defers reading the body so the caller
    can iter_content() it rather than this function buffering the whole
    thing first just to hand it back.

    timeout defaults to 30s, the ceiling every existing caller already
    ran under. bioassay.py's dose_response()/download_dose_response()
    pass a much larger value: measured live 2026-08-31 against AID 1851,
    PUG REST's full (non-concise) assay Data Table operation does not
    scale linearly -- 250 SIDs returned in under a second, 2000 SIDs took
    94s against the same server -- so a 30s timeout would fail requests
    that are genuinely still in progress, not stuck.
    """
    session = _get_session()
    last_exc: Exception | None = None

    for attempt in range(_MAX_RETRIES):
        _limiter.acquire()
        try:
            if method == "POST":
                response = session.post(url, data=params, timeout=timeout, stream=stream)
            else:
                response = session.get(url, params=params, timeout=timeout, stream=stream)
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
        return response

    raise PubChemError(f"PUG REST request failed after {_MAX_RETRIES} attempts: {last_exc}")


def _parse_json(response: requests.Response, url: str) -> dict[str, Any]:
    """response.json(), with one fallback for a real, observed PubChem
    quirk: a large `concise`/`assaysummary` response can carry a
    depositor's raw, unescaped ASCII control character inside a free-text
    field (an assay comment or description), which the standard decoder
    (strict=True, the default for both json.loads and Response.json())
    rejects outright with "Invalid control character" -- verified live
    2026-09-04 against a ~20MB response mid-batch. Not a transient/network
    problem (re-requesting the same identifier returns the identical
    bytes, verified by reparsing the same response text rather than
    refetching it), so this doesn't retry over the network for it.

    json.loads's own strict=False is the documented escape hatch for
    exactly this: it allows control characters (0x00-0x1F, including a
    literal tab/newline/null) inside a string instead of requiring them
    escaped as \\u00XX. Tried only after the strict parse fails, so a
    normal response pays no extra cost.
    """
    try:
        body: dict[str, Any] = response.json()
    except ValueError:
        try:
            body = json.loads(response.text, strict=False)
        except ValueError as exc:
            raise PubChemError(
                f"PUG REST response for {url} is not valid JSON, even with relaxed "
                f"control-character handling: {exc}"
            ) from exc
    return body


def request(
    path: str,
    params: dict[str, Any] | None = None,
    method: str = "GET",
    cacheable: bool = True,
    timeout: float = 30.0,
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

    url = f"{_BASE}{path}"
    response = _send_with_retry(method, url, params, timeout=timeout)
    body = _parse_json(response, url)
    if method == "GET" and cacheable:
        cache.put(path, params, body)
    return body


def request_text(
    path: str,
    params: dict[str, Any] | None = None,
    method: str = "GET",
    timeout: float = 30.0,
) -> str:
    """Like request(), for a PUG REST operation whose response is CSV, not
    JSON -- bioassay.py's dose_response()/download_dose_response(), which
    read the plain (non-concise) `/assay/aid/{aid}/CSV` operation. That
    operation returns a flat table directly, unlike the same AID's plain
    JSON/XML (PC_AssaySubmit, the deeply nested depositor record this
    package otherwise avoids -- verified live 2026-08-31), so CSV is the
    only usable format here.

    Not cached, for the same reason stream_to_file() isn't: a chunk of
    this data can run to single-digit MB on its own, and
    download_dose_response()'s on-disk progress tracking is already the
    thing that avoids repeat work across a resumed pull.
    """
    url = f"{_BASE}{path}"
    response = _send_with_retry(method, url, params, timeout=timeout)
    return response.text


def stream_to_file(path: str, dest: Path | str, params: dict[str, Any] | None = None) -> Path:
    """GET a PUG REST path and write the response body straight to dest,
    chunk by chunk, rather than buffering it in memory as request() does.

    For bioassay.download_assay_results(): a large screen's `concise` CSV
    can run to tens of MB (measured: a 69,000-row qHTS assay's concise CSV
    is ~9.5MB), and every other function in this package holds its PUG
    REST response in memory only as long as it takes to build the typed
    records this package returns. A download is meant to end up as a file
    on disk regardless, so there's nothing to gain by holding the whole
    body as one Python bytes object in between.

    Not cached (see request()'s cacheable=False for the same reasoning
    applied to search polling): the destination file already is the cached
    artifact.
    """
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    url = f"{_BASE}{path}"
    response = _send_with_retry("GET", url, params, stream=True)
    tmp = dest.with_suffix(f"{dest.suffix}.{uuid.uuid4().hex}.part")
    try:
        with tmp.open("wb") as f:
            for chunk in response.iter_content(chunk_size=1 << 16):
                if chunk:
                    f.write(chunk)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise
    else:
        os.replace(tmp, dest)
    return dest


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

    The polling calls themselves are never cached (request(..., cacheable=False)
    for each one): a ListKey is single-use, so there is nothing correct to
    key a cache entry on besides the search's own (path, params), and an
    intermediate "still waiting" response cached under that key would mean
    a later, unrelated call to the same search replays stale in-progress
    state instead of making a fresh request.

    The search itself is checked against the cache first, and the final
    resolved body (never an intermediate "Waiting" one) is written back
    under (path, params) once polling ends, the same key request() would
    use for any other idempotent GET. Before this, a search's own result
    was never cached at all, so re-running the same similarity/substructure
    query, e.g. resuming a batch job after a crash partway through, always
    repeated the full live search, including any async poll, even though
    resolve_many()'s follow-up lookup for the returned CIDs was already
    cached. Same 30-day TTL and staleness tradeoff as everything else (see
    cache.py); a caller that wants a guaranteed-fresh corpus search should
    disable_cache() or clear_cache() rather than expect this endpoint to be
    special-cased.
    """
    cached = cache.get(path, params)
    if cached is not None:
        return cached
    body = request(path, params, method="GET", cacheable=False)
    poll_wait = _POLL_INITIAL_SECONDS
    while "Waiting" in body and "ListKey" in body.get("Waiting", {}):
        time.sleep(poll_wait)
        poll_wait = min(poll_wait * _POLL_BACKOFF_FACTOR, _POLL_MAX_SECONDS)
        listkey = body["Waiting"]["ListKey"]
        body = request(f"/compound/listkey/{listkey}/cids/JSON", method="GET", cacheable=False)
    cache.put(path, params, body)
    return body
