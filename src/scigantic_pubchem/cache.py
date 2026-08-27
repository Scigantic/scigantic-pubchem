"""Local response cache -- ON by default, unlike scigantic-chembl's and
scigantic-bindingdb's caching (both default OFF).

That's a deliberate difference, not an inconsistency: those two packages
read from a public S3 mirror with no meaningful rate limit, so caching is a
pure convenience, opt-in because zero-setup-by-default is the point. This
package calls a rate-limited live API (PubChem's PUG REST) for every
lookup, so re-fetching the same identifier repeatedly in a notebook loop is
both slow and the exact thing that trips the throttling this package
otherwise works to avoid (see _client.py). Caching stays on unless turned
off explicitly.

Every PUG REST response is cached keyed by its exact request path + params,
so cache correctness doesn't depend on any function here knowing what a
"compound" or "cid" is -- see _client.py's request().

Cached entries expire after ttl_days (30 by default). This package's whole
argument for not mirroring PubChem is that a live query is more correct
than a stale snapshot -- an indefinitely-cached response would quietly
recreate that exact staleness inside this package instead. 30 days is long
enough to make a notebook session or a multi-day analysis fast without
repeatedly hitting PUG REST, short enough that the cache doesn't silently
diverge from PubChem for months. Pass ttl_days=None to disable expiry
entirely if that tradeoff is wrong for a specific use case.

Reads and writes of individual entries are safe to call concurrently (the
write path writes to a temp file and os.replace()s it into place, so a
reader never observes a partial write). enable_cache()/disable_cache()
themselves are not synchronized against concurrent reads -- like most
one-time configuration calls (comparable to mutating os.environ), they are
meant to be called once at the start of a script or session, not toggled
from multiple threads at once.
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from pathlib import Path
from typing import Any

_enabled = True
_cache_dir: Path | None = None
_ttl_seconds: float | None = 30 * 86400


def _default_cache_dir() -> Path:
    if sys.platform == "win32":
        base = os.environ.get("LOCALAPPDATA") or str(Path.home() / "AppData" / "Local")
    elif sys.platform == "darwin":
        base = str(Path.home() / "Library" / "Caches")
    else:
        base = os.environ.get("XDG_CACHE_HOME") or str(Path.home() / ".cache")
    return Path(base) / "scigantic-pubchem"


def _resolve_dir() -> Path:
    global _cache_dir
    if _cache_dir is None:
        env = os.environ.get("SCIGANTIC_PUBCHEM_CACHE")
        _cache_dir = Path(env) if env else _default_cache_dir()
        _cache_dir.mkdir(parents=True, exist_ok=True)
    return _cache_dir


def enable_cache(cache_dir: str | None = None, ttl_days: float | None = 30) -> Path:
    """Turn caching on (it already is, by default) and optionally point it
    at a specific directory and/or change how long an entry stays valid.

    ttl_days=None disables expiry -- entries are reused forever until
    clear_cache() or disable_cache(). Returns the resolved directory.
    """
    global _enabled, _cache_dir, _ttl_seconds
    if cache_dir is not None:
        _cache_dir = Path(cache_dir)
        _cache_dir.mkdir(parents=True, exist_ok=True)
    else:
        _resolve_dir()
    _ttl_seconds = ttl_days * 86400 if ttl_days is not None else None
    _enabled = True
    return _cache_dir  # type: ignore[return-value]


def disable_cache() -> None:
    """Turn caching off. Every call hits PUG REST fresh until re-enabled."""
    global _enabled
    _enabled = False


def is_cache_enabled() -> bool:
    return _enabled


def cache_dir() -> Path:
    return _resolve_dir()


def _key(path: str, params: dict[str, Any] | None) -> str:
    raw = json.dumps({"path": path, "params": params or {}}, sort_keys=True)
    return hashlib.sha256(raw.encode()).hexdigest()


def get(path: str, params: dict[str, Any] | None) -> dict[str, Any] | None:
    if not _enabled:
        return None
    file = _resolve_dir() / f"{_key(path, params)}.json"
    if not file.exists():
        return None
    try:
        entry = json.loads(file.read_text())
    except (json.JSONDecodeError, OSError):
        return None
    if _ttl_seconds is not None and time.time() - entry.get("cached_at", 0) > _ttl_seconds:
        file.unlink(missing_ok=True)
        return None
    return entry["value"]  # type: ignore[no-any-return]


def put(path: str, params: dict[str, Any] | None, value: dict[str, Any]) -> None:
    if not _enabled:
        return
    file = _resolve_dir() / f"{_key(path, params)}.json"
    tmp = file.with_suffix(".json.part")
    tmp.write_text(json.dumps({"cached_at": time.time(), "value": value}))
    os.replace(tmp, file)


def clear() -> int:
    """Delete every cached response. Returns how many files were removed."""
    d = _resolve_dir()
    n = 0
    for f in d.glob("*.json"):
        f.unlink()
        n += 1
    return n
