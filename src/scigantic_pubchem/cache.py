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
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
from pathlib import Path
from typing import Any

_enabled = True
_cache_dir: Path | None = None


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


def enable_cache(cache_dir: str | None = None) -> Path:
    """Turn caching on (it already is, by default) and optionally point it
    at a specific directory. Returns the resolved directory."""
    global _enabled, _cache_dir
    if cache_dir is not None:
        _cache_dir = Path(cache_dir)
        _cache_dir.mkdir(parents=True, exist_ok=True)
    else:
        _resolve_dir()
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
        return json.loads(file.read_text())  # type: ignore[no-any-return]
    except (json.JSONDecodeError, OSError):
        return None


def put(path: str, params: dict[str, Any] | None, value: dict[str, Any]) -> None:
    if not _enabled:
        return
    file = _resolve_dir() / f"{_key(path, params)}.json"
    tmp = file.with_suffix(".json.part")
    tmp.write_text(json.dumps(value))
    os.replace(tmp, file)


def clear() -> int:
    """Delete every cached response. Returns how many files were removed."""
    d = _resolve_dir()
    n = 0
    for f in d.glob("*.json"):
        f.unlink()
        n += 1
    return n
