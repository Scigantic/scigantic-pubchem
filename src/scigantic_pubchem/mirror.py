"""Read s3://scigantic-pubchem, a weekly-refreshed local parquet mirror of
PubChem's compound identifier/name/mass registry (smiles, titles,
iupac_names, inchi_keys, mass, parent, cid_sid, synonyms -- see that
bucket's _MIRROR.txt for the full table list and scope).

CID-KEYED ONLY. Every table in this mirror is exported in ascending CID
order, so a query filtered on cid gets real row-group pruning from a remote
parquet read (fast: only the matching row group is fetched over the
network). A query filtered on anything else -- InChIKey, name, SMILES --
gets none of that benefit and means scanning the whole column remotely, the
exact problem that ruled out a full-registry mirror in the first place (see
the archive's schema card). Use resolve()/xrefs() from this package's live
PUG REST wrapper for name/InChIKey/SMILES -> CID resolution; use this module
once you already have a CID (or a DataFrame keyed on one).

For real bulk/offline work -- loading this into your own DuckDB session
without paying network round trips per query -- use download() to pull the
parquet files down once, then query them locally with duckdb/pandas
directly. That is the direct fix for "I had to write my own script to
reformat a bulk file and load it into DuckDB myself": here, the file is
already a ready-made local parquet, no reformatting step.

Needs duckdb for the live httpfs functions: `pip install "scigantic-pubchem[bridge]"`.
download() only needs requests, already a base dependency.
"""

from __future__ import annotations

import os
from typing import TYPE_CHECKING, Any

from .bridge import _get_connection

if TYPE_CHECKING:
    import pandas as pd

_MIRROR_BUCKET = "scigantic-pubchem"

_TABLES = (
    "smiles",
    "titles",
    "iupac_names",
    "inchi_keys",
    "mass",
    "parent",
    "cid_sid",
    "synonyms",
)


def _table_path(table: str) -> str:
    if table not in _TABLES:
        raise ValueError(f"unknown mirror table {table!r}; known tables: {', '.join(_TABLES)}")
    return f"s3://{_MIRROR_BUCKET}/parquet/{table}.parquet"


def identifiers(cid: int) -> dict[str, Any] | None:
    """SMILES, title, IUPAC name, InChI/InChIKey (first listed), and mass
    for a single CID, joined live across the mirror's per-table files.

    Returns None if the CID is not in the mirror. A CID can have more than
    one InChIKey row (one per protonation state, see the archive's schema
    card); this returns the first one DuckDB happens to read, which is
    usually but not guaranteed to be the parent/neutral form -- use
    inchi_keys(cid) directly if you need all of them.
    """
    con = _get_connection()
    row = con.execute(
        f"""
        SELECT s.smiles, t.title, i.iupac_name, k.inchi, k.inchikey,
               m.molecular_formula, m.monoisotopic_mass, m.exact_mass
        FROM read_parquet('{_table_path("smiles")}') s
        LEFT JOIN read_parquet('{_table_path("titles")}') t USING (cid)
        LEFT JOIN read_parquet('{_table_path("iupac_names")}') i USING (cid)
        LEFT JOIN read_parquet('{_table_path("inchi_keys")}') k USING (cid)
        LEFT JOIN read_parquet('{_table_path("mass")}') m USING (cid)
        WHERE s.cid = ?
        LIMIT 1
        """,
        [cid],
    ).fetchone()
    if row is None:
        return None
    return {
        "cid": cid,
        "smiles": row[0],
        "title": row[1],
        "iupac_name": row[2],
        "inchi": row[3],
        "inchikey": row[4],
        "molecular_formula": row[5],
        "monoisotopic_mass": row[6],
        "exact_mass": row[7],
    }


def inchi_keys(cid: int) -> "pd.DataFrame":
    """All (InChI, InChIKey) rows for a CID -- one per protonation state."""
    con = _get_connection()
    return con.execute(
        f"SELECT inchi, inchikey FROM read_parquet('{_table_path('inchi_keys')}') WHERE cid = ?",
        [cid],
    ).df()


def synonyms(cid: int) -> list[str]:
    """All known names for a CID, from PubChem's filtered synonym list."""
    con = _get_connection()
    rows = con.execute(
        f"SELECT synonym FROM read_parquet('{_table_path('synonyms')}') WHERE cid = ?",
        [cid],
    ).fetchall()
    return [r[0] for r in rows]


def parent(cid: int) -> int | None:
    """This CID's parent CID (a CID may be its own parent, or have none)."""
    con = _get_connection()
    row = con.execute(
        f"SELECT parent_cid FROM read_parquet('{_table_path('parent')}') WHERE cid = ?",
        [cid],
    ).fetchone()
    return row[0] if row else None


def substance_ids(cid: int) -> "pd.DataFrame":
    """SIDs this CID derives from, with link_type (1=standardized form of
    the deposited substance, 2=component of the standardized form), per
    PubChem's own README-Extras."""
    con = _get_connection()
    return con.execute(
        f"SELECT sid, link_type FROM read_parquet('{_table_path('cid_sid')}') WHERE cid = ?",
        [cid],
    ).df()


def download(dest_dir: str, tables: "list[str] | None" = None) -> dict[str, str]:
    """Download this mirror's parquet files to dest_dir for real local/
    offline use -- open them directly with duckdb.read_parquet() or
    pandas.read_parquet() afterward, no network round trip per query.

    Defaults to all eight tables (~15.5 GB total); pass tables=[...] to
    fetch only the ones you need (e.g. tables=["smiles", "titles"] for a
    ~3.4 GB subset). Skips a file that already exists at the destination
    with a matching size, so re-running this after a partial download or a
    weekly mirror refresh only fetches what changed.

    Returns {table_name: local_path}.
    """
    import requests

    want = tables if tables is not None else list(_TABLES)
    unknown = [t for t in want if t not in _TABLES]
    if unknown:
        raise ValueError(f"unknown mirror table(s) {unknown}; known tables: {', '.join(_TABLES)}")

    os.makedirs(dest_dir, exist_ok=True)
    out: dict[str, str] = {}
    for table in want:
        url = f"https://{_MIRROR_BUCKET}.s3.amazonaws.com/parquet/{table}.parquet"
        local_path = os.path.join(dest_dir, f"{table}.parquet")
        head = requests.head(url, timeout=30)
        head.raise_for_status()
        remote_size = int(head.headers["Content-Length"])
        if os.path.exists(local_path) and os.path.getsize(local_path) == remote_size:
            out[table] = local_path
            continue
        with requests.get(url, stream=True, timeout=300) as resp:
            resp.raise_for_status()
            with open(local_path, "wb") as fh:
                for chunk in resp.iter_content(chunk_size=8 * 1024 * 1024):
                    fh.write(chunk)
        out[table] = local_path
    return out
