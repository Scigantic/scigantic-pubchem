"""Cross-reference a PubChem compound into scigantic-chembl and
scigantic-bindingdb -- both live, on-demand DuckDB httpfs queries against
those public mirrors, not a precomputed table shipped with this package.

Needs duckdb: `pip install "scigantic-pubchem[bridge]"`. Everything else in
this package works without it.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from .xrefs import chembl_id as _chembl_id

if TYPE_CHECKING:
    import duckdb
    import pandas as pd

_CHEMBL_BUCKET = "scigantic-chembl"
_BINDINGDB_BUCKET = "scigantic-bindingdb"


def _connect() -> "duckdb.DuckDBPyConnection":
    try:
        import duckdb
    except ImportError as exc:
        raise ImportError(
            "this needs duckdb. Install with: pip install 'scigantic-pubchem[bridge]'"
        ) from exc
    con = duckdb.connect()
    con.execute("SET enable_progress_bar=false")
    con.execute("INSTALL httpfs")
    con.execute("LOAD httpfs")
    con.execute("SET s3_region='us-east-1'")
    con.execute("CREATE OR REPLACE SECRET anon (TYPE s3, PROVIDER config, KEY_ID '', SECRET '')")
    return con


def chembl_context(
    identifier: str | int,
    namespace: str = "cid",
    chembl_release: str = "chembl_37",
) -> dict[str, Any] | None:
    """ChEMBL context for a PubChem compound: resolves the ChEMBL ID via a
    live PubChem cross-reference (see xrefs.chembl_id), then looks up that
    compound's pref_name, max_phase and molregno in scigantic-chembl.

    Returns None if PubChem has no ChEMBL cross-reference for this
    identifier, or if that ChEMBL ID isn't in this mirror's release.
    """
    cid = _chembl_id(identifier, namespace=namespace)
    if cid is None:
        return None
    con = _connect()
    mol = f"s3://{_CHEMBL_BUCKET}/{chembl_release}/parquet/molecule_dictionary.parquet"
    row = con.execute(
        f"SELECT molregno, chembl_id, pref_name, max_phase FROM read_parquet('{mol}') "
        "WHERE chembl_id = ?",
        [cid],
    ).fetchone()
    if row is None:
        return {"chembl_id": cid, "molregno": None, "pref_name": None, "max_phase": None}
    return {"molregno": row[0], "chembl_id": row[1], "pref_name": row[2], "max_phase": row[3]}


def bindingdb_measurements(
    cid: int,
    bindingdb_release: str = "202608",
) -> "pd.DataFrame":
    """Binding measurements scigantic-bindingdb has for this PubChem CID.

    BindingDB's own measurements table already carries a pubchem_cid column
    -- this is a live filter against it, not a precomputed join, so it
    reflects the mirror as it stands with no separate build step.
    """
    con = _connect()
    path = f"s3://{_BINDINGDB_BUCKET}/{bindingdb_release}/parquet/measurements.parquet"
    return con.execute(
        f"SELECT reactant_set_id, ligand_smiles, ki_nm_value, ic50_nm_value, "
        f"kd_nm_value, ec50_nm_value, curation_source "
        f"FROM read_parquet('{path}') WHERE pubchem_cid = ?",
        [cid],
    ).df()
