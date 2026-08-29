"""Shared parsing for PUG REST's `Table` response shape:
{"Table": {"Columns": {"Column": [...]}, "Row": [{"Cell": [...]}, ...]}}.

Every `concise`/`assaysummary` operation across the assay, compound, gene,
and protein domains returns this same shape, but not the same column set
(verified live 2026-08-29: assay-side `concise` has RNAi and Target
GeneID, gene-domain `concise` has Activity Qualifier but no Target GeneID,
protein-domain `concise` has Activity Qualifier but no Target Accession or
RNAi, and so on). Read by column name rather than positional index for
that reason; see bioassay.py and gene_protein.py, both of which parse this
into AssayResult rows.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from .models import AssayResult

if TYPE_CHECKING:
    from collections.abc import Sequence


def join_ids(ids: "int | str | Sequence[int | str]") -> str:
    """Comma-join one or more AIDs/CIDs/gene IDs/accessions for a batched
    PUG REST request, or pass a single one through unchanged."""
    if isinstance(ids, (int, str)):
        return str(ids)
    return ",".join(str(i) for i in ids)


def parse_table(body: dict[str, Any]) -> list[dict[str, str]]:
    table = body.get("Table", {})
    columns: list[str] = table.get("Columns", {}).get("Column", [])
    rows: list[dict[str, Any]] = table.get("Row", [])
    return [dict(zip(columns, row.get("Cell", []))) for row in rows]


def _cell_int(value: str | None) -> int | None:
    return int(value) if value else None


def _cell_float(value: str | None) -> float | None:
    if not value:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def _cell_str(value: str | None) -> str | None:
    return value if value else None


def row_to_result(row: dict[str, str]) -> AssayResult:
    return AssayResult(
        aid=int(row["AID"]),
        panel_member_id=_cell_str(row.get("Panel Member ID")),
        sid=_cell_int(row.get("SID")),
        cid=_cell_int(row.get("CID")),
        activity_outcome=_cell_str(row.get("Activity Outcome")),
        activity_qualifier=_cell_str(row.get("Activity Qualifier")),
        target_accession=_cell_str(row.get("Target Accession")),
        target_gene_id=_cell_str(row.get("Target GeneID")),
        activity_value_um=_cell_float(row.get("Activity Value [uM]")),
        activity_name=_cell_str(row.get("Activity Name")),
        assay_name=_cell_str(row.get("Assay Name")),
        assay_type=_cell_str(row.get("Assay Type")),
        pubmed_id=_cell_str(row.get("PubMed ID")),
        rnai=_cell_str(row.get("RNAi")),
    )
