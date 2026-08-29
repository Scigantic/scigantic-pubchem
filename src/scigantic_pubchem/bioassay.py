"""BioAssay data: assay summaries, bioactivity result tables in both
directions (per-assay and per-compound), AID/CID/SID cross-lookups, and a
streaming downloader for a whole assay's results.

PubChemPy's Assay/get_assays() only reaches PUG REST's `description`
operation, which returns the raw PC_AssayContainer record: a deeply nested
tree built to round-trip the depositor's original submission (protocol
text, full result-column schema, revision history, ...), not to be read
programmatically. Verified 2026-08-29 by reading pubchempy.py's source
directly (not assumed from its docs): the strings "concise", "assaysummary",
"summary", "cids", and "sids" never appear in it at all, under the `assay`
domain or otherwise. Those are the operations that actually carry the
tabular bioactivity data (AID/SID/CID/Activity Outcome/...) most callers
want, and `concise` specifically is PubChem's own purpose-built compact
format for that: the same data as the full record, as a table, meant for
bulk use. This module wraps those instead.

Every function here reads PUG REST's response table by column name rather
than positional index (see _parse_table): `concise` (per-assay input) and
`assaysummary` (per-compound input) return overlapping but not identical
column sets (assaysummary adds "Panel Member ID"; verified live 2026-08-29
against AID 1 and CID 2244), so a row's shape depends on which endpoint
produced it.
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from typing import TYPE_CHECKING, Any
from urllib.parse import quote

from . import _client
from .models import AssayResult, AssaySummary
from .resolve import Namespace

if TYPE_CHECKING:
    from collections.abc import Sequence
    from pathlib import Path

_POST_NAMESPACES = {"smiles", "inchi"}
_CHUNK_SIZE = 200
_MAX_CHUNK_WORKERS = 8


def _join_aids(aid: "int | str | Sequence[int | str]") -> str:
    if isinstance(aid, (int, str)):
        return str(aid)
    return ",".join(str(a) for a in aid)


def _parse_table(body: dict[str, Any]) -> list[dict[str, str]]:
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


def _row_to_result(row: dict[str, str]) -> AssayResult:
    return AssayResult(
        aid=int(row["AID"]),
        panel_member_id=_cell_str(row.get("Panel Member ID")),
        sid=_cell_int(row.get("SID")),
        cid=_cell_int(row.get("CID")),
        activity_outcome=_cell_str(row.get("Activity Outcome")),
        target_accession=_cell_str(row.get("Target Accession")),
        target_gene_id=_cell_str(row.get("Target GeneID")),
        activity_value_um=_cell_float(row.get("Activity Value [uM]")),
        activity_name=_cell_str(row.get("Activity Name")),
        assay_name=_cell_str(row.get("Assay Name")),
        assay_type=_cell_str(row.get("Assay Type")),
        pubmed_id=_cell_str(row.get("PubMed ID")),
        rnai=_cell_str(row.get("RNAi")),
    )


def _record_to_summary(rec: dict[str, Any]) -> AssaySummary:
    description = rec.get("Description")
    return AssaySummary(
        aid=int(rec["AID"]),
        name=rec.get("Name") or None,
        source_name=rec.get("SourceName") or None,
        description="\n".join(description) if description else None,
        method=rec.get("Method") or None,
        target=[
            {"accession": t.get("Accession", ""), "name": t.get("Name", "")}
            for t in rec.get("Target", [])
        ],
        cid_active=rec.get("CIDCountActive"),
        cid_inactive=rec.get("CIDCountInactive"),
        cid_total=rec.get("CIDCountAll"),
        sid_active=rec.get("SIDCountActive"),
        sid_inactive=rec.get("SIDCountInactive"),
        sid_total=rec.get("SIDCountAll"),
    )


def assay_summary(aid: int) -> AssaySummary | None:
    """The overview PubChem has on file for one assay: name, description,
    method, target(s) tested against, and active/inactive/total counts.
    None if the AID doesn't exist.
    """
    try:
        body = _client.request(f"/assay/aid/{aid}/summary/JSON")
    except _client.CompoundNotFoundError:
        return None
    records = body.get("AssaySummaries", {}).get("AssaySummary", [])
    if not records:
        return None
    return _record_to_summary(records[0])


def assay_results(aid: "int | str | Sequence[int | str]") -> list[AssayResult]:
    """Every (SID, CID, outcome) row PubChem recorded for one assay, or for
    several at once: PUG REST accepts a comma-separated AID list for
    `concise` in a single request (verified live 2026-08-29), so a small
    batch of related assays (e.g. every AID from aids_for_target()) costs
    one round trip rather than one per AID.

    Uses PUG REST's `concise` operation, the compact table built for this
    exact purpose rather than the full per-assay record. An assay with no
    live AIDs among those given returns an empty list rather than raising.
    """
    try:
        body = _client.request(f"/assay/aid/{_join_aids(aid)}/concise/JSON")
    except _client.CompoundNotFoundError:
        return []
    return [_row_to_result(row) for row in _parse_table(body)]


def assay_cids(aid: int) -> list[int]:
    """CIDs tested in one assay. Empty list if the AID doesn't exist."""
    try:
        body = _client.request(f"/assay/aid/{aid}/cids/JSON")
    except _client.CompoundNotFoundError:
        return []
    info = body.get("InformationList", {}).get("Information", [])
    return list(info[0].get("CID", [])) if info else []


def assay_sids(aid: int) -> list[int]:
    """SIDs (substances, PubChem's pre-standardization identifier) tested
    in one assay. Empty list if the AID doesn't exist."""
    try:
        body = _client.request(f"/assay/aid/{aid}/sids/JSON")
    except _client.CompoundNotFoundError:
        return []
    info = body.get("InformationList", {}).get("Information", [])
    return list(info[0].get("SID", [])) if info else []


def _assay_results_chunk(chunk: "Sequence[str]") -> list[AssayResult]:
    path = f"/compound/cid/{','.join(chunk)}/assaysummary/JSON"
    try:
        body = _client.request(path)
    except _client.CompoundNotFoundError:
        return []
    return [_row_to_result(row) for row in _parse_table(body)]


def compound_assay_results(identifier: str | int, namespace: Namespace = "cid") -> list[AssayResult]:
    """Every bioactivity row recorded across every assay that tested one
    compound: the reverse of assay_results(), keyed by compound instead of
    by assay. Uses PUG REST's `assaysummary` operation, which returns the
    same row shape as `concise` (plus a Panel Member ID column; see the
    module docstring) but selected by compound rather than by assay.
    """
    ident = quote(str(identifier), safe="")
    if namespace in _POST_NAMESPACES:
        path = f"/compound/{namespace}/assaysummary/JSON"
        params: dict[str, Any] | None = {namespace: str(identifier)}
        method = "POST"
    else:
        path = f"/compound/{namespace}/{ident}/assaysummary/JSON"
        params = None
        method = "GET"
    try:
        body = _client.request(path, params=params, method=method)
    except _client.CompoundNotFoundError:
        return []
    return [_row_to_result(row) for row in _parse_table(body)]


def compound_assay_results_many(cids: "Sequence[int | str]") -> dict[int, list[AssayResult]]:
    """Batched compound_assay_results() for CIDs, chunked and parallelized
    the same way resolve_many()/xrefs_many() are: verified live 2026-08-29
    that `assaysummary` accepts a comma-separated CID list in one request,
    so this chunks at 200 CIDs per call rather than one round trip per CID.
    A CID with no assay results at all is simply absent from the returned
    dict, the same way xrefs_many() omits a CID with no xrefs.
    """
    ids = [str(c) for c in cids]
    chunks = [ids[i : i + _CHUNK_SIZE] for i in range(0, len(ids), _CHUNK_SIZE)]

    if len(chunks) <= 1:
        rows = _assay_results_chunk(chunks[0]) if chunks else []
    else:
        rows = []
        with ThreadPoolExecutor(max_workers=min(_MAX_CHUNK_WORKERS, len(chunks))) as pool:
            for chunk_result in pool.map(_assay_results_chunk, chunks):
                rows.extend(chunk_result)

    results: dict[int, list[AssayResult]] = {}
    for row in rows:
        if row.cid is not None:
            results.setdefault(row.cid, []).append(row)
    return results


def aids_for_compound(identifier: str | int, namespace: Namespace = "cid") -> list[int]:
    """AIDs of every assay that tested one compound. Empty list if PubChem
    has no assay data for it."""
    ident = quote(str(identifier), safe="")
    if namespace in _POST_NAMESPACES:
        path = f"/compound/{namespace}/aids/JSON"
        params: dict[str, Any] | None = {namespace: str(identifier)}
        method = "POST"
    else:
        path = f"/compound/{namespace}/{ident}/aids/JSON"
        params = None
        method = "GET"
    try:
        body = _client.request(path, params=params, method=method)
    except _client.CompoundNotFoundError:
        return []
    info = body.get("InformationList", {}).get("Information", [])
    return list(info[0].get("AID", [])) if info else []


def aids_for_target(gene_symbol: str) -> list[int]:
    """AIDs of every assay run against a gene target, e.g. "EGFR". Empty
    list if PubChem has no assay data for that gene symbol."""
    ident = quote(gene_symbol, safe="")
    try:
        body = _client.request(f"/assay/target/genesymbol/{ident}/aids/JSON")
    except _client.CompoundNotFoundError:
        return []
    return list(body.get("IdentifierList", {}).get("AID", []))


def download_assay_results(
    aid: "int | str | Sequence[int | str]",
    dest: "Path | str",
    fmt: str = "csv",
) -> "Path":
    """Stream one or more assays' `concise` bioactivity table straight to
    a file, without ever holding the full response in memory the way
    assay_results() does.

    assay_results() is the right call for a typical assay (concise/JSON
    parses into a handful to a few thousand AssayResult records). This is
    for the other end of PubChem's size range: a large qHTS screen's
    concise table can run to hundreds of thousands of rows and tens of MB
    (measured: a 69,000-row assay's concise CSV, ~9.5MB) or more, and
    loading that into a Python list before a caller can even start
    processing it, only to often write it straight back out to disk or a
    dataframe anyway, is the exact kind of cost a dedicated download path
    should skip. fmt is "csv" (default, PubChem's own compact bulk format
    for this data) or "json".
    """
    if fmt not in ("csv", "json"):
        raise ValueError(f"fmt must be 'csv' or 'json', got {fmt!r}")
    path = f"/assay/aid/{_join_aids(aid)}/concise/{fmt.upper()}"
    return _client.stream_to_file(path, dest)
