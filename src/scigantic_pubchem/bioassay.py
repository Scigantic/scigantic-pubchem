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
than positional index (see _tables.parse_table): `concise` (per-assay
input) and `assaysummary` (per-compound input) return overlapping but not
identical column sets (assaysummary adds "Panel Member ID"; verified live
2026-08-29 against AID 1 and CID 2244), so a row's shape depends on which
endpoint produced it.

dose_response()/download_dose_response() reach past `concise` on
purpose, into the plain (non-concise) `/CSV` operation: PubChem's actual
raw Data Table, with every tested concentration's raw response, not the
outcome+potency summary `concise` reduces a qHTS screen to. Verified live
2026-08-31 against AID 1851 that this is a real, different gap: `concise`
carries no Max_Response, no curve-fit parameters, and no per-concentration
columns at all, even though PubChem's own public record for that AID
carries all of it, including for compounds with no fitted curve. This is
a structurally different operation from the rest of the module (a
different response format, a much wider and assay-dependent column set,
and its own hard identifier cap), so it gets its own CSV parsing rather
than reusing _tables.parse_table's JSON `Table` shape.
"""

from __future__ import annotations

import csv
import io
import json
import re
import warnings
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import TYPE_CHECKING, Any
from urllib.parse import quote

from . import _client
from ._tables import _cell_float, _cell_int, _cell_str, join_ids, parse_table, row_to_result
from .models import AssayResult, AssaySummary, DoseResponsePoint, RawAssayResult
from .resolve import Namespace

if TYPE_CHECKING:
    from collections.abc import Sequence

_POST_NAMESPACES = {"smiles", "inchi"}
_CHUNK_SIZE = 200
_MAX_CHUNK_WORKERS = 8

# dose_response()/download_dose_response(): PUG REST's raw (non-concise)
# assay Data Table, a structurally different operation from everything
# else in this module (see their docstrings for what it adds).
_CONCENTRATION_COLUMN_RE = re.compile(r"^Activity at ([\d.eE+-]+) uM$")
_RAW_KNOWN_COLUMNS = frozenset(
    {
        "PUBCHEM_RESULT_TAG",
        "PUBCHEM_SID",
        "PUBCHEM_CID",
        "PUBCHEM_EXT_DATASOURCE_SMILES",
        "PUBCHEM_ACTIVITY_OUTCOME",
        "PUBCHEM_ACTIVITY_SCORE",
        "PUBCHEM_ACTIVITY_URL",
        "PUBCHEM_ASSAYDATA_COMMENT",
        "Potency",
        "Curve_Description",
        "Fit_LogAC50",
        "Fit_HillSlope",
        "Fit_R2",
        "Fit_InfiniteActivity",
        "Fit_ZeroActivity",
        "Fit_CurveClass",
        "Excluded_Points",
        "Max_Response",
        "Panel ID",
        "Panel Name",
        "Panel Target",
    }
)
# dose_response(): PUG REST hard-caps the underlying operation at 10,000
# identifiers per request ("Assay record retrieval is limited to 10000
# SIDs", verified live 2026-08-31), but that cap is far past the point
# this is still a reasonable single call -- see the function's docstring.
_DOSE_RESPONSE_MAX_IDS = 200
# download_dose_response(): matches _CHUNK_SIZE above in name only. The
# right value here is much more conservative than that number implies;
# see the function's docstring for the live measurement behind it.
_DOSE_RESPONSE_CHUNK_SIZE = 200
_DOSE_RESPONSE_TIMEOUT = 180.0


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
        body = _client.request(f"/assay/aid/{join_ids(aid)}/concise/JSON")
    except _client.CompoundNotFoundError:
        return []
    return [row_to_result(row) for row in parse_table(body)]


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
    return [row_to_result(row) for row in parse_table(body)]


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
    return [row_to_result(row) for row in parse_table(body)]


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
    path = f"/assay/aid/{join_ids(aid)}/concise/{fmt.upper()}"
    return _client.stream_to_file(path, dest)


def _parse_dose_response_row(row: dict[str, str], aid: int, conc_columns: list[tuple[str, float]]) -> RawAssayResult:
    conc_names = {name for name, _ in conc_columns}
    dose_response = tuple(
        DoseResponsePoint(concentration_um=conc, response_percent=_cell_float(row.get(name)))
        for name, conc in conc_columns
    )
    extra = {k: v for k, v in row.items() if v and k not in _RAW_KNOWN_COLUMNS and k not in conc_names}
    return RawAssayResult(
        aid=aid,
        sid=int(row["PUBCHEM_SID"]),
        cid=_cell_int(row.get("PUBCHEM_CID")),
        smiles=_cell_str(row.get("PUBCHEM_EXT_DATASOURCE_SMILES")),
        activity_outcome=_cell_str(row.get("PUBCHEM_ACTIVITY_OUTCOME")),
        activity_score=_cell_int(row.get("PUBCHEM_ACTIVITY_SCORE")),
        activity_url=_cell_str(row.get("PUBCHEM_ACTIVITY_URL")),
        comment=_cell_str(row.get("PUBCHEM_ASSAYDATA_COMMENT")),
        panel_id=_cell_int(row.get("Panel ID")),
        panel_name=_cell_str(row.get("Panel Name")),
        panel_target=_cell_str(row.get("Panel Target")),
        potency_um=_cell_float(row.get("Potency")),
        curve_description=_cell_str(row.get("Curve_Description")),
        fit_logac50=_cell_float(row.get("Fit_LogAC50")),
        fit_hillslope=_cell_float(row.get("Fit_HillSlope")),
        fit_r2=_cell_float(row.get("Fit_R2")),
        fit_infinite_activity=_cell_float(row.get("Fit_InfiniteActivity")),
        fit_zero_activity=_cell_float(row.get("Fit_ZeroActivity")),
        fit_curveclass=_cell_str(row.get("Fit_CurveClass")),
        excluded_points=_cell_str(row.get("Excluded_Points")),
        max_response=_cell_float(row.get("Max_Response")),
        dose_response=dose_response,
        extra=extra,
    )


def _parse_dose_response_csv(text: str, aid: int) -> list[RawAssayResult]:
    reader = csv.DictReader(io.StringIO(text))
    conc_columns = [
        (name, float(m.group(1))) for name in (reader.fieldnames or []) if (m := _CONCENTRATION_COLUMN_RE.match(name))
    ]
    results = []
    for row in reader:
        tag = row.get("PUBCHEM_RESULT_TAG")
        if not tag or not tag.isdigit():
            continue  # metadata rows (RESULT_TYPE/RESULT_DESCR/RESULT_UNIT/...), not data
        results.append(_parse_dose_response_row(row, aid, conc_columns))
    return results


def dose_response(
    aid: int,
    sids: "Sequence[int | str] | None" = None,
    cids: "Sequence[int | str] | None" = None,
) -> list[RawAssayResult]:
    """The full, un-curated per-concentration qHTS readout for a bounded
    set of compounds in one assay: PUG REST's plain `/CSV` operation (not
    `concise`), the raw depositor Data Table rather than PubChem's
    outcome+potency summary. Includes Max_Response and every tested
    concentration's raw percent response, present even for a compound
    PubChem marked Inactive with no fitted curve at all -- verified live
    2026-08-31 against AID 1851 (a 5-isoform CYP inhibition qHTS panel,
    17,143 compounds).

    Pass sids or cids, not both. Capped at _DOSE_RESPONSE_MAX_IDS per
    call rather than chunked automatically: this operation is
    dramatically slower than concise/assaysummary at any real scale
    (measured live against AID 1851: 250 SIDs returned in under a
    second, 2000 SIDs took 94s against the same server), and PUG REST
    itself hard-caps it at 10,000 identifiers per request regardless
    ("Assay record retrieval is limited to 10000 SIDs", verified live).
    For a whole assay, use download_dose_response(), which chunks
    conservatively and writes to disk as it goes rather than risking one
    large, slow, unresumable request.
    """
    if (sids is None) == (cids is None):
        raise ValueError("pass exactly one of sids or cids")
    ids = [str(i) for i in (sids if sids is not None else cids)]  # type: ignore[union-attr]
    if len(ids) > _DOSE_RESPONSE_MAX_IDS:
        raise ValueError(
            f"{len(ids)} identifiers exceeds the {_DOSE_RESPONSE_MAX_IDS}-per-call cap of dose_response(); "
            "use download_dose_response() for a whole assay"
        )
    if not ids:
        return []
    param_name = "sid" if sids is not None else "cid"
    try:
        text = _client.request_text(
            f"/assay/aid/{aid}/CSV",
            params={param_name: ",".join(ids)},
            method="POST",
        )
    except _client.CompoundNotFoundError:
        return []
    return _parse_dose_response_csv(text, aid)


def _append_dose_response_chunk(dest: Path, text: str, include_preamble: bool) -> None:
    """Append one chunk's response to dest, keeping the header and
    PUG REST's own RESULT_TYPE/RESULT_DESCR/RESULT_UNIT/... metadata
    preamble only when include_preamble (the first chunk): every chunk's
    response repeats that preamble in full, not just the header, so a
    naive "skip the first line only" join would leave a duplicate
    4-line preamble sitting after every chunk boundary in the combined
    file -- caught by this module's own tests, not by inspection. Uses
    csv.reader/csv.writer rather than string splitting so a field
    containing an embedded comma or quote round-trips correctly.
    """
    rows = list(csv.reader(io.StringIO(text)))
    if not rows:
        return
    header, body_rows = rows[0], rows[1:]
    with dest.open("a", newline="") as f:
        writer = csv.writer(f)
        if include_preamble:
            writer.writerow(header)
            for row in body_rows:
                writer.writerow(row)
        else:
            for row in body_rows:
                if row and row[0].isdigit():
                    writer.writerow(row)


def download_dose_response(
    aid: int,
    dest: "Path | str",
    sids: "Sequence[int] | None" = None,
    chunk_size: int = _DOSE_RESPONSE_CHUNK_SIZE,
    resume: bool = True,
) -> "Path":
    """Stream one assay's full, un-curated dose-response Data Table to a
    CSV file: every tested concentration's raw response for every
    compound, the layer underneath download_assay_results()'s `concise`
    outcome+potency summary (see dose_response()'s docstring for what
    that adds and why it's capped per call).

    Chunked at chunk_size SIDs per request. The default is conservative
    on purpose: measured live against AID 1851 (17,143 SIDs), this
    operation's cost does not scale linearly with request size -- 250
    SIDs returned in under a second, 2000 SIDs took 94s against the same
    server -- so pulling a whole large panel assay is realistically many
    minutes of wall time even chunked, well before PUG REST's separate
    hard cap of 10,000 SIDs per request ever applies. Sequential, not
    parallelized the way compound_assay_results_many() is: firing
    several of these already-expensive requests at PubChem concurrently
    risks compounding the slowdown just measured, against a live public
    NIH service other callers share.

    Because a full pull can run to many minutes across many chunks, this
    is resumable by default: progress lives in a sidecar
    `{dest}.progress.json` (assay identity, chunk size, and total SID
    count, so a resume attempt against different arguments is detected
    and discarded with a warning rather than silently misapplied), plus
    the exact byte offset dest was known-good at. Re-running the same
    call after an interruption truncates dest back to that offset and
    continues from the last confirmed chunk, so a crash mid-write of one
    chunk can never leave a duplicated or corrupt file the way naively
    trusting "chunk N was requested" without confirming it was written
    would. resume=False starts over and overwrites dest.

    sids defaults to the assay's full SID list (assay_sids(aid)); pass an
    explicit subset to pull only part of a large assay.
    """
    dest = Path(dest)
    all_sids = list(sids) if sids is not None else assay_sids(aid)
    if not all_sids:
        raise ValueError(f"no SIDs found for AID {aid}")
    chunks = [all_sids[i : i + chunk_size] for i in range(0, len(all_sids), chunk_size)]
    fingerprint = {"aid": aid, "chunk_size": chunk_size, "total_sids": len(all_sids), "n_chunks": len(chunks)}

    progress_path = dest.with_name(dest.name + ".progress.json")
    state: dict[str, Any] | None = None
    if resume and progress_path.exists():
        try:
            saved = json.loads(progress_path.read_text())
        except (json.JSONDecodeError, OSError):
            saved = None
        if saved and saved.get("fingerprint") == fingerprint and dest.exists():
            state = saved
        elif saved is not None:
            warnings.warn(
                f"discarding stale resume state at {progress_path} (assay or chunking changed); starting over",
                stacklevel=2,
            )

    if state is None:
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text("")
        state = {"fingerprint": fingerprint, "done": 0, "byte_offset": 0}
        progress_path.write_text(json.dumps(state))

    with dest.open("r+b") as f:
        f.truncate(state["byte_offset"])

    for i in range(state["done"], len(chunks)):
        text = _client.request_text(
            f"/assay/aid/{aid}/CSV",
            params={"sid": ",".join(str(s) for s in chunks[i])},
            method="POST",
            timeout=_DOSE_RESPONSE_TIMEOUT,
        )
        _append_dose_response_chunk(dest, text, include_preamble=(i == 0))
        state = {"fingerprint": fingerprint, "done": i + 1, "byte_offset": dest.stat().st_size}
        progress_path.write_text(json.dumps(state))

    progress_path.unlink(missing_ok=True)
    return dest
