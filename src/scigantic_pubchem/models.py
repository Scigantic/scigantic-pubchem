"""Typed records returned by this package. PubChemPy returns loosely-typed
Compound objects backed by a raw dict; these are plain, typed dataclasses."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Compound:
    cid: int
    title: str | None
    smiles: str | None
    inchi: str | None
    inchi_key: str | None
    iupac_name: str | None
    molecular_formula: str | None
    molecular_weight: float | None


@dataclass(frozen=True)
class AssaySummary:
    """The flat, machine-readable overview PUG REST's assay `summary`
    operation returns for one AID: name, description text, and
    active/inactive/total counts. Deliberately not the raw `description`
    operation's PC_AssayContainer record, a deeply nested tree built to
    round-trip the depositor's original submission (protocol XML, full
    result-column schema, ...) rather than to be read programmatically;
    `summary` carries the same name/description/protocol/comment text
    plus the counts, in a shape meant to be consumed directly."""

    aid: int
    name: str | None
    source_name: str | None
    description: str | None
    method: str | None
    target: list[dict[str, str]]
    cid_active: int | None
    cid_inactive: int | None
    cid_total: int | None
    sid_active: int | None
    sid_inactive: int | None
    sid_total: int | None


@dataclass(frozen=True)
class AssayResult:
    """One row of bioactivity data: a (SID, CID) tested in one assay (AID),
    with its outcome. The row shape PUG REST's `concise` (per-assay) and
    `assaysummary` (per-compound) operations both return, just keyed off
    a different identifier; see bioassay.py."""

    aid: int
    panel_member_id: str | None
    sid: int | None
    cid: int | None
    activity_outcome: str | None
    target_accession: str | None
    target_gene_id: str | None
    activity_value_um: float | None
    activity_name: str | None
    assay_name: str | None
    assay_type: str | None
    pubmed_id: str | None
    rnai: str | None
