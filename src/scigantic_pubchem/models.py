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
    with its outcome. The row shape PUG REST's `concise` (per-assay),
    `assaysummary` (per-compound), and the gene/protein domains' own
    `concise` (per-gene, per-protein) operations all return, just keyed off
    a different identifier; see bioassay.py and gene_protein.py.

    Not every source populates every field: verified live 2026-08-29 that
    the four sources return overlapping but not identical column sets
    (e.g. gene-domain concise has no Target GeneID, since the gene is
    already the query; protein-domain concise has no Target Accession, no
    RNAi, for the same reason). Absent columns read as None, the same as
    a present-but-empty cell.
    """

    aid: int
    panel_member_id: str | None
    sid: int | None
    cid: int | None
    activity_outcome: str | None
    activity_qualifier: str | None
    target_accession: str | None
    target_gene_id: str | None
    activity_value_um: float | None
    activity_name: str | None
    assay_name: str | None
    assay_type: str | None
    pubmed_id: str | None
    rnai: str | None


@dataclass(frozen=True)
class DoseResponsePoint:
    """One titration point from a qHTS curve: a tested concentration and
    the raw percent response PubChem recorded there, before any curve
    fitting. Present even for a compound with no fitted curve at all --
    an Inactive row still carries its raw per-concentration readings,
    verified live 2026-08-31 against AID 1851."""

    concentration_um: float
    response_percent: float | None


@dataclass(frozen=True)
class RawAssayResult:
    """One (SID, panel member) row from an assay's full, un-curated Data
    Table -- PUG REST's plain `/assay/aid/{aid}/CSV` operation, not
    `concise`. `concise`/`assaysummary` (AssayResult above) strip a qHTS
    screen down to outcome plus a single fitted potency; this is the
    layer underneath: the raw per-concentration response PubChem's
    depositor actually submitted, including for compounds that never
    cleared curve-fitting into a potency at all.

    Verified live 2026-08-31 against AID 1851 (a 5-isoform CYP inhibition
    qHTS panel) that this is genuinely a different, much wider column set
    than `concise`'s fixed 12 columns, not an extension of it. Only
    PubChem's reserved PUBCHEM_* columns, the standard NCGC/NCATS qHTS
    curve-fit vocabulary (Potency, Curve_Description, Fit_*,
    Max_Response), and the Panel_* columns (present on multi-target panel
    assays) are modelled as named fields here -- the remainder of that
    AID's columns ("Inhibition Observed", "Approved Drug", "Collection",
    "Compound QC", ...) are depositor-specific commentary, not part of
    any fixed schema, so they land in `extra` rather than being hardcoded
    from one assay's shape. See dose_response()/download_dose_response()
    in bioassay.py.
    """

    aid: int
    sid: int
    cid: int | None
    smiles: str | None
    activity_outcome: str | None
    activity_score: int | None
    activity_url: str | None
    comment: str | None
    panel_id: int | None
    panel_name: str | None  # the mnemonic, e.g. "p450-cyp1a2" -- verified live 2026-08-31
    panel_target: str | None  # the protein accession, e.g. "NP_000752.2" -- easy to mix these two up
    potency_um: float | None
    curve_description: str | None
    fit_logac50: float | None
    fit_hillslope: float | None
    fit_r2: float | None
    fit_infinite_activity: float | None
    fit_zero_activity: float | None
    fit_curveclass: str | None
    excluded_points: str | None
    max_response: float | None
    dose_response: tuple[DoseResponsePoint, ...]
    extra: dict[str, str]


@dataclass(frozen=True)
class GeneInfo:
    gene_id: int
    symbol: str | None
    name: str | None
    taxonomy_id: int | None
    taxonomy: str | None
    description: str | None
    synonyms: list[str]


@dataclass(frozen=True)
class ProteinInfo:
    accession: str
    name: str | None
    taxonomy_id: int | None
    taxonomy: str | None
    synonyms: list[str]
