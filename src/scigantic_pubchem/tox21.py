"""The Tox21 Data Challenge panel: 12 qHTS nuclear-receptor and
stress-response assays (NR-AR, NR-AR-LBD, NR-AhR, NR-Aromatase, NR-ER,
NR-ER-LBD, NR-PPAR-gamma, SR-ARE, SR-ATAD5, SR-HSE, SR-MMP, SR-p53) -- the
same 12 columns MoleculeNet/DeepChem's tox21.csv exposes for ML
benchmarking.

The AID for each endpoint below was verified live 2026-08-29 three
independent ways, not taken from a single source: cross-checked against
Huang et al. 2016 (Frontiers in Environmental Science 3:85, Table 1),
confirmed each AID's own PubChem `summary` name names the right
target/pathway (e.g. AID 743122's name is literally "...activate the aryl
hydrocarbon receptor (AhR) signaling pathway"), and confirmed all 12 share
the same 8,099- or 7,329-CID count -- the same Tox21 10K compound library
screened across the whole panel, not 12 coincidentally similar assays.

This reconstructs the panel live from PubChem's own qHTS `concise` tables.
It is not a byte-for-bit copy of NCATS's original challenge distribution,
which did its own compound-level SMILES canonicalization and
train/leaderboard/test split -- use this for live access to the same 12
assays, and the original challenge SDF if bit-for-bit reproduction of a
published benchmark's rows matters.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from .bioassay import assay_results
from .models import AssayResult

if TYPE_CHECKING:
    from collections.abc import Sequence

TOX21_ENDPOINTS: dict[str, int] = {
    "NR-AR": 743040,
    "NR-AR-LBD": 743053,
    "NR-AhR": 743122,
    "NR-Aromatase": 743139,
    "NR-ER": 743079,
    "NR-ER-LBD": 743077,
    "NR-PPAR-gamma": 743140,
    "SR-ARE": 743219,
    "SR-ATAD5": 720516,
    "SR-HSE": 743228,
    "SR-MMP": 720637,
    "SR-p53": 720552,
}

_AID_TO_ENDPOINT = {aid: name for name, aid in TOX21_ENDPOINTS.items()}


def _resolve_endpoints(endpoints: "Sequence[str] | None") -> list[str]:
    names = list(endpoints) if endpoints is not None else list(TOX21_ENDPOINTS)
    unknown = [n for n in names if n not in TOX21_ENDPOINTS]
    if unknown:
        raise ValueError(f"Unknown Tox21 endpoint(s): {unknown}. Valid: {list(TOX21_ENDPOINTS)}")
    return names


def tox21_results(endpoints: "Sequence[str] | None" = None) -> list[AssayResult]:
    """Raw per-(SID, CID) bioactivity rows for one or more Tox21 panel
    endpoints, by name (e.g. "NR-AhR", "SR-p53"), or the full 12-assay
    panel if endpoints is None. Thin wrapper over assay_results(), which
    already batches multiple AIDs into one PUG REST request -- the full
    panel is one request, but its response runs to roughly 120,000 rows
    and 45MB (verified live 2026-08-29), so pass specific endpoints if you
    only need a few.
    """
    names = _resolve_endpoints(endpoints)
    return assay_results([TOX21_ENDPOINTS[name] for name in names])


def _consolidate(outcomes: "Sequence[str | None]") -> "int | None":
    """One compound's replicate rows for one endpoint -> a single binary
    label: Active if any replicate resolved Active, Inactive if every
    resolved replicate was Inactive, None if none resolved (Inconclusive
    or unscored only). A compound can carry more than one row per
    endpoint: verified live 2026-08-29 against SR-p53 (AID 720552) alone,
    1,873 of its CIDs have more than one SID, i.e. more than one physical
    sample tested."""
    if any(o == "Active" for o in outcomes):
        return 1
    if any(o == "Inactive" for o in outcomes):
        return 0
    return None


def tox21_matrix(endpoints: "Sequence[str] | None" = None) -> dict[int, dict[str, "int | None"]]:
    """The wide multi-task label matrix ML work actually wants: one row per
    CID, one column per requested Tox21 endpoint (default: all 12), valued
    1 (active), 0 (inactive), or None (tested but never resolved, or not
    tested in that endpoint at all -- Tox21 is a canonical benchmark with
    real missing labels, not every compound was screened against every
    endpoint). This assembly step, one label per compound per endpoint
    with replicate rows consolidated, is what raw PubChem access doesn't
    hand you; see _consolidate for the exact rule.
    """
    names = _resolve_endpoints(endpoints)
    rows = tox21_results(names)

    by_cid: dict[int, dict[str, list[str | None]]] = {}
    for row in rows:
        if row.cid is None:
            continue
        endpoint = _AID_TO_ENDPOINT[row.aid]
        by_cid.setdefault(row.cid, {}).setdefault(endpoint, []).append(row.activity_outcome)

    return {
        cid: {name: (_consolidate(outcomes[name]) if name in outcomes else None) for name in names}
        for cid, outcomes in by_cid.items()
    }
