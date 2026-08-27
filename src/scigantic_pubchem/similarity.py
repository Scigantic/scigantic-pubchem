"""Structure similarity and substructure search, live against PubChem's
entire corpus. No local fingerprint database needed.

scigantic-chembl's similar_compounds()/substructure_search() precompute
Morgan/Pattern fingerprints once and search them locally, because ChEMBL
has no live search API of its own. Fast, but bounded to the 1.68M
ChEMBL compounds that carry a comparable potency measurement. PubChem
runs this search itself, live, over its full ~120M-compound corpus
(fastsimilarity_2d / fastsubstructure), verified sub-second for a typical
query 2026-08-27. PubChemPy exposes the same PUG REST capability, but only
as a raw searchtype="similarity"/"substructure" parameter to its generic
get_compounds(), not a named, documented function.

Both can respond asynchronously for an expensive query (PubChem returns a
ListKey to poll rather than blocking the connection); see
_client.request_search() for that protocol, verified against PubChemPy's
own handling of it rather than assumed from docs.
"""

from __future__ import annotations

from . import _client
from .models import Compound
from .resolve import resolve_many


def similar_compounds(
    smiles: str,
    threshold: int = 90,
    max_records: int = 100,
    resolve: bool = True,
) -> list[Compound] | list[int]:
    """CIDs (or full Compound records) similar to a query structure.

    threshold is PubChem's 2D Tanimoto similarity cutoff, 0-100 (their
    default is 90). resolve=True (the default) fetches full Compound
    records for every match via resolve_many(), one extra batched round
    trip; resolve=False returns bare CIDs if that's all that's needed.

    PubChem does not return a similarity score alongside the matches, only
    the list of CIDs. Compute scored ranking yourself from the returned
    structures if you need it.
    """
    body = _client.request_search(
        "/compound/fastsimilarity_2d/smiles/cids/JSON",
        params={"smiles": smiles, "Threshold": threshold, "MaxRecords": max_records},
    )
    cids: list[int] = body.get("IdentifierList", {}).get("CID", [])
    return resolve_many(cids) if resolve else cids


def substructure_search(
    query: str,
    query_type: str = "smiles",
    max_records: int = 100,
    resolve: bool = True,
) -> list[Compound] | list[int]:
    """CIDs (or full Compound records) that contain a query substructure.

    query_type is "smiles" (default) or "smarts". PubChem runs these as
    genuinely different endpoints with different matching semantics
    (verified live: the same ring pattern given as SMILES vs the
    equivalent SMARTS returned overlapping but not identical CID lists),
    not two names for the same thing, so this does not try to guess which
    one a string is meant to be. resolve=True (the default) fetches full
    Compound records via resolve_many(); resolve=False returns bare CIDs.
    """
    if query_type not in ("smiles", "smarts"):
        raise ValueError(f"query_type must be 'smiles' or 'smarts', got {query_type!r}")
    body = _client.request_search(
        f"/compound/fastsubstructure/{query_type}/cids/JSON",
        params={query_type: query, "MaxRecords": max_records},
    )
    cids: list[int] = body.get("IdentifierList", {}).get("CID", [])
    return resolve_many(cids) if resolve else cids
