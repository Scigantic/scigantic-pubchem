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

from concurrent.futures import ThreadPoolExecutor
from typing import TYPE_CHECKING

from . import _client
from .models import Compound
from .resolve import resolve_many

if TYPE_CHECKING:
    from collections.abc import Sequence

# Lower than resolve_many()'s/compound_assay_results_many()'s 8: those chunk
# plain property/table lookups, while fastsimilarity_2d/fastsubstructure are
# PubChem's own "search" operation class, run against a corpus of ~120M
# compounds and documented as able to respond asynchronously under load (see
# request_search()). Measured live 2026-09-04: 8 concurrent
# fastsimilarity_2d searches against distinct drug-like queries, all
# resolved synchronously with no throttling and no errors -- but that is 8
# queries, not the hundreds a real batch (e.g. screening a few hundred test
# compounds against PubChem's live corpus) will actually run. Kept
# conservative pending a measurement at that scale; every request still
# passes through the same shared rate limiter regardless of max_workers, so
# this bounds how many searches are in flight at once, not the overall
# request rate.
_MAX_SEARCH_WORKERS = 5


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


def similar_compounds_many(
    smiles_list: "Sequence[str]",
    threshold: int = 90,
    max_records: int = 100,
    resolve: bool = True,
    max_workers: int = _MAX_SEARCH_WORKERS,
) -> "dict[str, list[Compound] | list[int]]":
    """similar_compounds() for many query structures at once: "find every
    close neighbor of each of these N compounds" as one call instead of a
    hand-rolled loop, the composable half of "find neighbors, then join
    assays" together with compound_assay_results_many().

    Each query still costs PubChem its own live fastsimilarity_2d search
    (there is no batched multi-query form of this endpoint on PubChem's
    side, unlike resolve_many()'s comma-separated CID list), but the
    searches run concurrently through a thread pool bounded by
    max_workers rather than one at a time -- the actual fix for N
    independent live searches paying N full round trips in sequence.
    Every request, from every worker, still passes through the same
    shared rate limiter as any other call in this package, so this cannot
    exceed PubChem's documented request rate regardless of max_workers;
    see _MAX_SEARCH_WORKERS's comment for why the default here is lower
    than resolve_many()'s.

    Returns a dict keyed by the input SMILES string (an empty list, not a
    missing key, if a query returns no hits), one query's matches at a
    time -- unlike similar_compounds() itself, this always returns bare
    CIDs or Compound records per query, resolve applies uniformly to all
    of them. A SMILES repeated in smiles_list collapses to one key, since
    an identical query always returns an identical result; pass distinct
    queries if that collision matters to a caller.

    If any one query fails after _client's own retries (a malformed
    SMILES, a real PubChem error), that exception propagates and aborts
    the whole batch, the same way a single failing chunk aborts
    resolve_many() today -- this does not add per-query error suppression
    that swallows a real failure into a silently empty result.
    """
    queries = list(smiles_list)
    if not queries:
        return {}

    def _search_one(smiles: str) -> "list[Compound] | list[int]":
        return similar_compounds(smiles, threshold=threshold, max_records=max_records, resolve=resolve)

    with ThreadPoolExecutor(max_workers=min(max_workers, len(queries))) as pool:
        matches = list(pool.map(_search_one, queries))
    return dict(zip(queries, matches))


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
