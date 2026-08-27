"""Resolve an arbitrary compound reference (name, SMILES, InChIKey, CID) to
a canonical structure and identifiers, in one PUG REST round trip.

Verified 2026-08-27: PUG REST's property endpoint accepts a comma-separated
property list and returns them all in a single response -- and normalizes
some names in the response (a requested "CanonicalSMILES" comes back keyed
"SMILES", "IsomericSMILES" comes back "ConnectivitySMILES"). This module
reads by whichever key is actually present rather than assuming the
requested name matches the response key.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Literal
from urllib.parse import quote

from . import _client
from .models import Compound

if TYPE_CHECKING:
    from collections.abc import Sequence

Namespace = Literal["name", "cid", "smiles", "inchikey", "inchi", "formula"]

# What PubChem actually calls each field back, whichever alias is used.
_PROPERTIES = "CanonicalSMILES,IsomericSMILES,InChI,InChIKey,IUPACName,MolecularFormula,MolecularWeight,Title"
_POST_NAMESPACES = {"smiles", "inchi"}


def _record_to_compound(rec: dict[str, Any]) -> Compound:
    smiles = rec.get("SMILES") or rec.get("CanonicalSMILES") or rec.get("ConnectivitySMILES")
    mw = rec.get("MolecularWeight")
    return Compound(
        cid=int(rec["CID"]),
        title=rec.get("Title"),
        smiles=smiles,
        inchi=rec.get("InChI"),
        inchi_key=rec.get("InChIKey"),
        iupac_name=rec.get("IUPACName"),
        molecular_formula=rec.get("MolecularFormula"),
        molecular_weight=float(mw) if mw is not None else None,
    )


def resolve(identifier: str | int, namespace: Namespace = "name") -> Compound | None:
    """Resolve one identifier to a Compound. None if PubChem has no match.

    namespace: "name" (default), "cid", "smiles", "inchikey", "inchi", or
    "formula". CID may be passed as an int with namespace="cid" or left as
    the default and passed as a string CID.
    """
    ident = str(identifier)
    path = f"/compound/{namespace}/property/{_PROPERTIES}/JSON"
    try:
        if namespace in _POST_NAMESPACES:
            # SMILES/InChI can contain characters unsafe for a URL path.
            body = _client.request(path, params={namespace: ident}, method="POST")
        else:
            path = f"/compound/{namespace}/{quote(ident, safe='')}/property/{_PROPERTIES}/JSON"
            body = _client.request(path)
    except _client.CompoundNotFoundError:
        # A miss is a real, expected outcome (verified live: PUG REST
        # returns a proper 404 for an unmatched identifier, not an empty
        # 200), not something every caller should have to except-handle.
        return None

    records = body.get("PropertyTable", {}).get("Properties", [])
    if not records:
        return None
    return _record_to_compound(records[0])


def resolve_many(cids: "Sequence[int | str]") -> list[Compound]:
    """Resolve a batch of CIDs in as few PUG REST round trips as possible.

    PUG REST accepts a comma-separated CID list in a single request;
    verified live 2026-08-27 with a 3-CID batch. Chunked at 200 CIDs per
    request, comfortably under PubChem's practical request-size limits,
    rather than sent as one unbounded URL for a very long list.
    """
    results: list[Compound] = []
    chunk_size = 200
    ids = [str(c) for c in cids]
    for i in range(0, len(ids), chunk_size):
        chunk = ids[i : i + chunk_size]
        path = f"/compound/cid/{','.join(chunk)}/property/{_PROPERTIES}/JSON"
        body = _client.request(path)
        records = body.get("PropertyTable", {}).get("Properties", [])
        results.extend(_record_to_compound(r) for r in records)
    return results
