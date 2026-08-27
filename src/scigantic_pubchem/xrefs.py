"""Cross-references PubChem already knows about, via PUG REST's xrefs
endpoint. Reachable through PubChemPy's low-level request()/get()
functions but not wrapped in a named, documented convenience there.

Verified live 2026-08-27: CID 2244 (aspirin)'s RegistryID xrefs include
"CHEMBL25" directly, alongside CAS numbers, ChEBI IDs, and many other
database identifiers PubChem's depositors have registered. This means a
ChEMBL cross-reference does not need a precomputed bridge table the way
scigantic-bindingdb's does for BindingDB. It can be read live, and it is
never stale, since it reflects whatever PubChem currently has on file.
"""

from __future__ import annotations

from urllib.parse import quote

from . import _client


def xrefs(identifier: str | int, xref_type: str = "RegistryID", namespace: str = "cid") -> list[str]:
    """Raw cross-reference IDs of one type PubChem has on file for a CID.

    xref_type is one of PubChem's own types: RegistryID (the broadest,
    covers most external database IDs, including ChEMBL), RN (CAS Registry
    Number), PubMedID, PatentID, and others documented at
    https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest#section=xrefs.
    """
    ident = quote(str(identifier), safe="")
    path = f"/compound/{namespace}/{ident}/xrefs/{xref_type}/JSON"
    try:
        body = _client.request(path)
    except _client.CompoundNotFoundError:
        # Verified live: PUG REST 404s when there are zero xrefs of this
        # type, not an empty 200. A real, expected outcome.
        return []
    info = body.get("InformationList", {}).get("Information", [])
    if not info:
        return []
    return list(info[0].get(xref_type, []))


def chembl_id(identifier: str | int, namespace: str = "cid") -> str | None:
    """The ChEMBL ID PubChem has on file for this compound, if any.

    A live lookup through PubChem's own RegistryID cross-references, not a
    precomputed table (see the module docstring). Returns the first match;
    in the rare case PubChem lists more than one CHEMBL id, use xrefs()
    directly to see all of them.
    """
    for ref in xrefs(identifier, xref_type="RegistryID", namespace=namespace):
        if ref.startswith("CHEMBL"):
            return ref
    return None
