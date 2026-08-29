"""Cross-references PubChem already knows about, via PUG REST's xrefs
endpoint. Reachable through PubChemPy's low-level request()/get()
functions but not wrapped in a named, documented convenience there.

Verified live 2026-08-27: CID 2244 (aspirin)'s RegistryID xrefs include
"CHEMBL25" directly, alongside CAS numbers, ChEBI IDs, and many other
database identifiers PubChem's depositors have registered. This means a
ChEMBL cross-reference does not need a precomputed bridge table the way
scigantic-bindingdb's does for BindingDB. It can be read live, and it is
never stale, since it reflects whatever PubChem currently has on file.

A different xrefs type, MMDBID, gives the same live-lookup treatment to
structural data: see pdb_structures() below.
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from typing import TYPE_CHECKING
from urllib.parse import quote

from . import _client

if TYPE_CHECKING:
    from collections.abc import Sequence

_CHUNK_SIZE = 200
_MAX_CHUNK_WORKERS = 8


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


def pdb_structures(identifier: str | int, namespace: str = "cid") -> list[int]:
    """MMDB IDs of deposited structures PubChem has on file with this
    compound bound as a ligand: NCBI's Molecular Modeling Database, its
    own mirror of PDB structure data with its own numeric ID space.

    Not the 4-character PDB ID itself. Verified live 2026-08-29 against
    PUG REST's own documented xrefs types: MMDBID is valid, "PDBID" is
    not (400 PUGREST.BadRequest, "Invalid xrefs type") -- PUG REST has no
    xrefs type that returns a PDB ID directly, only this MMDB proxy, and
    mapping an MMDB ID to its PDB ID needs a separate NCBI service PUG
    REST doesn't expose, out of scope for a package that only ever calls
    pubchem.ncbi.nlm.nih.gov/rest/pug. Each ID here is still a real,
    directly usable pointer to a deposited structure, viewable at
    https://www.ncbi.nlm.nih.gov/Structure/mmdb/mmdbsrv.cgi?uid={id}
    (confirmed live: CID 5291/imatinib resolves to 27 structures, CID
    2244/aspirin to 8; a compound never crystallized as a ligand, like
    CID 1/methane, returns an empty list rather than raising).
    """
    return [int(ref) for ref in xrefs(identifier, xref_type="MMDBID", namespace=namespace)]


def pdb_structures_many(cids: "Sequence[int | str]") -> dict[int, list[int]]:
    """Batched pdb_structures() for CIDs; see xrefs_many() for the
    round-trip savings. Every requested CID gets an entry, an empty list
    if it has no structures on file, the same way chembl_ids_many() gives
    every CID an entry regardless of xrefs_many()'s own chunk-dependent
    presence (see its docstring) -- this loops over the input CIDs
    explicitly rather than trusting which keys came back.
    """
    all_refs = xrefs_many(cids, xref_type="MMDBID")
    return {int(cid): [int(ref) for ref in all_refs.get(int(cid), [])] for cid in cids}


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


def _xrefs_chunk(chunk: "Sequence[str]", xref_type: str) -> dict[int, list[str]]:
    path = f"/compound/cid/{','.join(chunk)}/xrefs/{xref_type}/JSON"
    try:
        body = _client.request(path)
    except _client.CompoundNotFoundError:
        return {}
    info = body.get("InformationList", {}).get("Information", [])
    return {int(rec["CID"]): list(rec.get(xref_type, [])) for rec in info if "CID" in rec}


def xrefs_many(cids: "Sequence[int | str]", xref_type: str = "RegistryID") -> dict[int, list[str]]:
    """Cross-reference IDs of one type for a batch of CIDs, in as few PUG
    REST round trips as possible.

    xrefs() only takes one identifier, so cross-referencing N compounds
    meant N round trips with no way around it. Verified live 2026-08-27
    that PUG REST's xrefs endpoint accepts a comma-separated CID list the
    same way the property endpoint resolve_many() uses does (confirmed at
    200 CIDs in one request, not just a small batch), so this chunks and
    parallelizes the same way.

    Whether a CID with zero xrefs of this type appears in the result at
    all depends on its chunk-mates, not just on it: verified live
    2026-08-29 that a chunk where every queried CID has zero matches gets
    a genuine 404 (caught by _xrefs_chunk, contributing nothing to the
    result -- those CIDs are absent), but a chunk with at least one match
    returns 200 with an explicit empty-list entry for every CID in that
    chunk that had none. A caller that needs "was this CID even asked
    about" for a CID that might land in either kind of chunk should treat
    a missing key and an empty-list value the same way, not rely on one
    implying the other.

    Unlike xrefs(), this is CID-only: PUG REST does not support
    comma-separated batching for the name/smiles/inchikey namespaces
    (verified live: a multi-name request 404s), so there's no batched
    counterpart to build for those.
    """
    ids = [str(c) for c in cids]
    chunks = [ids[i : i + _CHUNK_SIZE] for i in range(0, len(ids), _CHUNK_SIZE)]

    if len(chunks) <= 1:
        return _xrefs_chunk(chunks[0], xref_type) if chunks else {}

    results: dict[int, list[str]] = {}
    with ThreadPoolExecutor(max_workers=min(_MAX_CHUNK_WORKERS, len(chunks))) as pool:
        for chunk_result in pool.map(lambda c: _xrefs_chunk(c, xref_type), chunks):
            results.update(chunk_result)
    return results


def chembl_ids_many(cids: "Sequence[int | str]") -> dict[int, str | None]:
    """The ChEMBL ID PubChem has on file for each CID in a batch, if any.

    Batched version of chembl_id(); see xrefs_many() for the round-trip
    savings. A CID with no ChEMBL id on file maps to None, same as
    chembl_id()'s return for a single miss. Unlike xrefs_many() itself
    (whose per-CID presence depends on its chunk-mates; see its
    docstring), every requested CID gets an entry here regardless, since
    this loops over the input CIDs explicitly rather than trusting which
    keys xrefs_many() happened to return.
    """
    all_refs = xrefs_many(cids, xref_type="RegistryID")
    result: dict[int, str | None] = {}
    for cid in cids:
        cid_int = int(cid)
        result[cid_int] = next(
            (ref for ref in all_refs.get(cid_int, []) if ref.startswith("CHEMBL")), None
        )
    return result
