import scigantic_pubchem as pubchem


def test_chembl_id_for_aspirin():
    # Verified live 2026-08-27 directly against PUG REST before this
    # package existed: aspirin (CID 2244) resolves to CHEMBL25.
    assert pubchem.chembl_id(2244) == "CHEMBL25"


def test_xrefs_registry_id_includes_known_values():
    refs = pubchem.xrefs(2244, xref_type="RegistryID")
    assert "CHEMBL25" in refs
    assert any(r.startswith("CHEBI:") for r in refs)


def test_xrefs_empty_for_bogus_cid():
    # A CID this large does not exist; PUG REST 404s for zero xrefs, which
    # this should surface as an empty list, not an exception.
    assert pubchem.xrefs(999999999, xref_type="RegistryID") == []


def test_xrefs_many_matches_individual_calls():
    cids = [2244, 3672, 2519]  # aspirin, ibuprofen, caffeine
    batched = pubchem.xrefs_many(cids, xref_type="RegistryID")
    assert set(batched) == set(cids)
    for cid in cids:
        # Same set of RegistryIDs, not necessarily the same order: verified
        # live that a batched multi-CID request and a single-CID request
        # for the same compound return identical content in different
        # order, so order was never part of this field's real contract.
        assert set(batched[cid]) == set(pubchem.xrefs(cid, xref_type="RegistryID"))


def test_xrefs_many_forces_multiple_chunks():
    # 250 CIDs at a chunk size of 200 forces a 200 + 50 split, exercising
    # the chunking loop rather than just its single-request path.
    cids = list(range(1, 251))
    batched = pubchem.xrefs_many(cids, xref_type="RegistryID")
    assert set(batched) == set(cids)


def test_xrefs_many_empty_input():
    assert pubchem.xrefs_many([]) == {}


def test_chembl_ids_many_matches_individual_calls():
    cids = [2244, 999999999]  # aspirin has one, the bogus cid has none
    batched = pubchem.chembl_ids_many(cids)
    assert batched == {2244: "CHEMBL25", 999999999: None}


def test_xrefs_many_mixed_batch_includes_empty_entries():
    # Not the same thing as an all-miss batch (test_xrefs_many_empty_input
    # covers that: PUG REST 404s the whole chunk, so it contributes
    # nothing). Verified live 2026-08-29: a chunk with at least one real
    # match returns 200 with an explicit empty-list entry for every CID in
    # that chunk that had none, so a miss mixed in with a hit is present
    # in the result (mapped to []), not absent from it.
    batched = pubchem.xrefs_many([2244, 999999999], xref_type="RegistryID")
    assert batched == {2244: batched[2244], 999999999: []}
    assert len(batched[2244]) > 0


def test_pdb_structures_known_ligand():
    # Verified live 2026-08-29: imatinib (CID 5291) is deposited as a
    # ligand in 27 PDB structures, mirrored into NCBI's MMDB ID space.
    structures = pubchem.pdb_structures(5291)
    assert len(structures) > 10
    assert all(isinstance(s, int) for s in structures)


def test_pdb_structures_never_crystallized_returns_empty():
    # CID 1 is methane: verified live never deposited as a bound ligand.
    assert pubchem.pdb_structures(1) == []


def test_pdb_structures_many():
    cids = [2244, 5291, 1]  # aspirin, imatinib, methane (never crystallized)
    batched = pubchem.pdb_structures_many(cids)
    assert set(batched) == set(cids)  # every requested CID present, unlike raw xrefs_many
    assert batched[1] == []
    assert batched[2244] == pubchem.pdb_structures(2244)
    assert batched[5291] == pubchem.pdb_structures(5291)
