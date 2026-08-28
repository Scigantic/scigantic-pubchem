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
