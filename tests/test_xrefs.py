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
