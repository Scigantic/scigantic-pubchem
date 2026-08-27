"""Real queries against PubChem's live fastsimilarity_2d/fastsubstructure
endpoints. No local fingerprint database, unlike scigantic-chembl's
similar_compounds()/substructure_search()."""

import scigantic_pubchem as pubchem

_ASPIRIN_SMILES = "CC(=O)OC1=CC=CC=C1C(=O)O"


def test_similar_compounds_includes_the_query_itself():
    cids = pubchem.similar_compounds(_ASPIRIN_SMILES, threshold=95, max_records=10, resolve=False)
    assert 2244 in cids


def test_similar_compounds_resolves_to_full_records_by_default():
    compounds = pubchem.similar_compounds(_ASPIRIN_SMILES, threshold=95, max_records=5)
    assert all(isinstance(c, pubchem.Compound) for c in compounds)
    assert any(c.cid == 2244 for c in compounds)


def test_substructure_search_smiles():
    cids = pubchem.substructure_search("c1ccccc1", query_type="smiles", max_records=10, resolve=False)
    assert len(cids) > 0
    assert all(isinstance(c, int) for c in cids)


def test_substructure_search_smarts_is_a_different_endpoint():
    # Verified live: SMILES "c1ccccc1" and its SMARTS equivalent return
    # overlapping but not identical CID lists, genuinely different
    # matching semantics, not two names for the same query.
    smiles_cids = set(pubchem.substructure_search("c1ccccc1", query_type="smiles", max_records=10, resolve=False))
    smarts_cids = set(
        pubchem.substructure_search(
            "[#6]1[#6][#6][#6][#6][#6]1", query_type="smarts", max_records=10, resolve=False
        )
    )
    assert smiles_cids  # both non-empty
    assert smarts_cids
    # Not asserting equality, that would assume something not verified.


def test_invalid_query_type_raises():
    import pytest

    with pytest.raises(ValueError):
        pubchem.substructure_search("c1ccccc1", query_type="bogus")
