"""Real queries against the live PUG REST API, no mocks. The whole point
of this package is that it answers real queries with no local setup."""

import scigantic_pubchem as pubchem


def test_resolve_by_name():
    aspirin = pubchem.resolve("aspirin")
    assert aspirin is not None
    assert aspirin.cid == 2244
    assert aspirin.inchi_key == "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
    assert aspirin.molecular_formula == "C9H8O4"
    assert aspirin.title == "Aspirin"


def test_resolve_by_cid():
    compound = pubchem.resolve(2244, namespace="cid")
    assert compound is not None
    assert compound.cid == 2244


def test_resolve_by_inchikey():
    compound = pubchem.resolve("BSYNRYMUTXBXSQ-UHFFFAOYSA-N", namespace="inchikey")
    assert compound is not None
    assert compound.cid == 2244


def test_resolve_unknown_returns_none():
    assert pubchem.resolve("this is not a real compound name xyzzy123") is None


def test_resolve_many():
    compounds = pubchem.resolve_many([2244, 2519, 1983])
    cids = {c.cid for c in compounds}
    assert cids == {2244, 2519, 1983}


def test_resolve_many_forces_multiple_chunks():
    # Bulk resolution, not just the single-request case: 250 CIDs at a
    # chunk size of 200 forces a 200 + 50 split, exercising the actual
    # chunking loop rather than just its trivial one-request path.
    # Verified live 2026-08-27 that CIDs 1-250 all resolve, so a short
    # result is a real signal, not upstream flakiness.
    cids = list(range(1, 251))
    compounds = pubchem.resolve_many(cids)
    resolved_cids = {c.cid for c in compounds}
    assert resolved_cids == set(cids)
    assert len(compounds) == 250
