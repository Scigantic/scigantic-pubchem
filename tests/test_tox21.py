"""Real queries against the live PUG REST API, no mocks, matching the rest
of this package's test suite. Every test below that hits the full 12-AID
panel calls tox21_results()/tox21_matrix() with the same endpoint order
(TOX21_ENDPOINTS' own), so the joined `concise` path is byte-identical
across tests and only the first pays a real network round trip; the rest
are cache hits (mirrors test_bioassay.py's AID-1/AID-3 reuse discipline).
"""

import pytest

import scigantic_pubchem as pubchem
from scigantic_pubchem.tox21 import TOX21_ENDPOINTS, _consolidate

_ALL_ENDPOINTS = {
    "NR-AR",
    "NR-AR-LBD",
    "NR-AhR",
    "NR-Aromatase",
    "NR-ER",
    "NR-ER-LBD",
    "NR-PPAR-gamma",
    "SR-ARE",
    "SR-ATAD5",
    "SR-HSE",
    "SR-MMP",
    "SR-p53",
}


def test_tox21_endpoints_has_all_twelve():
    assert set(TOX21_ENDPOINTS) == _ALL_ENDPOINTS


def test_tox21_results_single_endpoint():
    results = pubchem.tox21_results(["SR-ATAD5"])
    assert len(results) > 5000  # verified live 2026-08-29: 8,099 CIDs tested, some retested
    assert all(r.aid == TOX21_ENDPOINTS["SR-ATAD5"] for r in results)


def test_tox21_results_unknown_endpoint_raises():
    with pytest.raises(ValueError):
        pubchem.tox21_results(["NOT-A-REAL-ENDPOINT"])


def test_tox21_matrix_unknown_endpoint_raises():
    with pytest.raises(ValueError):
        pubchem.tox21_matrix(["NOT-A-REAL-ENDPOINT"])


def test_tox21_results_full_panel_covers_all_aids():
    # Verified live 2026-08-29: the full panel is ~123,000 rows across the
    # 12 AIDs in one batched `concise` request.
    results = pubchem.tox21_results()
    assert {r.aid for r in results} == set(TOX21_ENDPOINTS.values())
    assert len(results) > 100000


def test_tox21_matrix_full_panel_shape():
    matrix = pubchem.tox21_matrix()
    assert len(matrix) > 5000
    row = next(iter(matrix.values()))
    assert set(row) == _ALL_ENDPOINTS
    assert all(v in (0, 1, None) for v in row.values())


def test_tox21_matrix_subset_endpoints_only_has_those_columns():
    matrix = pubchem.tox21_matrix(["SR-ATAD5", "SR-p53"])
    assert len(matrix) > 0
    assert all(set(row) == {"SR-ATAD5", "SR-p53"} for row in matrix.values())


def test_consolidate_active_wins_over_inactive():
    assert _consolidate(["Inactive", "Active"]) == 1


def test_consolidate_all_inactive_is_zero():
    assert _consolidate(["Inactive", "Inactive"]) == 0


def test_consolidate_inconclusive_only_is_none():
    assert _consolidate(["Inconclusive"]) is None


def test_consolidate_empty_is_none():
    assert _consolidate([]) is None
