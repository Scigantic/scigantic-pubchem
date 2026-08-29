"""Real queries against the live PUG REST API, no mocks, matching the rest
of this package's test suite. AID 1 and AID 3 (DTP/NCI's NCI-H23 and
NCI-H226 growth inhibition screens) are used repeatedly here rather than
picking a fresh AID per test: both are large, stable, decades-old public
assays, and every non-download call here goes through the same JSON cache
every other function in this package does, so reusing them costs one real
network round trip per distinct (path, params), not one per test.

download_assay_results() is the exception: it deliberately bypasses the
cache (see _client.stream_to_file), so each download test below pays its
own real transfer.
"""

import csv
import json
from unittest import mock

import pytest

import scigantic_pubchem as pubchem
from scigantic_pubchem import _client


def test_assay_summary_known_assay():
    summary = pubchem.assay_summary(1)
    assert summary is not None
    assert summary.aid == 1
    assert summary.source_name == "DTP/NCI"
    assert summary.method is not None
    assert summary.description is not None
    # Verified live 2026-08-29: 55,532 total CIDs tested, 3,370 active.
    assert summary.cid_total is not None and summary.cid_total > 0
    assert summary.cid_active is not None and summary.cid_active > 0
    assert summary.cid_active <= summary.cid_total


def test_assay_summary_unknown_returns_none():
    assert pubchem.assay_summary(987654) is None


def test_assay_summary_target_panel():
    # AID 1433 is a multi-target kinase panel; verified live 2026-08-29 it
    # carries a non-empty Target list, unlike AID 1's cell-line phenotypic
    # screen (no molecular target of its own).
    summary = pubchem.assay_summary(1433)
    assert summary is not None
    assert len(summary.target) > 0
    assert all("accession" in t and "name" in t for t in summary.target)


def test_assay_results_known_assay():
    results = pubchem.assay_results(1)
    assert len(results) > 50000  # verified live 2026-08-29: 57,697 rows
    first = results[0]
    assert first.aid == 1
    assert first.sid is not None
    assert first.cid is not None
    assert first.activity_outcome in {"Active", "Inactive", "Inconclusive", "Unspecified", None}


def test_assay_results_unknown_returns_empty():
    assert pubchem.assay_results(987654) == []


def test_assay_results_batches_multiple_aids():
    # PUG REST accepts a comma-separated AID list for `concise` in one
    # request; verified live 2026-08-29 against AIDs 1 and 3 together.
    results = pubchem.assay_results([1, 3])
    aids_seen = {r.aid for r in results}
    assert aids_seen == {1, 3}


def test_assay_cids_and_sids_include_known_ids():
    # Verified live 2026-08-29: SID 66954 / CID 11122 is the first row of
    # AID 3's concise table.
    assert 11122 in pubchem.assay_cids(3)
    assert 66954 in pubchem.assay_sids(3)


def test_assay_cids_unknown_returns_empty():
    assert pubchem.assay_cids(987654) == []
    assert pubchem.assay_sids(987654) == []


def test_aids_for_compound_includes_known_assay():
    aids = pubchem.aids_for_compound(2244)  # aspirin
    assert 1 in aids
    assert 3 in aids


def test_aids_for_compound_unknown_returns_empty():
    assert pubchem.aids_for_compound(999999999) == []


def test_aids_for_target_egfr():
    aids = pubchem.aids_for_target("EGFR")
    assert 1433 in aids  # verified live 2026-08-29


def test_aids_for_target_unknown_returns_empty():
    assert pubchem.aids_for_target("NOTAREALGENE123") == []


def test_compound_assay_results_matches_aids_for_compound():
    results = pubchem.compound_assay_results(2244)
    assert len(results) > 0
    assert all(r.cid == 2244 for r in results)
    aids_from_results = {r.aid for r in results}
    aids_from_lookup = set(pubchem.aids_for_compound(2244))
    assert aids_from_results <= aids_from_lookup


def test_compound_assay_results_unknown_returns_empty():
    assert pubchem.compound_assay_results(999999999) == []


def test_compound_assay_results_many():
    cids = [2244, 3672, 999999999]  # aspirin, ibuprofen, a bogus CID
    batched = pubchem.compound_assay_results_many(cids)
    assert set(batched) == {2244, 3672}  # the bogus CID is absent, not an empty list
    assert all(r.cid == 2244 for r in batched[2244])
    assert all(r.cid == 3672 for r in batched[3672])


def test_download_assay_results_writes_csv(tmp_path):
    dest = tmp_path / "nested" / "aid3.csv"
    result = pubchem.download_assay_results(3, dest)
    assert result == dest
    assert dest.exists()
    with dest.open(newline="") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) > 50000  # verified live 2026-08-29: 54,003 rows
    # Row order isn't guaranteed to be stable across requests (verified live
    # 2026-08-29: CID 11122 was the first row in one request and buried
    # elsewhere in another), so check membership, not position.
    assert all(row["AID"] == "3" for row in rows)
    assert any(row["CID"] == "11122" for row in rows)


def test_download_assay_results_json_format(tmp_path):
    dest = tmp_path / "aid3.json"
    result = pubchem.download_assay_results(3, dest, fmt="json")
    body = json.loads(result.read_text())
    columns = body["Table"]["Columns"]["Column"]
    assert "AID" in columns and "CID" in columns


def test_download_assay_results_batches_multiple_aids(tmp_path):
    dest = tmp_path / "batch.csv"
    pubchem.download_assay_results([1, 3], dest)
    with dest.open(newline="") as f:
        aids_seen = {row["AID"] for row in csv.DictReader(f)}
    assert aids_seen == {"1", "3"}


def test_download_assay_results_unknown_aid_raises_and_leaves_no_file(tmp_path):
    dest = tmp_path / "bad.csv"
    with pytest.raises(_client.CompoundNotFoundError):
        pubchem.download_assay_results(987654, dest)
    assert not dest.exists()


def test_download_assay_results_rejects_bad_format(tmp_path):
    with pytest.raises(ValueError):
        pubchem.download_assay_results(3, tmp_path / "x.csv", fmt="xml")


def test_download_assay_results_leaves_no_partial_file_on_stream_error(tmp_path):
    # The one thing live testing can't exercise: a connection that drops
    # mid-transfer, after a 200 has already started streaming. Verifies
    # the temp-file-then-atomic-rename in _client.stream_to_file (the same
    # pattern cache.py's put() uses) actually cleans up rather than
    # leaving a half-written file where dest is supposed to be.
    import requests

    dest = tmp_path / "interrupted.csv"

    class _BrokenResponse:
        status_code = 200
        headers: dict = {}

        def iter_content(self, chunk_size):
            yield b"AID,SID,CID\n"
            raise requests.exceptions.ChunkedEncodingError("connection closed early")

    session = _client._get_session()
    with mock.patch.object(session, "get", return_value=_BrokenResponse()):
        with mock.patch.object(_client._limiter, "acquire", lambda: None):
            with pytest.raises(requests.exceptions.ChunkedEncodingError):
                pubchem.download_assay_results(3, dest)

    assert not dest.exists()
    assert list(tmp_path.glob("*.part")) == []
