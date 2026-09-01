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

dose_response()/download_dose_response() use AID 1851 (NCGC's 5-isoform
CYP inhibition qHTS panel, 17,143 SIDs), always scoped to an explicit
small sids= subset -- verified live 2026-08-31 that this operation does
not scale linearly (250 SIDs under a second, 2000 SIDs took 94s), so a
full-assay pull has no place in a test suite that runs on every push
across five Python versions. SID 842238 in particular is verified live to
carry rows across all 5 panel targets (a mix of Inactive/Inconclusive
outcomes), so it anchors the single-compound assertions below.
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


# All 6 verified live 2026-08-31 to carry exactly 5 panel rows each (30 total).
_AID_1851_SIDS = [4252546, 4253874, 4254135, 7977146, 11110647, 11111975]


def _data_rows(path):
    """Rows from a download_dose_response() CSV excluding PUG REST's own
    RESULT_TYPE/RESULT_DESCR/RESULT_UNIT/... metadata preamble (kept once,
    from the first chunk, not per data row)."""
    with path.open(newline="") as f:
        return [row for row in csv.DictReader(f) if row["PUBCHEM_RESULT_TAG"].isdigit()]


def test_dose_response_known_compound():
    results = pubchem.dose_response(1851, sids=[842238])
    assert len(results) == 5  # one row per panel target, verified live 2026-08-31
    # "Panel Name" carries the mnemonic (p450-cyp1a2); "Panel Target" is the
    # protein accession (NP_...) -- verified live 2026-08-31, easy to swap.
    names = {r.panel_name for r in results}
    assert names == {"p450-cyp1a2", "p450-cyp2c9", "p450-cyp2c19", "p450-cyp2d6", "p450-cyp3a4"}
    assert all(r.panel_target is not None and r.panel_target.startswith("NP_") for r in results)
    assert all(r.aid == 1851 and r.sid == 842238 for r in results)
    inactive = [r for r in results if r.activity_outcome == "Inactive"]
    assert inactive  # at least one panel target inactive for this SID
    # The whole point: an Inactive row (no fitted potency) still carries a
    # raw Max_Response and per-concentration readout.
    assert any(r.max_response is not None for r in inactive)
    assert any(r.potency_um is None for r in inactive)
    assert all(len(r.dose_response) > 0 for r in results)
    assert all(p.concentration_um > 0 for r in results for p in r.dose_response)


def test_dose_response_by_cid():
    results = pubchem.dose_response(1851, cids=[6602638])
    assert len(results) > 0
    assert all(r.cid == 6602638 for r in results)


def test_dose_response_requires_exactly_one_of_sids_cids():
    with pytest.raises(ValueError):
        pubchem.dose_response(1851)
    with pytest.raises(ValueError):
        pubchem.dose_response(1851, sids=[1], cids=[2])


def test_dose_response_rejects_too_many_ids():
    with pytest.raises(ValueError):
        pubchem.dose_response(1851, sids=list(range(201)))


def test_dose_response_unknown_aid_returns_empty():
    assert pubchem.dose_response(987654, sids=[842238]) == []


def test_download_dose_response_writes_and_resumes(tmp_path):
    dest = tmp_path / "dose_response.csv"
    result = pubchem.download_dose_response(1851, dest, sids=_AID_1851_SIDS, chunk_size=2)
    assert result == dest
    rows = _data_rows(dest)
    assert len(rows) == 5 * len(_AID_1851_SIDS)  # verified live: 5 panel rows per SID
    assert {row["PUBCHEM_SID"] for row in rows} == {str(s) for s in _AID_1851_SIDS}
    # PUG REST's metadata preamble (RESULT_TYPE/RESULT_DESCR/RESULT_UNIT/
    # RESULT_ATTR_CONC_MICROMOL) is repeated in every chunk's own response,
    # but should land in the combined file exactly once (from the first
    # chunk), not once per chunk.
    with dest.open(newline="") as f:
        all_rows = list(csv.DictReader(f))
    assert len(all_rows) - len(rows) == 4
    assert not dest.with_name(dest.name + ".progress.json").exists()  # cleaned up on completion


def test_download_dose_response_resumes_after_interruption(tmp_path):
    dest = tmp_path / "resumed.csv"
    progress_path = dest.with_name(dest.name + ".progress.json")
    real_request_text = _client.request_text

    calls = {"n": 0}

    def _flaky(*args, **kwargs):
        calls["n"] += 1
        if calls["n"] == 2:
            raise _client.PubChemError("simulated failure")
        return real_request_text(*args, **kwargs)

    with mock.patch("scigantic_pubchem.bioassay._client.request_text", side_effect=_flaky):
        with pytest.raises(_client.PubChemError):
            pubchem.download_dose_response(1851, dest, sids=_AID_1851_SIDS, chunk_size=2)

    assert dest.exists()
    assert progress_path.exists()
    partial_state = json.loads(progress_path.read_text())
    assert partial_state["done"] == 1  # first chunk succeeded, second raised

    result = pubchem.download_dose_response(1851, dest, sids=_AID_1851_SIDS, chunk_size=2)
    assert result == dest
    assert not progress_path.exists()
    rows = _data_rows(dest)
    # Exactly the expected row count: chunk 1 (the interrupted one) was
    # truncated away and redone cleanly, not double-applied, and the
    # 4-row metadata preamble was not repeated either.
    assert len(rows) == 5 * len(_AID_1851_SIDS)
    assert {row["PUBCHEM_SID"] for row in rows} == {str(s) for s in _AID_1851_SIDS}
    assert dest.read_text().count("PUBCHEM_RESULT_TAG") == 1  # one header line, not one per chunk


def test_download_dose_response_discards_mismatched_resume_state(tmp_path):
    dest = tmp_path / "mismatched.csv"
    real_request_text = _client.request_text

    calls = {"n": 0}

    def _fail_second(*args, **kwargs):
        calls["n"] += 1
        if calls["n"] == 2:
            raise _client.PubChemError("simulated failure")
        return real_request_text(*args, **kwargs)

    # Leave a genuine partial-progress state behind (chunk 1 of 3 done).
    with mock.patch("scigantic_pubchem.bioassay._client.request_text", side_effect=_fail_second):
        with pytest.raises(_client.PubChemError):
            pubchem.download_dose_response(1851, dest, sids=_AID_1851_SIDS[:3], chunk_size=1)
    progress_path = dest.with_name(dest.name + ".progress.json")
    assert progress_path.exists()

    with pytest.warns(UserWarning, match="stale resume state"):
        # Different sids -> different fingerprint -> discarded, not misapplied.
        result = pubchem.download_dose_response(1851, dest, sids=_AID_1851_SIDS[3:5], chunk_size=1)
    assert not progress_path.exists()
    rows = _data_rows(dest)
    assert {row["PUBCHEM_SID"] for row in rows} == {str(s) for s in _AID_1851_SIDS[3:5]}
    assert result == dest


def test_download_dose_response_unknown_aid_raises():
    with pytest.raises(ValueError):
        pubchem.download_dose_response(987654, "unused.csv")
