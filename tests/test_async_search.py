"""The one thing live testing could not verify: PUG REST's async job
protocol for a slow search (respond with {"Waiting": {"ListKey": ...}} and
require polling /compound/listkey/{key}/cids/JSON until it resolves).

Every real query tried against the live API during development, down to
a single-carbon substructure search (about as broad as a query gets),
resolved synchronously in under two seconds, so this path could not be
observed live. It is real: PubChemPy's own source implements the identical
protocol, and PubChem's PUG REST documentation describes it for expensive
searches. Scripted here with mocked responses so the polling logic itself
is verified deterministically rather than left checked only against
documentation and a code read of someone else's implementation.
"""

from unittest import mock

from scigantic_pubchem import _client


def _mock_response(status_code, json_body, headers=None):
    resp = mock.Mock()
    resp.status_code = status_code
    resp.headers = headers or {}
    resp.json.return_value = json_body
    resp.text = str(json_body)
    return resp


def test_polls_until_waiting_clears():
    waiting = _mock_response(200, {"Waiting": {"ListKey": "abc123"}})
    still_waiting = _mock_response(200, {"Waiting": {"ListKey": "abc123"}})
    done = _mock_response(200, {"IdentifierList": {"CID": [2244, 5161]}})

    session = _client._get_session()
    with mock.patch.object(session, "get", side_effect=[waiting, still_waiting, done]) as get:
        with mock.patch("time.sleep"):  # don't actually wait 2s x2 in a unit test
            body = _client.request_search("/compound/fastsimilarity_2d/smiles/cids/JSON", params={"smiles": "C"})

    assert body == {"IdentifierList": {"CID": [2244, 5161]}}
    assert get.call_count == 3
    # First call is the real search; the next two poll the listkey endpoint.
    first_url = get.call_args_list[0].args[0]
    poll_urls = [c.args[0] for c in get.call_args_list[1:]]
    assert "fastsimilarity_2d" in first_url
    assert all("listkey/abc123" in u for u in poll_urls)


def test_immediate_result_does_not_poll():
    done = _mock_response(200, {"IdentifierList": {"CID": [2244]}})
    session = _client._get_session()
    with mock.patch.object(session, "get", side_effect=[done]) as get:
        body = _client.request_search("/compound/fastsimilarity_2d/smiles/cids/JSON", params={"smiles": "C"})
    assert body == {"IdentifierList": {"CID": [2244]}}
    assert get.call_count == 1
