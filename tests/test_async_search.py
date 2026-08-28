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

import pytest

from scigantic_pubchem import _client


@pytest.fixture(autouse=True)
def _no_rate_limit_wait(monkeypatch):
    """Every test here mocks time.sleep to skip real waiting -- but
    _client.request() also calls the shared _RateLimiter.acquire() before
    each attempt, which itself calls time.sleep() while waiting for a
    token to refill. With time.sleep mocked to a no-op and time.monotonic
    left real, that refill wait can only be satisfied by real (tiny) wall
    time passing between spins, so a test making more calls than the
    limiter's capacity spins through thousands of fast, pointless
    iterations instead of the handful of sleeps it's actually checking.
    Rate limiting is its own concern, covered by test_rate_limit.py; it's
    bypassed here rather than making every test in this file also fake a
    clock granular enough for both concerns at once.
    """
    monkeypatch.setattr(_client._limiter, "acquire", lambda: None)


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


def test_poll_wait_doubles_then_caps():
    # 1 initial call + 6 polls (the first 5 still waiting, the 6th done)
    # = 7 GETs, so 6 sleeps: 0.5, 1.0, 2.0, 4.0, then capped at 5.0 rather
    # than growing to 8.0 and 10.0 for the last two.
    waiting = _mock_response(200, {"Waiting": {"ListKey": "abc123"}})
    done = _mock_response(200, {"IdentifierList": {"CID": [2244]}})
    responses = [waiting] * 6 + [done]

    session = _client._get_session()
    with mock.patch.object(session, "get", side_effect=responses):
        with mock.patch("time.sleep") as sleep:
            body = _client.request_search("/compound/fastsimilarity_2d/smiles/cids/JSON", params={"smiles": "C"})

    assert body == {"IdentifierList": {"CID": [2244]}}
    waited = [call.args[0] for call in sleep.call_args_list]
    assert waited == [0.5, 1.0, 2.0, 4.0, 5.0, 5.0]
