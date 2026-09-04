import json
from unittest import mock

import pytest

from scigantic_pubchem import _client, cache
from scigantic_pubchem._client import PubChemError, _parse_json, _parse_throttle_header


def _mock_response(status_code, text, headers=None):
    resp = mock.Mock()
    resp.status_code = status_code
    resp.headers = headers or {}
    resp.text = text
    resp.json.side_effect = lambda: json.loads(text)  # matches requests' real strict-by-default behavior
    return resp


def test_parse_json_recovers_from_an_unescaped_control_character():
    # A real, observed PubChem quirk: a depositor free-text field (an assay
    # comment/description) can carry a raw ASCII control character PubChem
    # never escaped, which json.loads/Response.json() reject outright
    # (strict=True, the default) as "Invalid control character".
    text = '{"a": "line1\x01line2"}'
    resp = _mock_response(200, text)
    body = _parse_json(resp, "https://example.invalid/x")
    assert body == {"a": "line1\x01line2"}


def test_parse_json_passes_through_valid_json_unchanged():
    resp = _mock_response(200, '{"a": 1}')
    assert _parse_json(resp, "https://example.invalid/x") == {"a": 1}


def test_parse_json_raises_a_clear_error_for_real_corruption():
    # Truncated mid-object: not a control-character problem, so the
    # strict=False fallback can't rescue it either. Should surface as this
    # package's own PubChemError, not a bare json.JSONDecodeError leaking
    # out of a request() caller who has no reason to expect one.
    resp = _mock_response(200, '{"a": "unterminated')
    with pytest.raises(PubChemError, match="not valid JSON"):
        _parse_json(resp, "https://example.invalid/x")


def test_request_recovers_the_relaxed_parse_end_to_end(tmp_path):
    text = '{"a": "line1' + chr(1) + 'line2"}'
    resp = _mock_response(200, text)
    session = _client._get_session()
    cache.enable_cache(cache_dir=str(tmp_path))
    try:
        with mock.patch.object(session, "get", return_value=resp):
            with mock.patch.object(_client._limiter, "acquire", lambda: None):
                body = _client.request("/some/path")
        assert body == {"a": "line1" + chr(1) + "line2"}
    finally:
        cache.enable_cache()


def test_parses_all_green():
    header = "Request Count status: Green (0%), Request Time status: Green (0%), Service status: Green (27%)"
    assert _parse_throttle_header(header) == "Green"


def test_yellow_beats_green():
    header = "Request Count status: Yellow (60%), Request Time status: Green (0%), Service status: Green (10%)"
    assert _parse_throttle_header(header) == "Yellow"


def test_red_beats_yellow():
    header = "Request Count status: Red (100%), Request Time status: Yellow (60%), Service status: Green (10%)"
    assert _parse_throttle_header(header) == "Red"


def test_missing_header_defaults_green():
    assert _parse_throttle_header(None) == "Green"
