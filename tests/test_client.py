from scigantic_pubchem._client import _parse_throttle_header


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
