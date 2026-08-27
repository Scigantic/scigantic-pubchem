from scigantic_pubchem.cli import main


def test_resolve_command(capsys):
    exit_code = main(["resolve", "aspirin"])
    assert exit_code == 0
    out = capsys.readouterr().out
    assert '"cid": 2244' in out


def test_chembl_id_command(capsys):
    exit_code = main(["chembl-id", "2244"])
    assert exit_code == 0
    out = capsys.readouterr().out
    assert "CHEMBL25" in out


def test_resolve_command_unknown(capsys):
    exit_code = main(["resolve", "this is not a real compound xyzzy123"])
    assert exit_code == 1
