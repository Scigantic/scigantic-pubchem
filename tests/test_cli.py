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


def test_assay_summary_command(capsys):
    exit_code = main(["assay-summary", "1"])
    assert exit_code == 0
    out = capsys.readouterr().out
    assert '"aid": 1' in out
    assert '"source_name": "DTP/NCI"' in out


def test_assay_summary_command_unknown(capsys):
    exit_code = main(["assay-summary", "987654"])
    assert exit_code == 1


def test_aids_for_target_command(capsys):
    exit_code = main(["aids-for-target", "EGFR"])
    assert exit_code == 0
    out = capsys.readouterr().out
    assert "1433" in out.splitlines()


def test_assay_download_command(capsys, tmp_path):
    dest = tmp_path / "aid3.csv"
    exit_code = main(["assay-download", "3", str(dest)])
    assert exit_code == 0
    assert dest.exists()
    assert capsys.readouterr().out.strip() == str(dest)
