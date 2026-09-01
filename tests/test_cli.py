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


def test_pdb_structures_command(capsys):
    exit_code = main(["pdb-structures", "5291"])
    assert exit_code == 0
    out = capsys.readouterr().out
    assert len(out.strip().splitlines()) > 10


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


def test_assay_cids_command(capsys):
    exit_code = main(["assay-cids", "3"])
    assert exit_code == 0
    out = capsys.readouterr().out
    assert "11122" in out.splitlines()


def test_assay_sids_command(capsys):
    exit_code = main(["assay-sids", "3"])
    assert exit_code == 0
    out = capsys.readouterr().out
    assert "66954" in out.splitlines()


def test_assay_download_command(capsys, tmp_path):
    dest = tmp_path / "aid3.csv"
    exit_code = main(["assay-download", "3", str(dest)])
    assert exit_code == 0
    assert dest.exists()
    assert capsys.readouterr().out.strip() == str(dest)


def test_dose_response_command(capsys):
    exit_code = main(["assay-dose-response", "1851", "--sid", "842238"])
    assert exit_code == 0
    out = capsys.readouterr().out
    assert '"aid": 1851' in out
    assert out.strip().count("\n") == 4  # 5 rows, one per panel target


def test_dose_response_download_command(capsys, tmp_path):
    dest = tmp_path / "dose_response.csv"
    exit_code = main(
        ["assay-dose-response-download", "1851", str(dest), "--sid", "842238", "--chunk-size", "2"]
    )
    assert exit_code == 0
    assert dest.exists()
    assert capsys.readouterr().out.strip() == str(dest)


def test_gene_info_command(capsys):
    exit_code = main(["gene-info", "EGFR"])
    assert exit_code == 0
    out = capsys.readouterr().out
    assert '"gene_id": 1956' in out


def test_gene_info_command_unknown(capsys):
    exit_code = main(["gene-info", "NOTAREALGENE123"])
    assert exit_code == 1


def test_protein_info_command(capsys):
    exit_code = main(["protein-info", "P00533"])
    assert exit_code == 0
    out = capsys.readouterr().out
    assert '"accession": "P00533"' in out


def test_tox21_results_command(capsys):
    exit_code = main(["tox21-results", "SR-ATAD5"])
    assert exit_code == 0
    out = capsys.readouterr().out
    assert '"aid": 720516' in out


def test_tox21_matrix_command(capsys):
    exit_code = main(["tox21-matrix", "SR-ATAD5"])
    assert exit_code == 0
    out = capsys.readouterr().out
    assert '"SR-ATAD5"' in out
