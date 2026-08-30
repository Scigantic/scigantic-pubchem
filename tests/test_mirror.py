import os

import scigantic_pubchem as pubchem


def test_identifiers_for_aspirin():
    rec = pubchem.identifiers(2244)
    assert rec is not None
    assert rec["cid"] == 2244
    assert rec["smiles"]
    assert rec["title"]
    assert rec["molecular_formula"]


def test_identifiers_missing_cid_returns_none():
    assert pubchem.identifiers(999_999_999_999) is None


def test_inchi_keys_multiple_protonation_states():
    df = pubchem.inchi_keys(2244)
    assert len(df) >= 1
    assert "inchikey" in df.columns


def test_synonyms_includes_common_name():
    names = pubchem.synonyms(2244)
    assert any("aspirin" in n.lower() for n in names)


def test_parent_returns_int_or_none():
    p = pubchem.parent(2244)
    assert p is None or isinstance(p, int)


def test_substance_ids_returns_dataframe_with_link_type():
    df = pubchem.substance_ids(2244)
    assert "sid" in df.columns
    assert "link_type" in df.columns
    if len(df) > 0:
        assert set(df["link_type"].unique()).issubset({1, 2})


def test_download_mirror_subset(tmp_path):
    paths = pubchem.download_mirror(str(tmp_path), tables=["smiles"])
    assert "smiles" in paths
    assert os.path.exists(paths["smiles"])
    assert os.path.getsize(paths["smiles"]) > 0


def test_download_mirror_skips_already_present_file(tmp_path):
    dest = tmp_path
    first = pubchem.download_mirror(str(dest), tables=["parent"])
    mtime_before = os.path.getmtime(first["parent"])
    second = pubchem.download_mirror(str(dest), tables=["parent"])
    mtime_after = os.path.getmtime(second["parent"])
    assert mtime_before == mtime_after  # re-fetched only if size differs
