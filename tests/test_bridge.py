import scigantic_pubchem as pubchem


def test_chembl_context_for_aspirin():
    ctx = pubchem.chembl_context(2244)
    assert ctx is not None
    assert ctx["chembl_id"] == "CHEMBL25"
    assert ctx["pref_name"] == "ASPIRIN"
    assert ctx["max_phase"] == 4


def test_chembl_context_no_xref_returns_none():
    ctx = pubchem.chembl_context(999999999)
    assert ctx is None


def test_bindingdb_measurements_returns_dataframe():
    df = pubchem.bindingdb_measurements(2244)
    assert "reactant_set_id" in df.columns
