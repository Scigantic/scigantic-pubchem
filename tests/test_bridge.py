import threading

import scigantic_pubchem as pubchem
import scigantic_pubchem.bridge as _bridge_module


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


def test_base_connection_is_reused_not_rebuilt_per_call():
    _bridge_module._base_con = None  # force a fresh base connection
    pubchem.chembl_context(2244)
    first = _bridge_module._base_con
    assert first is not None
    pubchem.bindingdb_measurements(2244)
    assert _bridge_module._base_con is first  # same base connection, not rebuilt


def test_concurrent_calls_return_correct_distinct_results():
    # Guards against a real bug, not a hypothetical one: a single shared
    # DuckDB Connection is not safe for concurrent execute()/fetchone()
    # calls from multiple threads. Verified directly against this exact
    # scenario (two threads racing execute()/fetchone() on one connection)
    # before this fix: some calls silently returned None instead of the
    # real row, not an exception, a corrupted result. _get_connection()
    # hands out a fresh cursor() per call specifically to avoid this.
    cids = [2244, 2519, 1983]
    sequential = {cid: pubchem.chembl_context(cid) for cid in cids}

    results = {}
    lock = threading.Lock()

    def worker(cid):
        ctx = pubchem.chembl_context(cid)
        with lock:
            results[cid] = ctx

    threads = [threading.Thread(target=worker, args=(cid,)) for cid in cids]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    assert results == sequential
