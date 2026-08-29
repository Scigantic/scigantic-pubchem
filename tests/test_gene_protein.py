"""Real queries against the live PUG REST API, no mocks, matching the rest
of this package's test suite. EGFR (Entrez Gene ID 1956, UniProt accession
P00533) is used throughout: a well-known, stable, heavily-annotated human
gene/protein pair, the same role AID 1/AID 3 and CID 2244 play in
test_bioassay.py.
"""

import scigantic_pubchem as pubchem


def test_gene_info_by_symbol():
    gene = pubchem.gene_info("EGFR")
    assert gene is not None
    assert gene.gene_id == 1956  # verified live 2026-08-29
    assert gene.symbol == "EGFR"
    assert gene.taxonomy_id == 9606
    assert gene.taxonomy is not None and "Homo sapiens" in gene.taxonomy
    assert gene.description is not None
    assert "HER1" in gene.synonyms


def test_gene_info_by_geneid():
    gene = pubchem.gene_info(1956, namespace="geneid")
    assert gene is not None
    assert gene.symbol == "EGFR"


def test_gene_info_unknown_returns_none():
    assert pubchem.gene_info("NOTAREALGENE123") is None


def test_protein_info_by_accession():
    protein = pubchem.protein_info("P00533")
    assert protein is not None
    assert protein.accession == "P00533"
    assert protein.taxonomy_id == 9606
    assert protein.name is not None and "growth factor" in protein.name.lower()


def test_protein_info_unknown_returns_none():
    assert pubchem.protein_info("NOTAREALACCESSION") is None


def test_gene_assay_results_by_symbol():
    results = pubchem.gene_assay_results("EGFR")
    assert len(results) > 1000  # verified live 2026-08-29: 2,930 rows
    first = results[0]
    assert first.target_gene_id is None  # gene-domain concise has no such column: the gene IS the query
    assert first.sid is not None
    assert first.cid is not None


def test_gene_assay_results_by_geneid_matches_symbol():
    by_symbol = {r.aid for r in pubchem.gene_assay_results("EGFR")}
    by_geneid = {r.aid for r in pubchem.gene_assay_results(1956, namespace="geneid")}
    assert by_symbol == by_geneid


def test_gene_assay_results_unknown_returns_empty():
    assert pubchem.gene_assay_results("NOTAREALGENE123") == []


def test_protein_assay_results():
    results = pubchem.protein_assay_results("P00533")
    assert len(results) > 10000  # verified live 2026-08-29: 53,532 rows
    first = results[0]
    assert first.target_accession is None  # protein-domain concise has no such column
    assert first.sid is not None
    assert first.cid is not None


def test_protein_assay_results_activity_qualifier_present():
    # The one column concise/assaysummary (assay- and compound-keyed) never
    # carry: verified live 2026-08-29 that PubChem's gene/protein-domain
    # concise tables add an "Activity Qualifier" column (=, <, >) neither
    # of the other two sources have.
    results = pubchem.protein_assay_results("P00533")
    qualifiers = {r.activity_qualifier for r in results}
    assert qualifiers & {"=", "<", ">"}


def test_protein_assay_results_unknown_returns_empty():
    assert pubchem.protein_assay_results("NOTAREALACCESSION") == []
