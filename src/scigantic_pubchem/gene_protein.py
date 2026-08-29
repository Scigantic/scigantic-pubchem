"""Gene and protein info, and the bioactivity table keyed by either: a
third direction alongside bioassay.py's per-assay and per-compound views,
this time keyed by the target itself.

PubChem's `gene` and `protein` domains have no PubChemPy equivalent at
all: verified 2026-08-29 by reading pubchempy.py's source directly, it has
no Gene or Protein class, no get_genes()/get_proteins(), and none of
"genesymbol", "geneid", or "ProteinAccession" appear in it anywhere.
AssayResult already carries target_accession/target_gene_id (bioassay.py),
so this module gives that field somewhere to resolve to: a human-readable
name, taxonomy, and description for the gene or protein a bioactivity row
is actually about, plus (via gene_assay_results()/protein_assay_results())
every bioactivity row recorded against that target directly, without
having to already know which AIDs to look in.
"""

from __future__ import annotations

from typing import Any, Literal
from urllib.parse import quote

from . import _client
from ._tables import parse_table, row_to_result
from .models import AssayResult, GeneInfo, ProteinInfo

GeneNamespace = Literal["genesymbol", "geneid"]


def _record_to_gene_info(rec: dict[str, Any]) -> GeneInfo:
    return GeneInfo(
        gene_id=int(rec["GeneID"]),
        symbol=rec.get("Symbol") or None,
        name=rec.get("Name") or None,
        taxonomy_id=rec.get("TaxonomyID"),
        taxonomy=rec.get("Taxonomy") or None,
        description=rec.get("Description") or None,
        synonyms=list(rec.get("Synonym") or []),
    )


def _record_to_protein_info(rec: dict[str, Any]) -> ProteinInfo:
    return ProteinInfo(
        accession=str(rec["ProteinAccession"]),
        name=rec.get("Name") or None,
        taxonomy_id=rec.get("TaxonomyID"),
        taxonomy=rec.get("Taxonomy") or None,
        synonyms=list(rec.get("Synonym") or []),
    )


def gene_info(identifier: str | int, namespace: GeneNamespace = "genesymbol") -> GeneInfo | None:
    """The overview PubChem has on file for one gene: symbol, full name,
    taxonomy, description, and known synonyms. namespace is "genesymbol"
    (default, e.g. "EGFR") or "geneid" (NCBI's numeric Entrez Gene ID).
    None if PubChem has no gene record for it.
    """
    ident = quote(str(identifier), safe="")
    try:
        body = _client.request(f"/gene/{namespace}/{ident}/summary/JSON")
    except _client.CompoundNotFoundError:
        return None
    records = body.get("GeneSummaries", {}).get("GeneSummary", [])
    if not records:
        return None
    return _record_to_gene_info(records[0])


def protein_info(accession: str) -> ProteinInfo | None:
    """The overview PubChem has on file for one protein: name, taxonomy,
    and known synonyms, keyed by its accession (e.g. UniProt "P00533").
    None if PubChem has no protein record for it.
    """
    ident = quote(accession, safe="")
    try:
        body = _client.request(f"/protein/accession/{ident}/summary/JSON")
    except _client.CompoundNotFoundError:
        return None
    records = body.get("ProteinSummaries", {}).get("ProteinSummary", [])
    if not records:
        return None
    return _record_to_protein_info(records[0])


def gene_assay_results(identifier: str | int, namespace: GeneNamespace = "genesymbol") -> list[AssayResult]:
    """Every bioactivity row recorded across every assay run against one
    gene target, read directly rather than requiring aids_for_target()
    plus a per-AID fetch first. Uses the gene domain's own `concise`
    operation, the same row shape assay_results()/compound_assay_results()
    return (see AssayResult), just selected by gene rather than by assay
    or compound.
    """
    ident = quote(str(identifier), safe="")
    try:
        body = _client.request(f"/gene/{namespace}/{ident}/concise/JSON")
    except _client.CompoundNotFoundError:
        return []
    return [row_to_result(row) for row in parse_table(body)]


def protein_assay_results(accession: str) -> list[AssayResult]:
    """Every bioactivity row recorded across every assay run against one
    protein target, keyed by accession. Same row shape as
    gene_assay_results(), from the protein domain's `concise` operation.
    """
    ident = quote(accession, safe="")
    try:
        body = _client.request(f"/protein/accession/{ident}/concise/JSON")
    except _client.CompoundNotFoundError:
        return []
    return [row_to_result(row) for row in parse_table(body)]
