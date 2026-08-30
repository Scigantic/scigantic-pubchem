"""Command-line interface: `scigantic-pubchem resolve` / `chembl-id` / `xrefs`."""

from __future__ import annotations

import argparse
import dataclasses
import json
import sys
from typing import Sequence, cast

from .bioassay import (
    aids_for_compound,
    aids_for_target,
    assay_results,
    assay_summary,
    compound_assay_results,
    download_assay_results,
)
from .gene_protein import (
    gene_assay_results,
    gene_info,
    protein_assay_results,
    protein_info,
)
from .mirror import download as mirror_download
from .mirror import identifiers as mirror_identifiers
from .resolve import resolve
from .similarity import similar_compounds, substructure_search
from .tox21 import tox21_matrix, tox21_results
from .xrefs import chembl_id, pdb_structures, xrefs


def _cmd_resolve(args: argparse.Namespace) -> int:
    compound = resolve(args.identifier, namespace=args.namespace)
    if compound is None:
        print(f"no match for {args.identifier!r}", file=sys.stderr)
        return 1
    print(json.dumps(dataclasses.asdict(compound), indent=2))
    return 0


def _cmd_chembl_id(args: argparse.Namespace) -> int:
    result = chembl_id(args.identifier, namespace=args.namespace)
    if result is None:
        print(f"no ChEMBL cross-reference for {args.identifier!r}", file=sys.stderr)
        return 1
    print(result)
    return 0


def _cmd_mirror_identifiers(args: argparse.Namespace) -> int:
    result = mirror_identifiers(args.cid)
    if result is None:
        print(f"CID {args.cid} not in the mirror", file=sys.stderr)
        return 1
    print(json.dumps(result, indent=2))
    return 0


def _cmd_mirror_download(args: argparse.Namespace) -> int:
    paths = mirror_download(args.dest, tables=args.table or None)
    for table, path in paths.items():
        print(f"{table}: {path}")
    return 0


def _cmd_pdb_structures(args: argparse.Namespace) -> int:
    for mmdb_id in pdb_structures(args.identifier, namespace=args.namespace):
        print(mmdb_id)
    return 0


def _cmd_xrefs(args: argparse.Namespace) -> int:
    for ref in xrefs(args.identifier, xref_type=args.type, namespace=args.namespace):
        print(ref)
    return 0


def _cmd_similar(args: argparse.Namespace) -> int:
    cids = cast(
        "list[int]",
        similar_compounds(args.smiles, threshold=args.threshold, max_records=args.max_records, resolve=False),
    )
    for cid in cids:
        print(cid)
    return 0


def _cmd_substructure(args: argparse.Namespace) -> int:
    cids = cast(
        "list[int]",
        substructure_search(args.query, query_type=args.query_type, max_records=args.max_records, resolve=False),
    )
    for cid in cids:
        print(cid)
    return 0


def _cmd_assay_summary(args: argparse.Namespace) -> int:
    summary = assay_summary(args.aid)
    if summary is None:
        print(f"no assay found for AID {args.aid}", file=sys.stderr)
        return 1
    print(json.dumps(dataclasses.asdict(summary), indent=2))
    return 0


def _cmd_assay_results(args: argparse.Namespace) -> int:
    for result in assay_results(args.aid):
        print(json.dumps(dataclasses.asdict(result)))
    return 0


def _cmd_compound_assay_results(args: argparse.Namespace) -> int:
    for result in compound_assay_results(args.identifier, namespace=args.namespace):
        print(json.dumps(dataclasses.asdict(result)))
    return 0


def _cmd_aids_for_compound(args: argparse.Namespace) -> int:
    for aid in aids_for_compound(args.identifier, namespace=args.namespace):
        print(aid)
    return 0


def _cmd_aids_for_target(args: argparse.Namespace) -> int:
    for aid in aids_for_target(args.gene_symbol):
        print(aid)
    return 0


def _cmd_assay_download(args: argparse.Namespace) -> int:
    dest = download_assay_results(args.aid, args.dest, fmt=args.format)
    print(dest)
    return 0


def _cmd_gene_info(args: argparse.Namespace) -> int:
    info = gene_info(args.identifier, namespace=args.namespace)
    if info is None:
        print(f"no gene found for {args.identifier!r}", file=sys.stderr)
        return 1
    print(json.dumps(dataclasses.asdict(info), indent=2))
    return 0


def _cmd_protein_info(args: argparse.Namespace) -> int:
    info = protein_info(args.accession)
    if info is None:
        print(f"no protein found for {args.accession!r}", file=sys.stderr)
        return 1
    print(json.dumps(dataclasses.asdict(info), indent=2))
    return 0


def _cmd_gene_assay_results(args: argparse.Namespace) -> int:
    for result in gene_assay_results(args.identifier, namespace=args.namespace):
        print(json.dumps(dataclasses.asdict(result)))
    return 0


def _cmd_protein_assay_results(args: argparse.Namespace) -> int:
    for result in protein_assay_results(args.accession):
        print(json.dumps(dataclasses.asdict(result)))
    return 0


def _cmd_tox21_results(args: argparse.Namespace) -> int:
    for result in tox21_results(args.endpoint or None):
        print(json.dumps(dataclasses.asdict(result)))
    return 0


def _cmd_tox21_matrix(args: argparse.Namespace) -> int:
    print(json.dumps(tox21_matrix(args.endpoint or None), indent=2))
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="scigantic-pubchem")
    subparsers = parser.add_subparsers(dest="command", required=True)

    resolve_parser = subparsers.add_parser("resolve", help="resolve a name/SMILES/InChIKey/CID to a compound record")
    resolve_parser.add_argument("identifier")
    resolve_parser.add_argument("--namespace", default="name", choices=["name", "cid", "smiles", "inchikey", "inchi", "formula"])
    resolve_parser.set_defaults(func=_cmd_resolve)

    chembl_parser = subparsers.add_parser("chembl-id", help="live ChEMBL cross-reference for a compound")
    chembl_parser.add_argument("identifier")
    chembl_parser.add_argument("--namespace", default="cid", choices=["name", "cid", "smiles", "inchikey", "inchi", "formula"])
    chembl_parser.set_defaults(func=_cmd_chembl_id)

    pdb_parser = subparsers.add_parser("pdb-structures", help="MMDB IDs of structures with this compound bound as a ligand")
    pdb_parser.add_argument("identifier")
    pdb_parser.add_argument("--namespace", default="cid", choices=["name", "cid", "smiles", "inchikey", "inchi", "formula"])
    pdb_parser.set_defaults(func=_cmd_pdb_structures)

    mirror_id_parser = subparsers.add_parser(
        "mirror-identifiers", help="CID-keyed identifier/name/mass lookup from the local parquet mirror (no live call)"
    )
    mirror_id_parser.add_argument("cid", type=int)
    mirror_id_parser.set_defaults(func=_cmd_mirror_identifiers)

    mirror_download_parser = subparsers.add_parser(
        "mirror-download", help="download the mirror's parquet files for local/offline DuckDB use"
    )
    mirror_download_parser.add_argument("dest", help="destination directory")
    mirror_download_parser.add_argument(
        "--table", action="append", help="table to download (repeatable; default: all eight)"
    )
    mirror_download_parser.set_defaults(func=_cmd_mirror_download)

    xrefs_parser = subparsers.add_parser("xrefs", help="raw cross-references PubChem has on file")
    xrefs_parser.add_argument("identifier")
    xrefs_parser.add_argument("--type", default="RegistryID")
    xrefs_parser.add_argument("--namespace", default="cid", choices=["name", "cid", "smiles", "inchikey", "inchi", "formula"])
    xrefs_parser.set_defaults(func=_cmd_xrefs)

    similar_parser = subparsers.add_parser("similar", help="2D similarity search over PubChem's full corpus")
    similar_parser.add_argument("smiles")
    similar_parser.add_argument("--threshold", type=int, default=90)
    similar_parser.add_argument("--max-records", type=int, default=100, dest="max_records")
    similar_parser.set_defaults(func=_cmd_similar)

    substructure_parser = subparsers.add_parser("substructure", help="substructure search over PubChem's full corpus")
    substructure_parser.add_argument("query")
    substructure_parser.add_argument("--query-type", default="smiles", choices=["smiles", "smarts"], dest="query_type")
    substructure_parser.add_argument("--max-records", type=int, default=100, dest="max_records")
    substructure_parser.set_defaults(func=_cmd_substructure)

    assay_summary_parser = subparsers.add_parser("assay-summary", help="overview and outcome counts for one assay")
    assay_summary_parser.add_argument("aid", type=int)
    assay_summary_parser.set_defaults(func=_cmd_assay_summary)

    assay_results_parser = subparsers.add_parser("assay-results", help="bioactivity rows for one or more assays")
    assay_results_parser.add_argument("aid", nargs="+", help="one or more AIDs")
    assay_results_parser.set_defaults(func=_cmd_assay_results)

    compound_results_parser = subparsers.add_parser(
        "compound-assay-results", help="every bioactivity row recorded for one compound"
    )
    compound_results_parser.add_argument("identifier")
    compound_results_parser.add_argument("--namespace", default="cid", choices=["name", "cid", "smiles", "inchikey", "inchi", "formula"])
    compound_results_parser.set_defaults(func=_cmd_compound_assay_results)

    aids_compound_parser = subparsers.add_parser("aids-for-compound", help="AIDs of every assay that tested a compound")
    aids_compound_parser.add_argument("identifier")
    aids_compound_parser.add_argument("--namespace", default="cid", choices=["name", "cid", "smiles", "inchikey", "inchi", "formula"])
    aids_compound_parser.set_defaults(func=_cmd_aids_for_compound)

    aids_target_parser = subparsers.add_parser("aids-for-target", help="AIDs of every assay run against a gene target")
    aids_target_parser.add_argument("gene_symbol")
    aids_target_parser.set_defaults(func=_cmd_aids_for_target)

    assay_download_parser = subparsers.add_parser("assay-download", help="stream an assay's bioactivity table to a file")
    assay_download_parser.add_argument("aid", nargs="+", help="one or more AIDs")
    assay_download_parser.add_argument("dest")
    assay_download_parser.add_argument("--format", default="csv", choices=["csv", "json"])
    assay_download_parser.set_defaults(func=_cmd_assay_download)

    gene_info_parser = subparsers.add_parser("gene-info", help="overview for one gene")
    gene_info_parser.add_argument("identifier")
    gene_info_parser.add_argument("--namespace", default="genesymbol", choices=["genesymbol", "geneid"])
    gene_info_parser.set_defaults(func=_cmd_gene_info)

    protein_info_parser = subparsers.add_parser("protein-info", help="overview for one protein")
    protein_info_parser.add_argument("accession")
    protein_info_parser.set_defaults(func=_cmd_protein_info)

    gene_results_parser = subparsers.add_parser(
        "gene-assay-results", help="every bioactivity row recorded against a gene target"
    )
    gene_results_parser.add_argument("identifier")
    gene_results_parser.add_argument("--namespace", default="genesymbol", choices=["genesymbol", "geneid"])
    gene_results_parser.set_defaults(func=_cmd_gene_assay_results)

    protein_results_parser = subparsers.add_parser(
        "protein-assay-results", help="every bioactivity row recorded against a protein target"
    )
    protein_results_parser.add_argument("accession")
    protein_results_parser.set_defaults(func=_cmd_protein_assay_results)

    tox21_results_parser = subparsers.add_parser(
        "tox21-results", help="raw bioactivity rows for the Tox21 Data Challenge panel"
    )
    tox21_results_parser.add_argument("endpoint", nargs="*", help="Tox21 endpoint(s), e.g. NR-AhR SR-p53 (default: all 12)")
    tox21_results_parser.set_defaults(func=_cmd_tox21_results)

    tox21_matrix_parser = subparsers.add_parser(
        "tox21-matrix", help="wide CID x endpoint label matrix for the Tox21 Data Challenge panel"
    )
    tox21_matrix_parser.add_argument("endpoint", nargs="*", help="Tox21 endpoint(s) (default: all 12)")
    tox21_matrix_parser.set_defaults(func=_cmd_tox21_matrix)

    args = parser.parse_args(argv)
    return cast(int, args.func(args))


if __name__ == "__main__":
    sys.exit(main())
