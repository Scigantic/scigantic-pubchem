"""Command-line interface: `scigantic-pubchem resolve` / `chembl-id` / `xrefs`."""

from __future__ import annotations

import argparse
import dataclasses
import json
import sys
from typing import Sequence, cast

from .resolve import resolve
from .similarity import similar_compounds, substructure_search
from .xrefs import chembl_id, xrefs


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

    args = parser.parse_args(argv)
    return cast(int, args.func(args))


if __name__ == "__main__":
    sys.exit(main())
