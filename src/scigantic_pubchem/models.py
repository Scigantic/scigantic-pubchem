"""Typed records returned by this package -- PubChemPy returns loosely-typed
Compound objects backed by a raw dict; these are plain, typed dataclasses."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Compound:
    cid: int
    title: str | None
    smiles: str | None
    inchi: str | None
    inchi_key: str | None
    iupac_name: str | None
    molecular_formula: str | None
    molecular_weight: float | None
