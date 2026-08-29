"""Query PubChem live via PUG REST, with caching and rate-limit-aware
retry PubChemPy doesn't have, plus live cross-references into
scigantic-chembl and scigantic-bindingdb. No mirror, no download."""

from importlib.metadata import PackageNotFoundError, version as _version

from .bioassay import (
    aids_for_compound,
    aids_for_target,
    assay_cids,
    assay_results,
    assay_sids,
    assay_summary,
    compound_assay_results,
    compound_assay_results_many,
    download_assay_results,
)
from .bridge import bindingdb_measurements, chembl_context
from .cache import cache_dir, disable_cache, enable_cache, is_cache_enabled
from .cache import clear as clear_cache
from .models import AssayResult, AssaySummary, Compound
from .resolve import resolve, resolve_many
from .similarity import similar_compounds, substructure_search
from .xrefs import chembl_id, chembl_ids_many, xrefs, xrefs_many

try:
    __version__ = _version("scigantic-pubchem")
except PackageNotFoundError:
    # Running from a source checkout with no install (editable or not).
    __version__ = "0.0.0"

__all__ = [
    "resolve",
    "resolve_many",
    "similar_compounds",
    "substructure_search",
    "xrefs",
    "xrefs_many",
    "chembl_id",
    "chembl_ids_many",
    "chembl_context",
    "bindingdb_measurements",
    "assay_summary",
    "assay_results",
    "assay_cids",
    "assay_sids",
    "compound_assay_results",
    "compound_assay_results_many",
    "aids_for_compound",
    "aids_for_target",
    "download_assay_results",
    "Compound",
    "AssaySummary",
    "AssayResult",
    "enable_cache",
    "disable_cache",
    "is_cache_enabled",
    "cache_dir",
    "clear_cache",
]
