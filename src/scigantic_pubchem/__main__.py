"""Enables `python -m scigantic_pubchem`, same commands as the `scigantic-pubchem` console script."""

import sys

from .cli import main

if __name__ == "__main__":
    sys.exit(main())
