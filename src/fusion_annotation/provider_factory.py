"""Shared ``DataProvider`` selection logic for the HTTP-facing entry points.

Both the MCP transport (``server/app.py``) and the REST API (``api/app.py``)
need to turn a ``(species, assembly)`` pair plus the
``FUSION_ANNOTATION_PROVIDER`` env var into a concrete provider instance. This
lives in one place so the two servers can't silently drift apart.

Requires the ``server`` extra (``pip install ".[server]"``) for
``GenomeNexusDataProvider``/``RestDataProvider``'s ``requests`` dependency;
importing this module without it will raise ``ImportError`` from those
submodules, same as before the refactor.
"""
from __future__ import annotations

import os

from .gn_provider import GenomeNexusDataProvider
from .rest_provider import RestDataProvider


def make_provider(species: str, assembly: str):
    """Instantiate the configured data provider.

    Reads ``FUSION_ANNOTATION_PROVIDER`` from the environment:
      - unset / "gn" / "genome_nexus" -> GenomeNexusDataProvider (default, fast)
      - "rest" / "ensembl"            -> RestDataProvider (fallback, slower)

    Non-human species always use ``RestDataProvider`` because
    ``GenomeNexusDataProvider`` only covers *Homo sapiens*.
    """
    backend = os.environ.get("FUSION_ANNOTATION_PROVIDER", "gn").strip().lower()
    if backend in ("rest", "ensembl") or (species or "").lower() not in ("homo_sapiens", "human"):
        return RestDataProvider(species=species, assembly=assembly)
    return GenomeNexusDataProvider(assembly=assembly)
