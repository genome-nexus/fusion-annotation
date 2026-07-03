"""Genome-build (GRCh38/GRCh37) selection and echo.

A genomic breakpoint only means something against the build it was called on, so
the assembly must be selectable and must never silently fall back to the wrong
build. These tests are offline: they check build normalization, which Ensembl
host each build maps to, the MCP provider's GRCh37 guard, and that the chosen
build is echoed back under resolved.genome_build.
"""
import json
import os

import pytest

from fusion_annotation import (
    Transcript, build_exon_cds_map, build_exon_genomic_map, annotate_fusion,
)
from fusion_annotation.providers import StaticProvider, MCPDataProvider

FIXTURE = os.path.join(os.path.dirname(__file__), "fixtures", "eml4_alk_fixture.json")


def _static_provider():
    fx = json.load(open(FIXTURE))
    txs = {}
    for key, t in fx["transcripts"].items():
        txs[key] = Transcript(
            gene_symbol=t["gene_symbol"], gene_id=t["gene_id"],
            transcript_id=t["transcript_id"], strand=t["strand"],
            cds=t["cds"], protein=t["protein"], uniprot=t["uniprot"],
            exon_cds=build_exon_cds_map(t["strand"], t["exons"], t["cds_g_start"], t["cds_g_end"]),
            exon_genomic=build_exon_genomic_map(t["strand"], t["exons"]),
            cds_g_start=t["cds_g_start"], cds_g_end=t["cds_g_end"], is_canonical=True)
    return StaticProvider(txs, domains=fx["domains"], knowledge=fx["knowledge"])


# ---- resolved.genome_build echo (offline, no server deps) -----------------
def test_resolved_defaults_to_grch38():
    provider = _static_provider()
    r = annotate_fusion(provider, "EML4", "ALK", five_exon=13, three_exon=20)
    assert r["resolved"]["genome_build"] == "GRCh38"


def test_resolved_echoes_provider_assembly():
    provider = _static_provider()
    provider.assembly = "GRCh37"   # a provider bound to GRCh37 must surface that
    r = annotate_fusion(provider, "EML4", "ALK", five_exon=13, three_exon=20)
    assert r["resolved"]["genome_build"] == "GRCh37"


# ---- MCP provider is GRCh38-only and must not silently downgrade ----------
def test_mcp_provider_accepts_grch38():
    p = MCPDataProvider(lambda *a, **k: None, assembly="GRCh38")
    assert p.assembly == "GRCh38"
    p2 = MCPDataProvider(lambda *a, **k: None, assembly=None)
    assert p2.assembly == "GRCh38"


@pytest.mark.parametrize("build", ["GRCh37", "hg19", "37"])
def test_mcp_provider_rejects_non_grch38(build):
    with pytest.raises(NotImplementedError):
        MCPDataProvider(lambda *a, **k: None, assembly=build)


# ---- REST provider build->host mapping (needs the `requests` server dep) ---
def test_rest_provider_build_selection():
    pytest.importorskip("requests")
    from fusion_annotation.rest_provider import (
        RestDataProvider, normalize_assembly,
        ENSEMBL_BASE_GRCH38, ENSEMBL_BASE_GRCH37,
    )

    assert RestDataProvider().assembly == "GRCh38"
    assert RestDataProvider().ensembl_base == ENSEMBL_BASE_GRCH38
    assert RestDataProvider(assembly="GRCh37").ensembl_base == ENSEMBL_BASE_GRCH37
    # aliases + case-insensitivity
    assert RestDataProvider(assembly="hg19").assembly == "GRCh37"
    assert RestDataProvider(assembly="grch38").assembly == "GRCh38"
    assert normalize_assembly("HG38") == "GRCh38"
    assert normalize_assembly(None) == "GRCh38"
    assert normalize_assembly("") == "GRCh38"
    with pytest.raises(ValueError):
        normalize_assembly("GRCh99")
