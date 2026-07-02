"""MCP tool wrapper — exposes fusion annotation as a single callable tool for
agentic / cBioPortal / Genome Nexus workflows.

The tool takes a compact fusion call (gene pair + exon boundaries) and returns the
full three-layer annotation as JSON. It is transport-agnostic: `annotate_fusion_tool`
is a plain function you can register with any MCP server framework (FastMCP,
the reference `mcp` SDK, or a custom dispatcher).

Example (FastMCP):

    from mcp.server.fastmcp import FastMCP
    from fusion_annotation.mcp_tool import annotate_fusion_tool, TOOL_SCHEMA

    app = FastMCP("fusion-annotation")

    @app.tool()
    def annotate_fusion(five_gene: str, three_gene: str,
                        five_exon: int, three_exon: int,
                        five_transcript: str = None, three_transcript: str = None) -> dict:
        # `host_mcp` is however your runtime exposes upstream MCP connectors
        return annotate_fusion_tool(host_mcp, five_gene, three_gene,
                                    five_exon, three_exon,
                                    five_transcript, three_transcript)

In a Claude Science repl cell the upstream connector callable is simply `host.mcp`.
"""
from __future__ import annotations
from typing import Callable, Optional

from .core import annotate_fusion
from .providers import MCPDataProvider


TOOL_NAME = "annotate_gene_fusion"

TOOL_DESCRIPTION = (
    "Annotate a gene fusion at the protein level. Given a 5' and 3' partner gene and "
    "the fused breakpoints, reconstruct the chimeric protein (frame, junction, "
    "hybrid codon, retained/lost domains), emit an HGVS.p-like junction string, and "
    "attach curated clinical knowledge (oncogenicity, therapies, evidence) for the "
    "categorical gene-pair fusion. Breakpoints may be given as exon numbers or, to "
    "pin the isoform unambiguously, as genomic positions."
)

TOOL_SCHEMA = {
    "name": TOOL_NAME,
    "description": TOOL_DESCRIPTION,
    "inputSchema": {
        "type": "object",
        "properties": {
            "five_gene": {"type": "string", "description": "5' partner gene symbol, e.g. EML4"},
            "three_gene": {"type": "string", "description": "3' partner gene symbol, e.g. ALK"},
            "five_exon": {"type": "integer", "description": "Last exon (1-based) contributed by the 5' partner. Give this or five_genomic."},
            "three_exon": {"type": "integer", "description": "First exon (1-based) contributed by the 3' partner. Give this or three_genomic."},
            "five_genomic": {"type": "string", "description": "Genomic breakpoint on the 5' partner (int, 'chr6:117324415', or HGVS 'g.117324415'). Pins the isoform; preferred over five_exon when both are given."},
            "three_genomic": {"type": "string", "description": "Genomic breakpoint on the 3' partner. Pins the isoform; preferred over three_exon when both are given."},
            "five_transcript": {"type": "string", "description": "Optional Ensembl transcript id for the 5' partner (defaults to canonical; echoed back under resolved)"},
            "three_transcript": {"type": "string", "description": "Optional Ensembl transcript id for the 3' partner"},
            "genome_build": {"type": "string", "enum": ["GRCh38", "GRCh37"], "description": "Genome assembly the coordinates/transcripts come from. Default GRCh38. Genomic breakpoints MUST match this build; echoed back under resolved.genome_build."},
        },
        "required": ["five_gene", "three_gene"],
    },
}


def annotate_fusion_tool(mcp: Callable,
                         five_gene: str, three_gene: str,
                         five_exon: Optional[int] = None, three_exon: Optional[int] = None,
                         five_transcript: Optional[str] = None,
                         three_transcript: Optional[str] = None,
                         five_genomic=None, three_genomic=None,
                         genome_build: str = "GRCh38",
                         species: str = "homo_sapiens") -> dict:
    """Backend for the MCP tool. `mcp` is a callable `mcp(server, method, **kwargs)`.

    Each partner's breakpoint is given as an exon number or a genomic position (the
    latter pins the isoform, interpreted against `genome_build`). Returns the
    three-layer annotation dict:
    {"interface": ..., "knowledge": ..., "resolved": ..., "warnings": ...}.

    NOTE: the MCP `genomes` connector is GRCh38-only, so a non-GRCh38
    `genome_build` raises here. The REST-backed public server (server/app.py)
    supports GRCh37.
    """
    provider = MCPDataProvider(mcp, species=species, assembly=genome_build)
    return annotate_fusion(
        provider, five_gene, three_gene,
        five_exon=five_exon, three_exon=three_exon,
        five_tx=five_transcript, three_tx=three_transcript,
        five_genomic=five_genomic, three_genomic=three_genomic)
