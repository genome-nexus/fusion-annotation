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
    "the fused exon boundaries, reconstruct the chimeric protein (frame, junction, "
    "hybrid codon, retained/lost domains), emit an HGVS.p-like junction string, and "
    "attach curated clinical knowledge (oncogenicity, therapies, evidence) for the "
    "categorical gene-pair fusion."
)

TOOL_SCHEMA = {
    "name": TOOL_NAME,
    "description": TOOL_DESCRIPTION,
    "inputSchema": {
        "type": "object",
        "properties": {
            "five_gene": {"type": "string", "description": "5' partner gene symbol, e.g. EML4"},
            "three_gene": {"type": "string", "description": "3' partner gene symbol, e.g. ALK"},
            "five_exon": {"type": "integer", "description": "Last exon (1-based) contributed by the 5' partner"},
            "three_exon": {"type": "integer", "description": "First exon (1-based) contributed by the 3' partner"},
            "five_transcript": {"type": "string", "description": "Optional Ensembl transcript id for the 5' partner"},
            "three_transcript": {"type": "string", "description": "Optional Ensembl transcript id for the 3' partner"},
        },
        "required": ["five_gene", "three_gene", "five_exon", "three_exon"],
    },
}


def annotate_fusion_tool(mcp: Callable,
                         five_gene: str, three_gene: str,
                         five_exon: int, three_exon: int,
                         five_transcript: Optional[str] = None,
                         three_transcript: Optional[str] = None,
                         species: str = "homo_sapiens") -> dict:
    """Backend for the MCP tool. `mcp` is a callable `mcp(server, method, **kwargs)`.

    Returns the three-layer annotation dict: {"interface": ..., "knowledge": ...}.
    """
    provider = MCPDataProvider(mcp, species=species)
    return annotate_fusion(
        provider, five_gene, three_gene,
        five_exon=five_exon, three_exon=three_exon,
        five_tx=five_transcript, three_tx=three_transcript)
