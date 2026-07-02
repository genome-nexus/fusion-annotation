"""Public MCP server for fusion-annotation.

Exposes the single `annotate_gene_fusion` tool over the streamable-HTTP MCP
transport, backed by `RestDataProvider` (direct calls to Ensembl, InterPro and
CIViC — no Claude Science dependency). Designed to run as a container on
Cloud Run and be added as a remote connector in Claude.ai / Claude Desktop.
"""
from __future__ import annotations

import os
import sys

from starlette.applications import Starlette
from starlette.responses import PlainTextResponse
from starlette.routing import Route

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from mcp.server.fastmcp import FastMCP  # noqa: E402
from mcp.server.transport_security import TransportSecuritySettings  # noqa: E402
from fusion_annotation.core import annotate_fusion  # noqa: E402
from fusion_annotation.rest_provider import RestDataProvider  # noqa: E402

# Host header allowlist for the MCP transport's DNS-rebinding protection.
# Set FUSION_ANNOTATION_ALLOWED_HOSTS to a comma-separated list (e.g. the
# Cloud Run service hostname) in production. Defaults cover local testing
# and leave the door open for any *.run.app hostname if unset.
_allowed_hosts_env = os.environ.get("FUSION_ANNOTATION_ALLOWED_HOSTS", "").strip()
ALLOWED_HOSTS = (
    [h.strip() for h in _allowed_hosts_env.split(",") if h.strip()]
    if _allowed_hosts_env else ["localhost", "127.0.0.1", "testserver"]
)

mcp = FastMCP(
    "fusion-annotation",
    instructions=(
        "Annotate a gene fusion at the protein level: reconstructs the chimeric "
        "protein from two partner genes and their exon breakpoints, reports "
        "reading frame / junction residue / retained-or-lost protein domains, "
        "emits an HGVS.p-like junction string, and attaches curated fusion "
        "knowledge (therapies, evidence, disease context) from CIViC. "
        "This is a research/informatics tool, not a diagnostic device."
    ),
    stateless_http=True,
    json_response=True,
    transport_security=TransportSecuritySettings(
        enable_dns_rebinding_protection=bool(_allowed_hosts_env),
        allowed_hosts=ALLOWED_HOSTS,
        allowed_origins=["*"],
    ),
)


@mcp.tool()
def annotate_gene_fusion(
    five_gene: str,
    three_gene: str,
    five_exon: int | None = None,
    three_exon: int | None = None,
    five_genomic: int | str | None = None,
    three_genomic: int | str | None = None,
    five_transcript: str | None = None,
    three_transcript: str | None = None,
    species: str = "homo_sapiens",
) -> dict:
    """Annotate a gene fusion at the protein level.

    Given a 5' and 3' partner gene and the fused breakpoints, reconstruct the
    chimeric protein (frame, junction, hybrid codon, retained/lost domains),
    emit an HGVS.p-like junction string, and attach curated clinical knowledge
    (therapies, evidence, disease context) for the categorical gene-pair fusion.

    Each partner's breakpoint may be given either as an exon number or as a
    genomic position. A genomic position pins the transcript isoform and removes
    exon-numbering ambiguity between overlapping isoforms; when both are supplied
    for a partner the genomic position wins. The transcript actually used and how
    each breakpoint was interpreted are echoed back under ``resolved``, and a
    ``warnings`` list flags a known oncogenic pair that comes back out-of-frame.

    Args:
        five_gene: 5' partner gene symbol, e.g. "EML4".
        three_gene: 3' partner gene symbol, e.g. "ALK".
        five_exon: last exon (1-based) contributed by the 5' partner.
        three_exon: first exon (1-based) contributed by the 3' partner.
        five_genomic: genomic breakpoint on the 5' partner — an int position, a
            "chr6:117324415" form, or an HGVS "g.117324415" term.
        three_genomic: genomic breakpoint on the 3' partner (same forms).
        five_transcript: optional Ensembl transcript id for the 5' partner
            (defaults to its canonical transcript; echoed back under ``resolved``).
        three_transcript: optional Ensembl transcript id for the 3' partner.
        species: Ensembl species (default "homo_sapiens").
    """
    provider = RestDataProvider(species=species)
    return annotate_fusion(
        provider, five_gene, three_gene,
        five_exon=five_exon, three_exon=three_exon,
        five_tx=five_transcript, three_tx=three_transcript,
        five_genomic=five_genomic, three_genomic=three_genomic)


async def healthz(request) -> PlainTextResponse:
    return PlainTextResponse("ok")


def build_app() -> Starlette:
    inner = mcp.streamable_http_app()
    inner.router.routes.append(Route("/healthz", healthz))
    return inner


app = build_app()

if __name__ == "__main__":
    import uvicorn
    port = int(os.environ.get("PORT", "8080"))
    uvicorn.run(app, host="0.0.0.0", port=port)
