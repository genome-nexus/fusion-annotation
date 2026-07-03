"""Public MCP server for fusion-annotation.

Exposes the single `annotate_gene_fusion` tool over the streamable-HTTP MCP
transport, backed by ``GenomeNexusDataProvider`` by default (Genome Nexus +
UCSC + CIViC — no Ensembl REST, no Claude Science dependency, ~1–2 s
per annotation). Set ``FUSION_ANNOTATION_PROVIDER=rest`` to fall back to the
legacy Ensembl-backed ``RestDataProvider``.

Designed to run as a container on Cloud Run and be added as a remote
connector in Claude.ai / Claude Desktop.
"""
from __future__ import annotations

import os
import sys
from urllib.parse import urlparse

from starlette.responses import PlainTextResponse

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from mcp.server.fastmcp import FastMCP  # noqa: E402
from mcp.server.transport_security import TransportSecuritySettings  # noqa: E402
from fusion_annotation.core import annotate_fusion  # noqa: E402
from fusion_annotation.gn_provider import GenomeNexusDataProvider  # noqa: E402
from fusion_annotation.rest_provider import RestDataProvider  # noqa: E402


def _make_provider(species: str, assembly: str):
    """Instantiate the configured data provider.

    Reads ``FUSION_ANNOTATION_PROVIDER`` from the environment:
      - unset / "gn" / "genome_nexus" → GenomeNexusDataProvider (default, fast)
      - "rest" / "ensembl"            → RestDataProvider (fallback, slower)

    Non-human species always use ``RestDataProvider`` because
    ``GenomeNexusDataProvider`` only covers *Homo sapiens*.
    """
    backend = os.environ.get("FUSION_ANNOTATION_PROVIDER", "gn").strip().lower()
    if backend in ("rest", "ensembl") or (species or "").lower() not in ("homo_sapiens", "human"):
        return RestDataProvider(species=species, assembly=assembly)
    return GenomeNexusDataProvider(assembly=assembly)


def normalize_allowed_host(entry: str) -> str:
    """Extract a bare ``host[:port]`` from an allowlist entry.

    Defensive against a full URL (scheme and/or path) being pasted in by
    mistake -- e.g. "https://foo.run.app/mcp" instead of "foo.run.app" --
    which would otherwise never string-match the bare Host header Starlette
    sees on incoming requests and would silently reject 100% of legitimate
    traffic with a 421 (this exact misconfiguration took the deployed server
    down: FUSION_ANNOTATION_ALLOWED_HOSTS was set to the full MCP endpoint
    URL instead of just its hostname).
    """
    s = entry.strip()
    if "//" not in s:
        s = "//" + s   # force urlparse to treat a bare "host[:port]" as a netloc, not a path
    return urlparse(s).netloc or entry.strip()


# Host header allowlist for the MCP transport's DNS-rebinding protection.
# Set FUSION_ANNOTATION_ALLOWED_HOSTS to a comma-separated list of bare
# hostnames (e.g. the Cloud Run service hostname, NOT a full URL) in
# production. Defaults cover local testing and leave the door open for any
# *.run.app hostname if unset.
_allowed_hosts_env = os.environ.get("FUSION_ANNOTATION_ALLOWED_HOSTS", "").strip()
ALLOWED_HOSTS = (
    [normalize_allowed_host(h) for h in _allowed_hosts_env.split(",") if h.strip()]
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
    genome_build: str = "GRCh38",
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
            "chr6:117324415" form, or an HGVS "g.117324415" term. Interpreted
            against `genome_build`.
        three_genomic: genomic breakpoint on the 3' partner (same forms).
        five_transcript: optional Ensembl transcript id for the 5' partner
            (defaults to its canonical transcript; echoed back under ``resolved``).
        three_transcript: optional Ensembl transcript id for the 3' partner.
        genome_build: human genome assembly the coordinates/transcripts come from
            — "GRCh38" (default) or "GRCh37" (aliases hg38/hg19). Genomic
            breakpoints MUST match this build; it is echoed back under
            ``resolved.genome_build``.
        species: species identifier (default "homo_sapiens"). Currently only
            human is supported by the default ``GenomeNexusDataProvider``; the
            fallback ``RestDataProvider`` (``FUSION_ANNOTATION_PROVIDER=rest``)
            forwards this to Ensembl REST.
    """
    provider = _make_provider(species=species, assembly=genome_build)
    return annotate_fusion(
        provider, five_gene, three_gene,
        five_exon=five_exon, three_exon=three_exon,
        five_tx=five_transcript, three_tx=three_transcript,
        five_genomic=five_genomic, three_genomic=three_genomic)


async def health(request) -> PlainTextResponse:
    return PlainTextResponse("ok")


class HealthzMiddleware:
    """Intercept /health before FastMCP's middleware stack sees it.

    FastMCP's DNS-rebinding protection and session manager run as ASGI
    middleware around the inner app; appending /health to the inner router
    means it gets caught by that middleware first and never reaches the route.
    Wrapping at the ASGI level ensures /health is handled unconditionally.

    Note: Cloud Run intercepts GET /healthz at the infrastructure level and
    returns a Google-branded 404 before the request reaches the container, so
    we use /health instead.
    """

    def __init__(self, app):
        self.app = app

    async def __call__(self, scope, receive, send):
        if scope.get("type") == "http" and scope.get("path", "").rstrip("/") == "/health":
            await PlainTextResponse("ok")(scope, receive, send)
        else:
            await self.app(scope, receive, send)


def build_app():
    return HealthzMiddleware(mcp.streamable_http_app())


app = build_app()

if __name__ == "__main__":
    import uvicorn
    port = int(os.environ.get("PORT", "8080"))
    uvicorn.run(app, host="0.0.0.0", port=port)
