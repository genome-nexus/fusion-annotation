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

import base64
import logging
import os
import sys
from urllib.parse import urlparse

from starlette.responses import PlainTextResponse

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from mcp.server.fastmcp import FastMCP  # noqa: E402
from mcp.server.transport_security import TransportSecuritySettings  # noqa: E402
from mcp.types import CallToolResult, ImageContent, TextContent  # noqa: E402
from fusion_annotation.core import annotate_fusion  # noqa: E402
from fusion_annotation.provider_factory import make_provider as _make_provider  # noqa: E402

log = logging.getLogger(__name__)


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
    include_diagram: bool = True,
) -> CallToolResult:
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
        include_diagram: if True (default), attach a rendered PNG diagram of
            the 5' partner, 3' partner, and fusion protein domain tracks
            (same rendering as docs/fusion_domain_map.png and the web UI's
            domain diagram) alongside the JSON annotation. Set False to skip
            rendering and get a text/structured-only response (faster, no
            matplotlib import).
    """
    provider = _make_provider(species=species, assembly=genome_build)
    result = annotate_fusion(
        provider, five_gene, three_gene,
        five_exon=five_exon, three_exon=three_exon,
        five_tx=five_transcript, three_tx=three_transcript,
        five_genomic=five_genomic, three_genomic=three_genomic)

    web_link = _build_web_link(
        five_gene, three_gene, five_exon, three_exon, 
        five_genomic, three_genomic, five_transcript, three_transcript, 
        genome_build)

    content: list[TextContent | ImageContent] = [
        TextContent(
            type="text", 
            text=f"View in web UI: {web_link}\n\n{_json_dumps(result)}"
        )
    ]
    if include_diagram:
        try:
            from fusion_annotation.domain_diagram import render_domain_diagram_png
            png_bytes = render_domain_diagram_png(result)
            content.append(ImageContent(
                type="image",
                data=base64.b64encode(png_bytes).decode("ascii"),
                mimeType="image/png",
            ))
        except Exception:
            # Diagram rendering is best-effort — never fail the whole tool
            # call (e.g. missing matplotlib, or an unusual domain layout)
            # just because the image couldn't be produced.
            log.warning("domain diagram rendering failed for %s::%s", five_gene, three_gene,
                        exc_info=True)

    return CallToolResult(content=content, structuredContent=result)


@mcp.tool()
def curate_fusion(
    five_gene: str,
    three_gene: str,
    five_exon: int | None = None,
    three_exon: int | None = None,
    five_genomic: int | str | None = None,
    three_genomic: int | str | None = None,
    genome_build: str = "GRCh38",
    tumor_type: str | None = None,
    force_gene_curation: bool = False,
) -> CallToolResult:
    """Curate a gene fusion with literature and clinical knowledge.

    Retrieves PubMed literature, OncoKB gene data, and CIViC evidence for
    the fusion and its partner genes, then uses Claude to synthesise a
    structured curation summary covering:
    - Whether the fusion has been observed in the literature / in cancer
    - Functional oncogenicity and therapeutic response
    - Per-gene cancer role, mutation and expression profiles
    - Supporting PMIDs and high-impact journal flags (★)

    Results are cached locally (file-based by default) so repeat queries for
    the same fusion and tumor type skip the LLM call entirely.

    Requires ANTHROPIC_API_KEY to be set in the server environment.

    Args:
        five_gene: 5' partner gene symbol, e.g. "EML4".
        three_gene: 3' partner gene symbol, e.g. "ALK".
        five_exon: last exon (1-based) contributed by the 5' partner.
        three_exon: first exon (1-based) contributed by the 3' partner.
        five_genomic: genomic breakpoint on the 5' partner.
        three_genomic: genomic breakpoint on the 3' partner.
        genome_build: genome assembly — "GRCh38" (default) or "GRCh37".
        tumor_type: optional tumor/cancer type context (e.g. "lung adenocarcinoma")
            used to bias literature retrieval toward the disease under review.
        force_gene_curation: if True, run per-gene literature curation even
            when sufficient fusion-level literature was found (default False).
    """
    from types import SimpleNamespace
    from fusion_annotation.gene_curation import curate_fusion_genes, GeneCurationUnavailable

    fusion = SimpleNamespace(
        five_gene=five_gene,
        three_gene=three_gene,
        five_exon=five_exon,
        three_exon=three_exon,
        five_genomic=str(five_genomic) if five_genomic is not None else None,
        three_genomic=str(three_genomic) if three_genomic is not None else None,
        genome_build=genome_build,
        tumor_type=tumor_type,
    )
    try:
        result = curate_fusion_genes(
            [fusion],
            annotation_results=None,
            force_gene_curation=force_gene_curation,
        )
    except GeneCurationUnavailable as exc:
        return CallToolResult(
            content=[TextContent(type="text", text=f"Curation unavailable: {exc}")],
            isError=True,
        )
    return CallToolResult(content=[TextContent(type="text", text=_json_dumps(result))])


def _json_dumps(result: dict) -> str:
    import json
    return json.dumps(result, indent=2)


def _build_web_link(five_gene: str, three_gene: str,
                   five_exon: int | None, three_exon: int | None,
                   five_genomic: int | str | None, three_genomic: int | str | None,
                   five_transcript: str | None, three_transcript: str | None,
                   genome_build: str) -> str:
    """Generate a shareable link to the production web UI with these parameters."""
    from urllib.parse import urlencode
    params = {
        "five_gene": five_gene,
        "three_gene": three_gene,
        "genome_build": genome_build,
    }
    if five_exon is not None:
        params["five_exon"] = str(five_exon)
    if three_exon is not None:
        params["three_exon"] = str(three_exon)
    if five_genomic is not None:
        params["five_genomic"] = str(five_genomic)
    if three_genomic is not None:
        params["three_genomic"] = str(three_genomic)
    if five_transcript is not None:
        params["five_transcript"] = five_transcript
    if three_transcript is not None:
        params["three_transcript"] = three_transcript
    
    query_string = urlencode(params)
    base_url = os.environ.get("FUSION_ANNOTATION_WEB_URL", 
                              "https://genome-nexus.github.io/fusion-annotation/")
    return f"{base_url}?{query_string}"


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
