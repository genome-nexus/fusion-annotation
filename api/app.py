"""Public REST API for fusion-annotation.

A thin FastAPI wrapper around the same annotation engine used by the MCP
server (``server/app.py``): ``fusion_annotation.core.annotate_fusion``, backed
by ``GenomeNexusDataProvider`` by default (set
``FUSION_ANNOTATION_PROVIDER=rest`` to fall back to the Ensembl-backed
``RestDataProvider``, same convention as the MCP server).

This is the backend for the ``web/`` React UI, but is also usable directly —
``GET /api/annotate`` accepts the fusion inputs as query params (so the full
request URL doubles as a shareable, stateless permalink: reopening it just
re-runs the annotation, nothing is persisted server-side) and
``POST /api/annotate`` accepts the same fields as a JSON body for
programmatic callers. Interactive OpenAPI docs are served at ``/api/docs``.

Designed to run as a container on Cloud Run, alongside (not replacing) the
MCP server.
"""
from __future__ import annotations

import os
import sys
from typing import Optional

import requests
from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import PlainTextResponse
from pydantic import BaseModel, Field
from slowapi import Limiter, _rate_limit_exceeded_handler
from slowapi.errors import RateLimitExceeded
from slowapi.util import get_remote_address

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from fusion_annotation.core import annotate_fusion  # noqa: E402
from fusion_annotation.provider_factory import make_provider  # noqa: E402


# Per-IP rate limit, e.g. "30/minute" — protects the public deployment from
# abuse without requiring API keys/accounts. Override in production via env var.
RATE_LIMIT = os.environ.get("FUSION_ANNOTATION_RATE_LIMIT", "30/minute")

# Comma-separated list of allowed CORS origins for the web/ SPA. Defaults to
# "*" (open) since this is a public, unauthenticated demo API; lock this down
# to the deployed SPA's origin in production via env var.
_cors_origins_env = os.environ.get("FUSION_ANNOTATION_CORS_ORIGINS", "*").strip()
CORS_ORIGINS = [o.strip() for o in _cors_origins_env.split(",") if o.strip()] or ["*"]

limiter = Limiter(key_func=get_remote_address)

app = FastAPI(
    title="fusion-annotation API",
    description=(
        "Annotate a gene fusion at the protein level: reconstructs the chimeric "
        "protein from two partner genes and their exon breakpoints, reports "
        "reading frame / junction residue / retained-or-lost protein domains, "
        "emits an HGVS.p-like junction string, and attaches curated fusion "
        "knowledge (therapies, evidence, disease context) from CIViC. "
        "This is a research/informatics tool, not a diagnostic device."
    ),
    version="0.1.0",
    docs_url="/api/docs",
    openapi_url="/api/openapi.json",
)
app.state.limiter = limiter
app.add_exception_handler(RateLimitExceeded, _rate_limit_exceeded_handler)
app.add_middleware(
    CORSMiddleware,
    allow_origins=CORS_ORIGINS,
    allow_methods=["GET", "POST"],
    allow_headers=["*"],
)


class AnnotateRequest(BaseModel):
    """Same inputs as the MCP `annotate_gene_fusion` tool (server/app.py)."""

    five_gene: str = Field(..., description="5' partner gene symbol, e.g. EML4")
    three_gene: str = Field(..., description="3' partner gene symbol, e.g. ALK")
    five_exon: Optional[int] = Field(
        None, description="Last exon (1-based) contributed by the 5' partner. Give this or five_genomic.")
    three_exon: Optional[int] = Field(
        None, description="First exon (1-based) contributed by the 3' partner. Give this or three_genomic.")
    five_genomic: Optional[str] = Field(
        None, description="Genomic breakpoint on the 5' partner (int, 'chr6:117324415', or HGVS 'g.117324415'). "
                           "Pins the isoform; preferred over five_exon when both are given.")
    three_genomic: Optional[str] = Field(
        None, description="Genomic breakpoint on the 3' partner. Pins the isoform; preferred over three_exon "
                           "when both are given.")
    five_transcript: Optional[str] = Field(
        None, description="Optional Ensembl transcript id for the 5' partner (defaults to canonical).")
    three_transcript: Optional[str] = Field(
        None, description="Optional Ensembl transcript id for the 3' partner.")
    genome_build: str = Field(
        "GRCh38", description="Genome assembly the coordinates/transcripts come from. GRCh38 (default) or GRCh37.")
    species: str = Field(
        "homo_sapiens", description="Species identifier. Non-human species always use RestDataProvider.")


class AnnotateResponse(BaseModel):
    interface: dict
    knowledge: dict
    resolved: dict
    warnings: list[str]


class BatchAnnotateRequest(BaseModel):
    fusions: list[AnnotateRequest] = Field(
        ...,
        min_length=1,
        description="Fusion annotations to run in one request.",
    )


class BatchAnnotateItemResult(BaseModel):
    input: AnnotateRequest
    result: Optional[AnnotateResponse] = None
    error: Optional[str] = None


class BatchAnnotateResponse(BaseModel):
    results: list[BatchAnnotateItemResult]


def _annotate_with_provider(provider, params: AnnotateRequest) -> dict:
    return annotate_fusion(
        provider, params.five_gene, params.three_gene,
        five_exon=params.five_exon, three_exon=params.three_exon,
        five_tx=params.five_transcript, three_tx=params.three_transcript,
        five_genomic=params.five_genomic, three_genomic=params.three_genomic)


def _run_annotation(params: AnnotateRequest) -> dict:
    """Shared handler for GET/POST /api/annotate — builds a provider and calls
    the core engine, translating engine-level errors into HTTP responses."""
    try:
        provider = make_provider(species=params.species, assembly=params.genome_build)
        return _annotate_with_provider(provider, params)
    except (ValueError, KeyError) as exc:
        # Bad/unresolvable input: unknown gene, non-coding exon, breakpoint
        # that doesn't map into the CDS, etc.
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except requests.exceptions.RequestException as exc:
        # Upstream annotation source (Genome Nexus/Ensembl/InterPro/CIViC)
        # unreachable or erroring — not the caller's fault.
        raise HTTPException(status_code=502, detail=f"upstream annotation source error: {exc}") from exc


@app.get("/health")
async def health() -> PlainTextResponse:
    return PlainTextResponse("ok")


@app.get("/api/annotate", response_model=AnnotateResponse)
@limiter.limit(RATE_LIMIT)
def annotate_get(
    request: Request,  # noqa: ARG001 - required by slowapi's Limiter
    five_gene: str = Query(..., description="5' partner gene symbol, e.g. EML4"),
    three_gene: str = Query(..., description="3' partner gene symbol, e.g. ALK"),
    five_exon: Optional[int] = Query(None),
    three_exon: Optional[int] = Query(None),
    five_genomic: Optional[str] = Query(None),
    three_genomic: Optional[str] = Query(None),
    five_transcript: Optional[str] = Query(None),
    three_transcript: Optional[str] = Query(None),
    genome_build: str = Query("GRCh38"),
    species: str = Query("homo_sapiens"),
) -> dict:
    """Annotate a gene fusion. Inputs as query params, so this URL is a
    shareable permalink — reopening it re-runs the annotation (nothing is
    stored server-side).

    Plain ``def``, not ``async def``: annotate_fusion/make_provider do
    blocking network I/O via `requests`, so a sync endpoint lets FastAPI run
    it in its worker thread pool instead of blocking the event loop.
    """
    params = AnnotateRequest(
        five_gene=five_gene, three_gene=three_gene,
        five_exon=five_exon, three_exon=three_exon,
        five_genomic=five_genomic, three_genomic=three_genomic,
        five_transcript=five_transcript, three_transcript=three_transcript,
        genome_build=genome_build, species=species)
    return _run_annotation(params)


@app.post("/api/annotate", response_model=AnnotateResponse)
@limiter.limit(RATE_LIMIT)
def annotate_post(request: Request, params: AnnotateRequest) -> dict:  # noqa: ARG001
    """Annotate a gene fusion. Same fields as GET /api/annotate, as a JSON body.

    Plain ``def`` for the same reason as annotate_get above.
    """
    return _run_annotation(params)


@app.post("/api/annotate/batch", response_model=BatchAnnotateResponse)
@limiter.limit(RATE_LIMIT)
def annotate_batch(request: Request, params: BatchAnnotateRequest) -> BatchAnnotateResponse:  # noqa: ARG001
    """Annotate multiple fusions in one request.

    The batch path reuses one provider instance per request so repeated Genome
    Nexus/CIViC setup work is not repeated for every row. Individual bad inputs
    are returned as per-item errors instead of failing the whole batch.
    """
    try:
        provider = make_provider(
            species=params.fusions[0].species,
            assembly=params.fusions[0].genome_build,
        )
    except requests.exceptions.RequestException as exc:
        raise HTTPException(status_code=502, detail=f"upstream annotation source error: {exc}") from exc

    results = []
    for item in params.fusions:
        try:
            if item.species != params.fusions[0].species or item.genome_build != params.fusions[0].genome_build:
                provider = make_provider(species=item.species, assembly=item.genome_build)
            result = _annotate_with_provider(provider, item)
            results.append(BatchAnnotateItemResult(input=item, result=result))
        except (ValueError, KeyError) as exc:
            results.append(BatchAnnotateItemResult(input=item, error=str(exc)))
        except requests.exceptions.RequestException as exc:
            results.append(
                BatchAnnotateItemResult(
                    input=item,
                    error=f"upstream annotation source error: {exc}",
                )
            )

    return BatchAnnotateResponse(results=results)


if __name__ == "__main__":
    import uvicorn
    port = int(os.environ.get("PORT", "8080"))
    # proxy_headers + forwarded_allow_ips make uvicorn trust X-Forwarded-For
    # from its (only) ingress. Required for get_remote_address-based rate
    # limiting to actually key on the real client IP rather than Cloud Run's
    # internal proxy address — without this every caller behind the proxy
    # would share one rate-limit bucket. "*" is the standard setting Cloud
    # Run's own docs recommend, since the container is never reachable except
    # through that trusted proxy. Override via FUSION_ANNOTATION_TRUSTED_PROXIES
    # (comma-separated) if deploying behind a different/untrusted-by-default setup.
    forwarded_allow_ips = os.environ.get("FUSION_ANNOTATION_TRUSTED_PROXIES", "*")
    uvicorn.run(app, host="0.0.0.0", port=port,
                proxy_headers=True, forwarded_allow_ips=forwarded_allow_ips)
