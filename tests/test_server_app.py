"""Guards against a regression that took down the deployed MCP server: the
FUSION_ANNOTATION_ALLOWED_HOSTS env var was set to a full URL (scheme + path)
instead of a bare hostname, so the DNS-rebinding Host-header check never
matched the real Host header and rejected 100% of traffic with 421.
server.app.normalize_allowed_host() defensively extracts the bare host[:port]
from whatever is configured.

Uses StaticProvider seeded from the EML4::ALK fixture (same convention as
tests/test_api_app.py), monkeypatched in place of make_provider() so these
tests never hit a real annotation source.
"""
import json
import os
import sys

import pytest

# server/app.py pulls in the server extra (starlette, mcp, uvicorn), which isn't
# installed for the zero-dep core test matrix — skip cleanly there. The dedicated
# `test-server` CI job installs those deps so these assertions actually run.
pytest.importorskip("starlette")
pytest.importorskip("mcp")

sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir, "server"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir, "src"))

from fusion_annotation import Transcript, build_exon_cds_map  # noqa: E402
from fusion_annotation.providers import StaticProvider  # noqa: E402

import app  # noqa: E402
from app import normalize_allowed_host  # noqa: E402

FIXTURE = os.path.join(os.path.dirname(__file__), "fixtures", "eml4_alk_fixture.json")


@pytest.fixture(scope="module")
def static_provider():
    fx = json.load(open(FIXTURE))
    txs = {}
    for key, t in fx["transcripts"].items():
        exon_cds = build_exon_cds_map(t["strand"], t["exons"], t["cds_g_start"], t["cds_g_end"])
        txs[key] = Transcript(
            gene_symbol=t["gene_symbol"], gene_id=t["gene_id"],
            transcript_id=t["transcript_id"], strand=t["strand"],
            cds=t["cds"], protein=t["protein"], uniprot=t["uniprot"], exon_cds=exon_cds)
    return StaticProvider(txs, domains={k: v for k, v in fx["domains"].items()},
                          knowledge=fx["knowledge"])


@pytest.fixture()
def _use_static_provider(static_provider, monkeypatch):
    monkeypatch.setattr(app, "_make_provider", lambda species, assembly: static_provider)


@pytest.mark.parametrize("entry,expected", [
    ("my-service-abc123.a.run.app", "my-service-abc123.a.run.app"),
    ("https://my-service-abc123.a.run.app/mcp",
     "my-service-abc123.a.run.app"),
    ("http://foo.run.app", "foo.run.app"),
    ("foo.run.app/mcp/", "foo.run.app"),
    ("localhost:8080", "localhost:8080"),
    ("  foo.run.app  ", "foo.run.app"),
])
def test_normalize_allowed_host(entry, expected):
    assert normalize_allowed_host(entry) == expected


def test_annotate_gene_fusion_attaches_diagram_image(_use_static_provider):
    pytest.importorskip("matplotlib")

    from mcp.types import CallToolResult

    result = app.annotate_gene_fusion("EML4", "ALK", five_exon=13, three_exon=20)

    assert isinstance(result, CallToolResult)
    assert result.structuredContent["interface"]["five_gene"] == "EML4"
    types = [c.type for c in result.content]
    assert types == ["text", "image"]
    image = result.content[1]
    assert image.mimeType == "image/png"
    assert len(image.data) > 0


def test_annotate_gene_fusion_include_diagram_false_skips_image(_use_static_provider):
    from mcp.types import CallToolResult

    result = app.annotate_gene_fusion(
        "EML4", "ALK", five_exon=13, three_exon=20, include_diagram=False)

    assert isinstance(result, CallToolResult)
    assert [c.type for c in result.content] == ["text"]

