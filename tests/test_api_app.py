"""Tests for the public REST API (api/app.py) — offline, no network.

Uses StaticProvider seeded from the EML4::ALK fixture (same one
tests/test_eml4_alk.py validates against Ensembl primary data), monkeypatched
in place of make_provider() so these tests never hit a real annotation
source.
"""
import importlib.util
import json
import os
import sys

import pytest

# api/app.py pulls in the api extra (fastapi, slowapi, requests), which isn't
# installed for the zero-dep core test matrix — skip cleanly there, same
# convention as tests/test_server_app.py for the MCP server.
pytest.importorskip("fastapi")
pytest.importorskip("slowapi")

sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir, "src"))

from fastapi.testclient import TestClient  # noqa: E402

from fusion_annotation import Transcript, build_exon_cds_map  # noqa: E402
from fusion_annotation.providers import StaticProvider  # noqa: E402

# Loaded via importlib (not sys.path + `import app`) under a unique module
# name: server/app.py is also a top-level module literally named `app`, and
# tests/test_server_app.py imports it the same way — a plain `import app`
# here would collide with it in sys.modules depending on collection order.
_API_APP_PATH = os.path.join(os.path.dirname(__file__), os.pardir, "api", "app.py")
_spec = importlib.util.spec_from_file_location("fusion_annotation_api_app", _API_APP_PATH)
api_app = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = api_app
_spec.loader.exec_module(api_app)

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
def client(static_provider, monkeypatch):
    monkeypatch.setattr(api_app, "make_provider", lambda species, assembly: static_provider)
    return TestClient(api_app.app)


def test_health(client):
    r = client.get("/health")
    assert r.status_code == 200
    assert r.text == "ok"


def test_annotate_get_eml4_alk(client):
    r = client.get("/api/annotate", params={
        "five_gene": "EML4", "three_gene": "ALK", "five_exon": 13, "three_exon": 20})
    assert r.status_code == 200
    body = r.json()
    assert body["interface"]["in_frame"] is True
    assert body["interface"]["categorical_key"] == "EML4::ALK"
    assert "hgvsp_like" in body["interface"]
    assert body["knowledge"]["oncogenic"] == "Oncogenic"
    assert body["resolved"]["five"]["gene"] == "EML4"


def test_annotate_post_matches_get(client):
    r = client.post("/api/annotate", json={
        "five_gene": "EML4", "three_gene": "ALK", "five_exon": 13, "three_exon": 20})
    assert r.status_code == 200
    assert r.json()["interface"]["categorical_key"] == "EML4::ALK"


def test_annotate_batch_returns_per_fusion_results(client):
    r = client.post("/api/annotate/batch", json={
        "fusions": [
            {"five_gene": "EML4", "three_gene": "ALK", "five_exon": 13, "three_exon": 20},
            {"five_gene": "NOPE_NOT_A_GENE", "three_gene": "ALK", "five_exon": 13, "three_exon": 20},
        ]
    })

    assert r.status_code == 200
    body = r.json()
    assert len(body["results"]) == 2
    assert body["results"][0]["result"]["interface"]["categorical_key"] == "EML4::ALK"
    assert body["results"][0]["error"] is None
    assert body["results"][1]["result"] is None
    assert "NOPE_NOT_A_GENE" in body["results"][1]["error"]


def test_annotate_missing_required_field(client):
    r = client.get("/api/annotate", params={"five_gene": "EML4"})
    assert r.status_code == 422  # FastAPI request validation, not our handler


def test_annotate_unknown_gene_is_400(client):
    r = client.get("/api/annotate", params={
        "five_gene": "NOPE_NOT_A_GENE", "three_gene": "ALK", "five_exon": 13, "three_exon": 20})
    assert r.status_code == 400
    assert "NOPE_NOT_A_GENE" in r.json()["detail"]


def test_annotate_no_breakpoint_given_is_400(client):
    r = client.get("/api/annotate", params={"five_gene": "EML4", "three_gene": "ALK"})
    assert r.status_code == 400


def test_cors_headers_present(client):
    r = client.get("/api/annotate", params={
        "five_gene": "EML4", "three_gene": "ALK", "five_exon": 13, "three_exon": 20},
        headers={"Origin": "https://example.org"})
    assert r.headers.get("access-control-allow-origin") == "*"


def test_rate_limit_returns_429(static_provider, monkeypatch):
    # RATE_LIMIT is read from the env at module import time, so load a fresh
    # module instance (distinct sys.modules key) with a tiny limit to exercise
    # the 429 path without waiting a minute or issuing 30+ requests.
    monkeypatch.setenv("FUSION_ANNOTATION_RATE_LIMIT", "2/minute")
    spec = importlib.util.spec_from_file_location("fusion_annotation_api_app_ratelimited", _API_APP_PATH)
    limited_app = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = limited_app
    spec.loader.exec_module(limited_app)
    monkeypatch.setattr(limited_app, "make_provider", lambda species, assembly: static_provider)

    limited_client = TestClient(limited_app.app)
    params = {"five_gene": "EML4", "three_gene": "ALK", "five_exon": 13, "three_exon": 20}
    assert limited_client.get("/api/annotate", params=params).status_code == 200
    assert limited_client.get("/api/annotate", params=params).status_code == 200
    r = limited_client.get("/api/annotate", params=params)
    assert r.status_code == 429
