import sys
from types import SimpleNamespace

import pytest

from fusion_annotation import gene_curation
from fusion_annotation.gene_curation import (
    FileCacheBackend,
    FusionCurationContext,
    NullCacheBackend,
    PubMedRecord,
    RedisCacheBackend,
)


@pytest.fixture(autouse=True)
def reset_cache_backend(monkeypatch):
    monkeypatch.setattr(gene_curation, "_CACHE_BACKEND", None)
    monkeypatch.setattr(gene_curation, "_CACHE_BACKEND_CONFIG", None)


def test_cache_backend_can_be_disabled(monkeypatch):
    monkeypatch.setenv("FUSION_GENE_CURATION_CACHE", "0")

    backend = gene_curation._cache_backend()

    assert isinstance(backend, NullCacheBackend)


def test_cache_backend_defaults_to_file(tmp_path, monkeypatch):
    monkeypatch.setenv("FUSION_GENE_CURATION_CACHE_DIR", str(tmp_path))
    monkeypatch.delenv("FUSION_GENE_CURATION_CACHE_BACKEND", raising=False)
    backend = gene_curation._cache_backend()

    assert isinstance(backend, FileCacheBackend)
    backend.write("ns", {"gene": "ALK"}, {"value": 1})
    assert backend.read("ns", {"gene": "ALK"}) == {"value": 1}


def test_redis_cache_backend_round_trips_json(monkeypatch):
    store = {}
    expirations = {}

    class FakeRedisClient:
        def get(self, key):
            return store.get(key)

        def set(self, key, value, **kwargs):
            store[key] = value
            expirations[key] = kwargs.get("ex")

    class FakeRedis:
        @staticmethod
        def from_url(url):
            assert url == "redis://cache.example/0"
            return FakeRedisClient()

    monkeypatch.setitem(sys.modules, "redis", SimpleNamespace(Redis=FakeRedis))

    backend = RedisCacheBackend("redis://cache.example/0", "test-prefix")
    backend.write("gene-curation", {"gene": "ALK"}, {"cached": True}, ttl_seconds=60)

    assert backend.read("gene-curation", {"gene": "ALK"}) == {"cached": True}
    assert list(expirations.values()) == [60]


def test_cache_backend_reuses_redis_client_for_same_config(monkeypatch):
    clients = []

    class FakeRedis:
        @staticmethod
        def from_url(url):
            client = SimpleNamespace(url=url, get=lambda key: None)
            clients.append(client)
            return client

    monkeypatch.setitem(sys.modules, "redis", SimpleNamespace(Redis=FakeRedis))
    monkeypatch.setenv("FUSION_GENE_CURATION_CACHE_BACKEND", "redis")
    monkeypatch.setenv("FUSION_GENE_CURATION_REDIS_URL", "redis://cache.example/0")

    first = gene_curation._cache_backend()
    second = gene_curation._cache_backend()

    assert first is second
    assert len(clients) == 1


def test_retrieve_pubmed_records_caches_empty_results(tmp_path, monkeypatch):
    calls = 0

    class FakeResponse:
        def raise_for_status(self):
            return None

        def json(self):
            return {"esearchresult": {"idlist": []}}

    def fake_get(url, **kwargs):
        nonlocal calls
        calls += 1
        assert url == gene_curation.ESEARCH_URL
        return FakeResponse()

    monkeypatch.setenv("FUSION_GENE_CURATION_CACHE_DIR", str(tmp_path))
    monkeypatch.setattr(gene_curation.requests, "get", fake_get)

    first = gene_curation.retrieve_pubmed_records("NOEVIDENCE")
    second = gene_curation.retrieve_pubmed_records("NOEVIDENCE")

    assert first == []
    assert second == []
    assert calls == 1


def test_curate_gene_skips_model_when_pubmed_is_empty(monkeypatch):
    monkeypatch.setattr(gene_curation, "retrieve_pubmed_records", lambda *args, **kwargs: [])
    monkeypatch.setitem(sys.modules, "anthropic", None)

    result = gene_curation.curate_gene("NOEVIDENCE", anthropic_api_key="configured")

    assert result == {
        "gene": "NOEVIDENCE",
        "cancer_associated": None,
        "rationale": "No PubMed abstracts were retrieved for this gene.",
        "supporting_pmids": [],
        "retrieved_pmids": [],
        "fusion_contexts": [],
        "insufficient_evidence": True,
    }


def test_curate_gene_reuses_persistent_cache(tmp_path, monkeypatch):
    monkeypatch.setenv("FUSION_GENE_CURATION_CACHE_DIR", str(tmp_path))
    monkeypatch.setattr(
        gene_curation,
        "retrieve_pubmed_records",
        lambda *args, **kwargs: [
            PubMedRecord(
                pmid="123",
                title="ALK fusion evidence",
                abstract="ALK fusions are oncogenic in lung cancer.",
            )
        ],
    )

    calls = 0

    class FakeMessages:
        def create(self, **kwargs):
            nonlocal calls
            calls += 1
            return SimpleNamespace(
                content=[
                    SimpleNamespace(
                        type="text",
                        text=(
                            '{"gene":"ALK","cancer_associated":true,'
                            '"rationale":"ALK fusions are oncogenic.",'
                            '"supporting_pmids":["123"],'
                            '"retrieved_pmids":["123"],'
                            '"insufficient_evidence":false}'
                        ),
                    )
                ]
            )

    class FakeAnthropic:
        def __init__(self, api_key):
            self.api_key = api_key
            self.messages = FakeMessages()

    monkeypatch.setitem(
        sys.modules,
        "anthropic",
        SimpleNamespace(Anthropic=FakeAnthropic),
    )

    first = gene_curation.curate_gene("ALK", anthropic_api_key="configured")
    second = gene_curation.curate_gene("ALK", anthropic_api_key="configured")

    assert first == second
    assert first["supporting_pmids"] == ["123"]
    assert calls == 1


def test_curate_gene_prompt_requests_concise_curator_rationale(monkeypatch):
    monkeypatch.setenv("FUSION_GENE_CURATION_CACHE", "0")
    monkeypatch.setattr(
        gene_curation,
        "retrieve_pubmed_records",
        lambda *args, **kwargs: [
            PubMedRecord(
                pmid="123",
                title="ALK fusion evidence",
                abstract="ALK fusions are oncogenic in lung cancer.",
            )
        ],
    )
    seen = {}

    class FakeMessages:
        def create(self, **kwargs):
            seen["system"] = kwargs["system"]
            seen["prompt"] = kwargs["messages"][0]["content"]
            return SimpleNamespace(
                content=[
                    SimpleNamespace(
                        type="text",
                        text=(
                            '{"gene":"ALK","cancer_associated":true,'
                            '"rationale":"ALK has concise functional cancer evidence.",'
                            '"supporting_pmids":["123"],'
                            '"retrieved_pmids":["123"],'
                            '"insufficient_evidence":false}'
                        ),
                    )
                ]
            )

    class FakeAnthropic:
        def __init__(self, api_key):
            self.messages = FakeMessages()

    monkeypatch.setitem(
        sys.modules,
        "anthropic",
        SimpleNamespace(Anthropic=FakeAnthropic),
    )

    gene_curation.curate_gene("ALK", anthropic_api_key="configured")

    assert "1-2 short sentences" in seen["prompt"]
    assert "40-75 words" in seen["prompt"]
    assert "Do not enumerate every paper" in seen["prompt"]
    assert "fast curator review" in seen["system"]


def test_curate_gene_prompt_includes_genome_nexus_fusion_context(monkeypatch):
    monkeypatch.setenv("FUSION_GENE_CURATION_CACHE", "0")
    monkeypatch.setattr(
        gene_curation,
        "retrieve_pubmed_records",
        lambda *args, **kwargs: [
            PubMedRecord(
                pmid="123",
                title="ALK fusion evidence",
                abstract="ALK fusions are oncogenic in lung cancer.",
            )
        ],
    )
    seen = {}

    class FakeMessages:
        def create(self, **kwargs):
            seen["prompt"] = kwargs["messages"][0]["content"]
            return SimpleNamespace(
                content=[
                    SimpleNamespace(
                        type="text",
                        text=(
                            '{"gene":"ALK","cancer_associated":true,'
                            '"rationale":"ALK has fusion-specific evidence.",'
                            '"supporting_pmids":["123"],'
                            '"retrieved_pmids":["123"],'
                            '"insufficient_evidence":false}'
                        ),
                    )
                ]
            )

    class FakeAnthropic:
        def __init__(self, api_key):
            self.messages = FakeMessages()

    monkeypatch.setitem(
        sys.modules,
        "anthropic",
        SimpleNamespace(Anthropic=FakeAnthropic),
    )

    context = FusionCurationContext(
        gene="ALK",
        fusion="EML4::ALK",
        side="three_prime",
        partner_gene="EML4",
        three_transcript="ENST00000389048",
        three_exon="20",
        three_genomic="chr2:29446394",
        three_protein_breakpoint="p.1058",
        retained_domains=("Protein kinase domain (1116-1383)",),
        kinase_gene="ALK",
        kinase_gene_side="three_prime",
        kinase_domain_status="retained",
    )

    result = gene_curation.curate_gene(
        "ALK",
        anthropic_api_key="configured",
        fusion_contexts=[context],
    )

    assert "Genome Nexus fusion-position context" in seen["prompt"]
    assert "EML4::ALK" in seen["prompt"]
    assert "ENST00000389048" in seen["prompt"]
    assert "domain_status=retained" in seen["prompt"]
    assert result["fusion_contexts"][0]["fusion"] == "EML4::ALK"


def test_curate_fusion_genes_passes_token_controls(monkeypatch):
    class Fusion:
        five_gene = "EML4"
        three_gene = "ALK"

    seen = {}

    def fake_curate_gene(gene, **kwargs):
        seen[gene] = kwargs
        return {"gene": gene, "insufficient_evidence": True}

    monkeypatch.setenv("ANTHROPIC_API_KEY", "configured")
    monkeypatch.setenv("NCBI_API_KEY", "ncbi")
    monkeypatch.setenv("FUSION_GENE_CURATION_MODEL", "curation-model")
    monkeypatch.setenv("FUSION_GENE_CURATION_MAX_RESULTS", "3")
    monkeypatch.setenv("FUSION_GENE_CURATION_ABSTRACT_CHARS", "600")
    monkeypatch.setenv("FUSION_GENE_CURATION_WORKERS", "1")
    monkeypatch.setattr(gene_curation, "curate_gene", fake_curate_gene)

    result = gene_curation.curate_fusion_genes([Fusion()])

    assert result["genes"] == [
        {"gene": "ALK", "insufficient_evidence": True},
        {"gene": "EML4", "insufficient_evidence": True},
    ]
    assert seen["ALK"]["model"] == "curation-model"
    assert seen["ALK"]["ncbi_api_key"] == "ncbi"
    assert seen["ALK"]["max_results"] == 3
    assert seen["ALK"]["abstract_chars"] == 600
    assert seen["ALK"]["fusion_contexts"][0].fusion == "EML4::ALK"
