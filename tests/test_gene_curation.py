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
        status_code = 200
        headers = {}

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


def test_retrieve_pubmed_records_preserves_mixed_content_xml(tmp_path, monkeypatch):
    class SearchResponse:
        status_code = 200
        headers = {}
        text = ""

        def raise_for_status(self):
            return None

        def json(self):
            return {"esearchresult": {"idlist": ["123"]}}

    class FetchResponse:
        status_code = 200
        headers = {}
        text = """
        <PubmedArticleSet>
          <PubmedArticle>
            <MedlineCitation>
              <PMID>123</PMID>
              <Article>
                <ArticleTitle>ALK <i>fusion</i> evidence</ArticleTitle>
                <Abstract>
                  <AbstractText>Kinase <b>domain</b> retained.</AbstractText>
                  <AbstractText>Functional cancer evidence.</AbstractText>
                </Abstract>
              </Article>
            </MedlineCitation>
          </PubmedArticle>
        </PubmedArticleSet>
        """

        def raise_for_status(self):
            return None

    def fake_get(url, **kwargs):
        if url == gene_curation.ESEARCH_URL:
            return SearchResponse()
        assert url == gene_curation.EFETCH_URL
        return FetchResponse()

    monkeypatch.setenv("FUSION_GENE_CURATION_NCBI_MIN_INTERVAL_SECONDS", "0")
    monkeypatch.setenv("FUSION_GENE_CURATION_CACHE_DIR", str(tmp_path))
    monkeypatch.setattr(gene_curation.requests, "get", fake_get)

    records = gene_curation.retrieve_pubmed_records("ALK")

    assert records == [
        PubMedRecord(
            pmid="123",
            title="ALK fusion evidence",
            abstract="Kinase domain retained. Functional cancer evidence.",
        )
    ]


def test_retrieve_pubmed_records_retries_ncbi_429(monkeypatch):
    calls = []

    class RateLimitedResponse:
        status_code = 429
        headers = {}
        text = ""

        def raise_for_status(self):
            raise AssertionError("429 should be retried before raise_for_status")

    class SearchResponse:
        status_code = 200
        headers = {}
        text = ""

        def raise_for_status(self):
            return None

        def json(self):
            return {"esearchresult": {"idlist": ["123"]}}

    class FetchResponse:
        status_code = 200
        headers = {}
        text = """
        <PubmedArticleSet>
          <PubmedArticle>
            <MedlineCitation>
              <PMID>123</PMID>
              <Article>
                <ArticleTitle>ALK fusion evidence</ArticleTitle>
                <Abstract>
                  <AbstractText>ALK fusions are oncogenic.</AbstractText>
                </Abstract>
              </Article>
            </MedlineCitation>
          </PubmedArticle>
        </PubmedArticleSet>
        """

        def raise_for_status(self):
            return None

    def fake_get(url, **kwargs):
        calls.append(url)
        if url == gene_curation.ESEARCH_URL and calls.count(gene_curation.ESEARCH_URL) == 1:
            return RateLimitedResponse()
        if url == gene_curation.ESEARCH_URL:
            return SearchResponse()
        assert url == gene_curation.EFETCH_URL
        return FetchResponse()

    monkeypatch.setenv("FUSION_GENE_CURATION_CACHE", "0")
    monkeypatch.setattr(gene_curation.requests, "get", fake_get)
    client = gene_curation.NcbiClient(
        min_interval_seconds=0,
        max_retries=1,
        backoff_seconds=0,
    )

    records = gene_curation.retrieve_pubmed_records("ALK", ncbi_client=client)

    assert calls == [
        gene_curation.ESEARCH_URL,
        gene_curation.ESEARCH_URL,
        gene_curation.EFETCH_URL,
    ]
    assert records[0].pmid == "123"


def test_curate_fusion_genes_reports_pubmed_rate_limit_as_gene_error(monkeypatch):
    class Fusion:
        five_gene = "EML4"
        three_gene = "ALK"

    def fake_curate_fusion(*args, **kwargs):
        return {
            "fusion": "EML4::ALK",
            "fusion_literature_identified": False,
            "insufficient_evidence": True,
        }

    def fake_retrieve_pubmed_records(*args, **kwargs):
        raise gene_curation.PubMedRateLimitError("NCBI PubMed rate limit was reached.")

    monkeypatch.setenv("ANTHROPIC_API_KEY", "configured")
    monkeypatch.setenv("FUSION_GENE_CURATION_WORKERS", "1")
    monkeypatch.setattr(gene_curation, "curate_fusion", fake_curate_fusion)
    monkeypatch.setattr(gene_curation, "retrieve_pubmed_records", fake_retrieve_pubmed_records)

    result = gene_curation.curate_fusion_genes([Fusion()])

    assert result["genes"][0]["insufficient_evidence"] is True
    assert "NCBI PubMed rate limit" in result["genes"][0]["error"]


def test_unique_genes_from_fusions_accepts_dicts_and_missing_values():
    assert gene_curation.unique_genes_from_fusions([
        {"five_gene": "EML4", "three_gene": "ALK"},
        {"five_gene": "ALK"},
        None,
        SimpleNamespace(five_gene="CD74", three_gene="ROS1"),
    ]) == ["EML4", "ALK", "CD74", "ROS1"]


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

def test_curate_gene_uses_oncokb_before_pubmed_or_model(monkeypatch):
    def fake_retrieve_pubmed_records(*args, **kwargs):
        raise AssertionError("PubMed should not be queried when OncoKB has a curated gene record")

    monkeypatch.setattr(gene_curation, "retrieve_pubmed_records", fake_retrieve_pubmed_records)
    monkeypatch.setitem(sys.modules, "anthropic", None)

    result = gene_curation.curate_gene(
        "ALK",
        anthropic_api_key="",
        oncokb_genes_by_symbol={
            "ALK": {
                "hugoSymbol": "ALK",
                "geneType": "ONCOGENE",
                "summary": "ALK is a receptor tyrosine kinase altered in multiple cancers.",
                "highestSensitiveLevel": "LEVEL_1",
            }
        },
    )

    assert result["curation_source"] == "OncoKB"
    assert result["cancer_associated"] is True
    assert result["oncokb_gene_type"] == "ONCOGENE"
    assert result["oncokb_highest_sensitive_level"] == "LEVEL_1"
    assert result["supporting_pmids"] == []
    assert "OncoKB curates ALK as oncogene" in result["rationale"]


def test_curate_fusion_genes_fetches_oncokb_once_for_gene_fallback(monkeypatch):
    class Fusion:
        five_gene = "EML4"
        three_gene = "ALK"

    calls = 0

    def fake_fetch_oncokb_curated_gene_index(api_token):
        nonlocal calls
        calls += 1
        assert api_token == "oncokb-token"
        return {
            "EML4": {
                "hugoSymbol": "EML4",
                "geneType": "NEITHER",
                "summary": "EML4 is curated by OncoKB.",
            },
            "ALK": {
                "hugoSymbol": "ALK",
                "geneType": "ONCOGENE",
                "summary": "ALK is a receptor tyrosine kinase altered in cancer.",
            },
        }

    def fake_curate_fusion(fusion, **kwargs):
        return {
            "fusion": fusion,
            "fusion_literature_identified": False,
            "insufficient_evidence": True,
        }

    def fake_retrieve_pubmed_records(*args, **kwargs):
        raise AssertionError("PubMed should not be queried when OncoKB covers fallback genes")

    monkeypatch.setenv("ANTHROPIC_API_KEY", "configured")
    monkeypatch.setenv("ONCOKB_API_TOKEN", "oncokb-token")
    monkeypatch.setenv("FUSION_GENE_CURATION_WORKERS", "1")
    monkeypatch.setattr(gene_curation, "curate_fusion", fake_curate_fusion)
    monkeypatch.setattr(gene_curation, "fetch_oncokb_curated_gene_index", fake_fetch_oncokb_curated_gene_index)
    monkeypatch.setattr(gene_curation, "retrieve_pubmed_records", fake_retrieve_pubmed_records)

    result = gene_curation.curate_fusion_genes([Fusion()])

    assert calls == 1
    assert [item["gene"] for item in result["genes"]] == ["ALK", "EML4"]
    assert [item["curation_source"] for item in result["genes"]] == ["OncoKB", "OncoKB"]


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


def test_fetch_oncokb_curated_gene_index_reuses_persistent_cache(tmp_path, monkeypatch):
    monkeypatch.setenv("FUSION_GENE_CURATION_CACHE_DIR", str(tmp_path))
    calls = 0

    class FakeResponse:
        status_code = 200

        def raise_for_status(self):
            return None

        def json(self):
            return [
                {
                    "hugoSymbol": "ALK",
                    "geneType": "ONCOGENE",
                    "summary": "ALK is curated by OncoKB.",
                }
            ]

    def fake_get(url, **kwargs):
        nonlocal calls
        calls += 1
        assert url.endswith("/utils/allCuratedGenes")
        assert kwargs["params"]["includeEvidence"] == "true"
        assert kwargs["headers"]["Authorization"] == "Bearer oncokb-token"
        return FakeResponse()

    monkeypatch.setattr(gene_curation.requests, "get", fake_get)

    first = gene_curation.fetch_oncokb_curated_gene_index("oncokb-token")
    second = gene_curation.fetch_oncokb_curated_gene_index("oncokb-token")

    assert first == second
    assert first["ALK"]["geneType"] == "ONCOGENE"
    assert calls == 1


def test_curate_fusion_reuses_persistent_cache(tmp_path, monkeypatch):
    monkeypatch.setenv("FUSION_GENE_CURATION_CACHE_DIR", str(tmp_path))
    monkeypatch.setattr(
        gene_curation,
        "retrieve_pubmed_records_for_queries",
        lambda *args, **kwargs: [
            PubMedRecord(
                pmid="123",
                title="LMNA NTRK1 fusion evidence",
                abstract="LMNA::NTRK1 fusions are oncogenic in soft tissue tumors.",
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
                            '{"fusion":"LMNA::NTRK1","fusion_literature_identified":true,'
                            '"cancer_associated":true,'
                            '"rationale":"LMNA::NTRK1 has fusion-specific evidence.",'
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

    first = gene_curation.curate_fusion("LMNA::NTRK1", anthropic_api_key="configured")
    second = gene_curation.curate_fusion("LMNA::NTRK1", anthropic_api_key="configured")

    assert first == second
    assert first["supporting_pmids"] == ["123"]
    assert calls == 1


def test_curate_gene_parses_json_wrapped_in_markdown_fence(monkeypatch):
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

    class FakeMessages:
        def create(self, **kwargs):
            return SimpleNamespace(
                content=[
                    SimpleNamespace(
                        type="text",
                        text=(
                            "```json\n"
                            '{"gene":"ALK","cancer_associated":true,'
                            '"rationale":"ALK has functional cancer evidence.",'
                            '"supporting_pmids":["123"],'
                            '"retrieved_pmids":["123"],'
                            '"insufficient_evidence":false}'
                            "\n```"
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

    result = gene_curation.curate_gene("ALK", anthropic_api_key="configured")

    assert result["cancer_associated"] is True
    assert result["supporting_pmids"] == ["123"]


def test_curate_fusion_parses_json_from_unclosed_markdown_fence(monkeypatch):
    monkeypatch.setenv("FUSION_GENE_CURATION_CACHE", "0")
    monkeypatch.setattr(
        gene_curation,
        "retrieve_pubmed_records_for_queries",
        lambda *args, **kwargs: [
            PubMedRecord(
                pmid="123",
                title="EML4 ALK fusion evidence",
                abstract="EML4 ALK fusions are oncogenic in lung cancer.",
            )
        ],
    )

    class FakeMessages:
        def create(self, **kwargs):
            return SimpleNamespace(
                content=[
                    SimpleNamespace(
                        type="text",
                        text=(
                            "```json\n"
                            '{"fusion":"EML4::ALK","fusion_literature_identified":true,'
                            '"cancer_associated":true,'
                            '"rationale":"EML4::ALK has fusion-specific evidence.",'
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

    result = gene_curation.curate_fusion("EML4::ALK", anthropic_api_key="configured")

    assert result["fusion_literature_identified"] is True
    assert result["rationale"] == "EML4::ALK has fusion-specific evidence."
    assert result["supporting_pmids"] == ["123"]


def test_curate_fusion_does_not_store_raw_json_as_fallback_rationale(monkeypatch):
    monkeypatch.setenv("FUSION_GENE_CURATION_CACHE", "0")
    monkeypatch.setattr(
        gene_curation,
        "retrieve_pubmed_records_for_queries",
        lambda *args, **kwargs: [
            PubMedRecord(
                pmid="123",
                title="EML4 ALK fusion evidence",
                abstract="EML4 ALK fusions are oncogenic in lung cancer.",
            )
        ],
    )

    class FakeMessages:
        def create(self, **kwargs):
            return SimpleNamespace(
                content=[
                    SimpleNamespace(
                        type="text",
                        text='{"fusion":"EML4::ALK","rationale":"missing closing brace"',
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

    result = gene_curation.curate_fusion("EML4::ALK", anthropic_api_key="configured")

    assert result["insufficient_evidence"] is True
    assert result["rationale"].startswith("The curation model returned a malformed")
    assert '"fusion"' not in result["rationale"]


def test_cached_json_blob_rationale_is_repaired(tmp_path, monkeypatch):
    monkeypatch.setenv("FUSION_GENE_CURATION_CACHE_DIR", str(tmp_path))
    monkeypatch.setattr(
        gene_curation,
        "retrieve_pubmed_records_for_queries",
        lambda *args, **kwargs: [
            PubMedRecord(
                pmid="123",
                title="EML4 ALK fusion evidence",
                abstract="EML4 ALK fusions are oncogenic in lung cancer.",
            )
        ],
    )

    class FakeMessages:
        def create(self, **kwargs):
            return SimpleNamespace(
                content=[
                    SimpleNamespace(
                        type="text",
                        text=(
                            '{"fusion":"EML4::ALK","fusion_literature_identified":true,'
                            '"cancer_associated":true,'
                            '"rationale":"{\\"fusion\\":\\"EML4::ALK\\",'
                            '\\"rationale\\":\\"Recovered rationale.\\",'
                            '\\"supporting_pmids\\":[\\"123\\"]}",'
                            '"supporting_pmids":[],"retrieved_pmids":["123"],'
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

    result = gene_curation.curate_fusion("EML4::ALK", anthropic_api_key="configured")

    assert result["rationale"] == "Recovered rationale."
    assert result["supporting_pmids"] == ["123"]


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


def test_fusion_context_marks_gene_pair_only_when_annotation_unavailable():
    class Fusion:
        five_gene = "CD74"
        three_gene = "ROS1"

    contexts = gene_curation.fusion_contexts_by_gene(
        [Fusion()],
        annotation_results=[{"input": Fusion(), "error": "no breakpoint given"}],
    )

    ros1 = contexts["ROS1"][0]
    assert ros1.fusion == "CD74::ROS1"
    assert ros1.breakpoint_context_available is False
    assert ros1.fusion_specificity == "gene_pair_only"
    assert ros1.kinase_domain_status is None
    assert "Exact Genome Nexus breakpoint context was unavailable" in ros1.limitations[0]


def test_fusion_context_marks_gene_pair_only_when_breakpoint_unknown():
    class Fusion:
        five_gene = "CD74"
        three_gene = "ROS1"

    annotation_result = {
        "interface": {
            "frame_status": "unknown",
            "five_last_aa": None,
            "three_first_aa": None,
            "domains": [
                {
                    "gene": "ROS1",
                    "name": "Protein kinase domain",
                    "type": "domain",
                    "accession": "IPR000719",
                    "start": 1946,
                    "end": 2222,
                    "status": "UNKNOWN",
                }
            ],
        },
        "resolved": {
            "five": {"breakpoint": {"type": "unknown"}},
            "three": {"breakpoint": {"type": "unknown"}},
        },
    }

    contexts = gene_curation.fusion_contexts_by_gene(
        [Fusion()],
        annotation_results=[{"input": Fusion(), "result": annotation_result, "error": None}],
    )

    ros1 = contexts["ROS1"][0]
    assert ros1.breakpoint_context_available is False
    assert ros1.fusion_specificity == "gene_pair_only"
    assert ros1.kinase_domain_status == "unknown"
    assert "Exact Genome Nexus breakpoint context was unavailable" in ros1.limitations[0]


def test_curate_gene_prompt_warns_against_position_claims_for_gene_pair_only(monkeypatch):
    monkeypatch.setenv("FUSION_GENE_CURATION_CACHE", "0")
    monkeypatch.setattr(
        gene_curation,
        "retrieve_pubmed_records",
        lambda *args, **kwargs: [
            PubMedRecord(
                pmid="1",
                title="CD74 ROS1 fusion in lung cancer",
                abstract="CD74 ROS1 fusions are reported in lung cancer.",
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
                            '{"gene":"ROS1","cancer_associated":true,'
                            '"rationale":"ROS1 fusions are reported.",'
                            '"supporting_pmids":["1"],'
                            '"retrieved_pmids":["1"],'
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
        gene="ROS1",
        fusion="CD74::ROS1",
        side="three_prime",
        partner_gene="CD74",
        limitations=("Exact Genome Nexus breakpoint context was unavailable.",),
    )

    result = gene_curation.curate_gene(
        "ROS1",
        anthropic_api_key="configured",
        fusion_contexts=[context],
    )

    assert "Specificity: gene_pair_only" in seen["prompt"]
    assert "avoid claims about the exact exon" in seen["prompt"]
    assert result["fusion_contexts"][0]["breakpoint_context_available"] is False


def test_curate_fusion_genes_passes_token_controls(monkeypatch):
    class Fusion:
        five_gene = "EML4"
        three_gene = "ALK"

    seen = {}

    def fake_curate_gene(gene, **kwargs):
        seen[gene] = kwargs
        return {"gene": gene, "insufficient_evidence": True}

    def fake_curate_fusion(fusion, **kwargs):
        seen[fusion] = kwargs
        return {
            "fusion": fusion,
            "fusion_literature_identified": False,
            "insufficient_evidence": True,
        }

    monkeypatch.setenv("ANTHROPIC_API_KEY", "configured")
    monkeypatch.setenv("NCBI_API_KEY", "ncbi")
    monkeypatch.setenv("FUSION_GENE_CURATION_MODEL", "curation-model")
    monkeypatch.setenv("FUSION_GENE_CURATION_MAX_RESULTS", "3")
    monkeypatch.setenv("FUSION_GENE_CURATION_ABSTRACT_CHARS", "600")
    monkeypatch.setenv("FUSION_GENE_CURATION_WORKERS", "1")
    monkeypatch.setattr(gene_curation, "curate_fusion", fake_curate_fusion)
    monkeypatch.setattr(gene_curation, "curate_gene", fake_curate_gene)

    result = gene_curation.curate_fusion_genes([Fusion()])

    assert result["fusions"] == [
        {"fusion": "EML4::ALK", "fusion_literature_identified": False, "insufficient_evidence": True}
    ]
    assert result["genes"] == [
        {"gene": "ALK", "insufficient_evidence": True},
        {"gene": "EML4", "insufficient_evidence": True},
    ]
    assert seen["EML4::ALK"]["model"] == "curation-model"
    assert seen["ALK"]["model"] == "curation-model"
    assert seen["ALK"]["ncbi_api_key"] == "ncbi"
    assert seen["ALK"]["max_results"] == 3
    assert seen["ALK"]["abstract_chars"] == 600
    assert seen["ALK"]["fusion_contexts"][0].fusion == "EML4::ALK"


def test_curate_fusion_genes_skips_gene_calls_when_fusion_is_sufficient(monkeypatch):
    class Fusion:
        five_gene = "LMNA"
        three_gene = "NTRK1"

    def fake_curate_fusion(fusion, **kwargs):
        return {
            "fusion": fusion,
            "fusion_literature_identified": True,
            "cancer_associated": True,
            "rationale": "LMNA::NTRK1 is described in cancer literature.",
            "supporting_pmids": ["1"],
            "retrieved_pmids": ["1"],
            "insufficient_evidence": False,
        }

    def fake_curate_gene(gene, **kwargs):
        raise AssertionError(f"gene curation should be skipped for {gene}")

    monkeypatch.setenv("ANTHROPIC_API_KEY", "configured")
    monkeypatch.setenv("FUSION_GENE_CURATION_WORKERS", "1")
    monkeypatch.setattr(gene_curation, "curate_fusion", fake_curate_fusion)
    monkeypatch.setattr(gene_curation, "curate_gene", fake_curate_gene)

    result = gene_curation.curate_fusion_genes([Fusion()])

    assert result["fusions"][0]["fusion"] == "LMNA::NTRK1"
    assert result["genes"] == []


def test_curate_fusion_genes_force_gene_calls_when_requested(monkeypatch):
    class Fusion:
        five_gene = "LMNA"
        three_gene = "NTRK1"

    called_genes = []

    def fake_curate_fusion(fusion, **kwargs):
        return {
            "fusion": fusion,
            "fusion_literature_identified": True,
            "supporting_pmids": ["1"],
            "insufficient_evidence": False,
        }

    def fake_curate_gene(gene, **kwargs):
        called_genes.append(gene)
        return {"gene": gene, "insufficient_evidence": True}

    monkeypatch.setenv("ANTHROPIC_API_KEY", "configured")
    monkeypatch.setenv("FUSION_GENE_CURATION_WORKERS", "1")
    monkeypatch.setattr(gene_curation, "curate_fusion", fake_curate_fusion)
    monkeypatch.setattr(gene_curation, "curate_gene", fake_curate_gene)

    result = gene_curation.curate_fusion_genes([Fusion()], force_gene_curation=True)

    assert called_genes == ["LMNA", "NTRK1"]
    assert [item["gene"] for item in result["genes"]] == ["LMNA", "NTRK1"]
