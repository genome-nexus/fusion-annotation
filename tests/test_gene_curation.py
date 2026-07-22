import sys
from types import SimpleNamespace

from fusion_annotation import gene_curation
from fusion_annotation.gene_curation import (
    FusionCurationContext,
    PubMedRecord,
)


def test_retrieve_pubmed_records_preserves_mixed_content_xml(monkeypatch):
    class SearchResponse:
        text = ""

        def raise_for_status(self):
            return None

        def json(self):
            return {"esearchresult": {"idlist": ["123"]}}

    class FetchResponse:
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

    monkeypatch.setattr(gene_curation.requests, "get", fake_get)

    records = gene_curation.retrieve_pubmed_records("ALK")

    assert records == [
        PubMedRecord(
            pmid="123",
            title="ALK fusion evidence",
            abstract="Kinase domain retained. Functional cancer evidence.",
        )
    ]


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


def test_curate_gene_parses_json_wrapped_in_markdown_fence(monkeypatch):
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


def test_curate_gene_prompt_requests_concise_curator_rationale(monkeypatch):
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
