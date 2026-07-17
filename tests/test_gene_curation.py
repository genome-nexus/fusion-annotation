import sys
from types import SimpleNamespace

from fusion_annotation import gene_curation
from fusion_annotation.gene_curation import PubMedRecord


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

    assert result == {
        "genes": [
            {"gene": "ALK", "insufficient_evidence": True},
            {"gene": "EML4", "insufficient_evidence": True},
        ]
    }
    assert seen["ALK"]["model"] == "curation-model"
    assert seen["ALK"]["ncbi_api_key"] == "ncbi"
    assert seen["ALK"]["max_results"] == 3
    assert seen["ALK"]["abstract_chars"] == 600
