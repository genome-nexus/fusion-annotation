"""Server-side literature curation for genes found in fusion batches."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import os
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from typing import Iterable

import requests


ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


@dataclass(frozen=True)
class PubMedRecord:
    pmid: str
    title: str
    abstract: str


class GeneCurationUnavailable(RuntimeError):
    """Raised when server-side curation is not configured."""


def unique_genes_from_fusions(fusions: Iterable[object]) -> list[str]:
    genes: list[str] = []
    seen = set()
    for fusion in fusions:
        for gene in (getattr(fusion, "five_gene"), getattr(fusion, "three_gene")):
            normalized = str(gene).strip().upper()
            if normalized and normalized not in seen:
                seen.add(normalized)
                genes.append(normalized)
    return genes


def retrieve_pubmed_records(
    gene: str,
    *,
    ncbi_api_key: str = "",
    max_results: int = 8,
) -> list[PubMedRecord]:
    params = {
        "db": "pubmed",
        "term": f'"{gene}"[Title/Abstract] AND (cancer OR tumor OR tumour OR carcinoma)',
        "retmode": "json",
        "retmax": str(max_results),
        "sort": "relevance",
    }
    if ncbi_api_key:
        params["api_key"] = ncbi_api_key
    search = requests.get(ESEARCH_URL, params=params, timeout=15)
    search.raise_for_status()
    pmids = search.json().get("esearchresult", {}).get("idlist", [])
    if not pmids:
        return []

    fetch_params = {
        "db": "pubmed",
        "id": ",".join(pmids),
        "rettype": "abstract",
        "retmode": "xml",
    }
    if ncbi_api_key:
        fetch_params["api_key"] = ncbi_api_key
    fetch = requests.get(EFETCH_URL, params=fetch_params, timeout=30)
    fetch.raise_for_status()

    records: list[PubMedRecord] = []
    root = ET.fromstring(fetch.text)
    for article in root.findall(".//PubmedArticle"):
        pmid_el = article.find(".//PMID")
        pmid = pmid_el.text if pmid_el is not None else None
        title_el = article.find(".//ArticleTitle")
        title = (title_el.text or "").strip() if title_el is not None else ""
        abstract_parts = article.findall(".//AbstractText")
        abstract = " ".join(
            (part.text or "").strip() for part in abstract_parts if part.text
        ).strip()
        if pmid and abstract:
            records.append(PubMedRecord(pmid=pmid, title=title, abstract=abstract))
    return records


def curate_gene(
    gene: str,
    *,
    anthropic_api_key: str,
    ncbi_api_key: str = "",
    model: str = "claude-3-5-haiku-latest",
) -> dict:
    if not anthropic_api_key:
        raise GeneCurationUnavailable("ANTHROPIC_API_KEY is not configured for server-side curation.")

    records = retrieve_pubmed_records(gene, ncbi_api_key=ncbi_api_key)
    context = "\n\n".join(
        f"PMID {record.pmid}\nTitle: {record.title}\nAbstract: {record.abstract[:1200]}"
        for record in records
    ) or "No PubMed abstracts retrieved."
    prompt = f"""\
Gene: {gene}

PubMed context:
{context}

Return one JSON object with:
- gene
- cancer_associated: true/false/null
- rationale: concise curator-facing explanation grounded only in the PubMed context
- supporting_pmids: up to 4 PMIDs from the context
- retrieved_pmids: all PMIDs provided in the context
- insufficient_evidence: true when the context is too sparse
"""

    import anthropic

    client = anthropic.Anthropic(api_key=anthropic_api_key)
    response = client.messages.create(
        model=model,
        max_tokens=1024,
        system=(
            "You are a cancer genomics literature curator. "
            "Use only the provided PubMed context. Return valid JSON only."
        ),
        messages=[{"role": "user", "content": prompt}],
    )
    text = "".join(
        block.text for block in response.content
        if getattr(block, "type", None) == "text"
    ).strip()
    try:
        payload = json.loads(text)
    except json.JSONDecodeError:
        payload = {
            "gene": gene,
            "cancer_associated": None,
            "rationale": text,
            "supporting_pmids": [],
            "retrieved_pmids": [record.pmid for record in records],
            "insufficient_evidence": True,
        }
    payload.setdefault("gene", gene)
    payload.setdefault("retrieved_pmids", [record.pmid for record in records])
    return payload


def curate_fusion_genes(fusions: Iterable[object]) -> dict:
    anthropic_api_key = os.environ.get("ANTHROPIC_API_KEY", "")
    if not anthropic_api_key:
        raise GeneCurationUnavailable("ANTHROPIC_API_KEY is not configured for server-side curation.")

    genes = unique_genes_from_fusions(fusions)
    ncbi_api_key = os.environ.get("NCBI_API_KEY", "")
    model = os.environ.get("FUSION_GENE_CURATION_MODEL", "claude-3-5-haiku-latest")
    max_workers = max(1, int(os.environ.get("FUSION_GENE_CURATION_WORKERS", "3")))
    results = []

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                curate_gene,
                gene,
                anthropic_api_key=anthropic_api_key,
                ncbi_api_key=ncbi_api_key,
                model=model,
            ): gene
            for gene in genes
        }
        for future in as_completed(futures):
            gene = futures[future]
            try:
                results.append(future.result())
            except Exception as exc:
                results.append({
                    "gene": gene,
                    "error": str(exc),
                    "insufficient_evidence": True,
                    "supporting_pmids": [],
                    "retrieved_pmids": [],
                })

    results.sort(key=lambda item: item.get("gene", ""))
    return {"genes": results}
