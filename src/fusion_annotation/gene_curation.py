"""Server-side literature curation for genes found in fusion batches."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
import os
from pathlib import Path
import tempfile
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from time import time
from typing import Any
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


def _cache_enabled() -> bool:
    return os.environ.get("FUSION_GENE_CURATION_CACHE", "1").lower() not in {"0", "false", "no"}


def _cache_dir() -> Path:
    override = os.environ.get("FUSION_GENE_CURATION_CACHE_DIR")
    if override:
        return Path(override)
    return Path(tempfile.gettempdir()) / "fusion-annotation" / "gene-curation-cache"


def _stable_cache_key(payload: dict[str, Any]) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _text_digest(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _cache_path(namespace: str, key_payload: dict[str, Any]) -> Path:
    safe_namespace = "".join(
        char if char.isalnum() or char in {"-", "_"} else "_"
        for char in namespace
    )
    return _cache_dir() / safe_namespace / f"{_stable_cache_key(key_payload)}.json"


def _read_cache(namespace: str, key_payload: dict[str, Any], ttl_seconds: int = 0) -> Any | None:
    if not _cache_enabled():
        return None
    path = _cache_path(namespace, key_payload)
    if not path.exists():
        return None
    if ttl_seconds > 0 and time() - path.stat().st_mtime > ttl_seconds:
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError, TypeError):
        return None


def _write_cache(namespace: str, key_payload: dict[str, Any], value: Any) -> None:
    if not _cache_enabled():
        return
    path = _cache_path(namespace, key_payload)
    tmp_path = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp_path.write_text(json.dumps(value, sort_keys=True) + "\n", encoding="utf-8")
        os.replace(tmp_path, path)
    except OSError:
        try:
            tmp_path.unlink(missing_ok=True)
        except OSError:
            pass


def _record_to_cache(record: PubMedRecord) -> dict[str, str]:
    return {"pmid": record.pmid, "title": record.title, "abstract": record.abstract}


def _record_from_cache(value: dict[str, str]) -> PubMedRecord:
    return PubMedRecord(
        pmid=str(value.get("pmid", "")),
        title=str(value.get("title", "")),
        abstract=str(value.get("abstract", "")),
    )


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
    cache_key = {
        "version": "pubmed-records-v1",
        "gene": gene.upper(),
        "max_results": max_results,
    }
    ttl_seconds = max(0, int(os.environ.get("FUSION_GENE_CURATION_PUBMED_CACHE_TTL_SECONDS", "86400")))
    cached = _read_cache("pubmed-records", cache_key, ttl_seconds=ttl_seconds)
    if isinstance(cached, list):
        return [
            _record_from_cache(item)
            for item in cached
            if isinstance(item, dict) and item.get("pmid")
        ]

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
    _write_cache("pubmed-records", cache_key, [_record_to_cache(record) for record in records])
    return records


def _curation_cache_key(gene: str, records: list[PubMedRecord], model: str) -> dict[str, Any]:
    return {
        "version": "gene-curation-v1",
        "gene": gene.upper(),
        "model": model,
        "records": [
            {
                "pmid": record.pmid,
                "title_sha256": _text_digest(record.title),
                "abstract_sha256": _text_digest(record.abstract),
            }
            for record in records
        ],
    }


def _no_pubmed_evidence_result(gene: str) -> dict:
    return {
        "gene": gene,
        "cancer_associated": None,
        "rationale": "No PubMed abstracts were retrieved for this gene.",
        "supporting_pmids": [],
        "retrieved_pmids": [],
        "insufficient_evidence": True,
    }


def curate_gene(
    gene: str,
    *,
    anthropic_api_key: str,
    ncbi_api_key: str = "",
    model: str = "claude-3-5-haiku-latest",
    max_results: int = 8,
    abstract_chars: int = 1200,
) -> dict:
    if not anthropic_api_key:
        raise GeneCurationUnavailable("ANTHROPIC_API_KEY is not configured for server-side curation.")

    records = retrieve_pubmed_records(gene, ncbi_api_key=ncbi_api_key, max_results=max_results)
    if not records:
        return _no_pubmed_evidence_result(gene)

    cache_key = _curation_cache_key(gene, records, model)
    cache_key["abstract_chars"] = abstract_chars
    cached = _read_cache("gene-curation", cache_key)
    if isinstance(cached, dict):
        return cached

    context = "\n\n".join(
        f"PMID {record.pmid}\nTitle: {record.title}\nAbstract: {record.abstract[:abstract_chars]}"
        for record in records
    )
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
    _write_cache("gene-curation", cache_key, payload)
    return payload


def curate_fusion_genes(fusions: Iterable[object]) -> dict:
    anthropic_api_key = os.environ.get("ANTHROPIC_API_KEY", "")
    if not anthropic_api_key:
        raise GeneCurationUnavailable("ANTHROPIC_API_KEY is not configured for server-side curation.")

    genes = unique_genes_from_fusions(fusions)
    ncbi_api_key = os.environ.get("NCBI_API_KEY", "")
    model = os.environ.get("FUSION_GENE_CURATION_MODEL", "claude-3-5-haiku-latest")
    max_results = max(1, int(os.environ.get("FUSION_GENE_CURATION_MAX_RESULTS", "8")))
    abstract_chars = max(200, int(os.environ.get("FUSION_GENE_CURATION_ABSTRACT_CHARS", "1200")))
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
                max_results=max_results,
                abstract_chars=abstract_chars,
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
