"""Server-side literature curation for Genome Nexus fusion batches."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
import logging
import os
from pathlib import Path
import tempfile
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from dataclasses import asdict
from time import time
from typing import Any
from typing import Iterable
from typing import Optional
from typing import Protocol

import requests


ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class PubMedRecord:
    pmid: str
    title: str
    abstract: str


@dataclass(frozen=True)
class FusionCurationContext:
    gene: str
    fusion: str
    side: str
    partner_gene: str
    five_transcript: Optional[str] = None
    three_transcript: Optional[str] = None
    five_exon: Optional[str] = None
    three_exon: Optional[str] = None
    five_genomic: Optional[str] = None
    three_genomic: Optional[str] = None
    five_protein_breakpoint: Optional[str] = None
    three_protein_breakpoint: Optional[str] = None
    retained_domains: tuple[str, ...] = ()
    lost_domains: tuple[str, ...] = ()
    disrupted_domains: tuple[str, ...] = ()
    kinase_gene: Optional[str] = None
    kinase_gene_side: Optional[str] = None
    kinase_domain_status: Optional[str] = None
    annotation_error: Optional[str] = None


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


def _safe_namespace(namespace: str) -> str:
    return "".join(
        char if char.isalnum() or char in {"-", "_"} else "_"
        for char in namespace
    )


class CacheBackend(Protocol):
    def read(
        self,
        namespace: str,
        key_payload: dict[str, Any],
        ttl_seconds: int = 0,
    ) -> Optional[Any]:
        """Return a cached value, or None on miss/stale/error."""

    def write(
        self,
        namespace: str,
        key_payload: dict[str, Any],
        value: Any,
        ttl_seconds: int = 0,
    ) -> None:
        """Persist a JSON-serializable cached value."""


class NullCacheBackend:
    def read(
        self,
        namespace: str,
        key_payload: dict[str, Any],
        ttl_seconds: int = 0,
    ) -> Optional[Any]:
        return None

    def write(
        self,
        namespace: str,
        key_payload: dict[str, Any],
        value: Any,
        ttl_seconds: int = 0,
    ) -> None:
        return None


class FileCacheBackend:
    def __init__(self, cache_dir: Path):
        self.cache_dir = cache_dir

    def _path(self, namespace: str, key_payload: dict[str, Any]) -> Path:
        return (
            self.cache_dir
            / _safe_namespace(namespace)
            / f"{_stable_cache_key(key_payload)}.json"
        )

    def read(
        self,
        namespace: str,
        key_payload: dict[str, Any],
        ttl_seconds: int = 0,
    ) -> Optional[Any]:
        path = self._path(namespace, key_payload)
        if not path.exists():
            return None
        if ttl_seconds > 0 and time() - path.stat().st_mtime > ttl_seconds:
            return None
        try:
            return json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError, TypeError):
            return None

    def write(
        self,
        namespace: str,
        key_payload: dict[str, Any],
        value: Any,
        ttl_seconds: int = 0,
    ) -> None:
        path = self._path(namespace, key_payload)
        tmp_path = path.with_name(f".{path.name}.{os.getpid()}.tmp")
        try:
            path.parent.mkdir(parents=True, exist_ok=True)
            tmp_path.write_text(
                json.dumps(value, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            os.replace(tmp_path, path)
        except OSError:
            try:
                tmp_path.unlink(missing_ok=True)
            except OSError:
                pass


class RedisCacheBackend:
    def __init__(self, redis_url: str, prefix: str):
        try:
            import redis
        except ImportError as exc:
            raise RuntimeError(
                "Redis cache backend requires the optional `redis` package."
            ) from exc
        self.client = redis.Redis.from_url(redis_url)
        self.prefix = prefix

    def _key(self, namespace: str, key_payload: dict[str, Any]) -> str:
        return f"{self.prefix}:{_safe_namespace(namespace)}:{_stable_cache_key(key_payload)}"

    def read(
        self,
        namespace: str,
        key_payload: dict[str, Any],
        ttl_seconds: int = 0,
    ) -> Optional[Any]:
        try:
            raw = self.client.get(self._key(namespace, key_payload))
        except Exception as exc:  # pragma: no cover - depends on Redis runtime
            logger.warning("Redis curation cache read failed: %s", exc)
            return None
        if not raw:
            return None
        try:
            cached = json.loads(raw)
        except (ValueError, TypeError):
            return None
        if (
            ttl_seconds > 0
            and time() - float(cached.get("stored_at", 0)) > ttl_seconds
        ):
            return None
        return cached.get("value")

    def write(
        self,
        namespace: str,
        key_payload: dict[str, Any],
        value: Any,
        ttl_seconds: int = 0,
    ) -> None:
        payload = {"stored_at": time(), "value": value}
        try:
            kwargs = {"ex": ttl_seconds} if ttl_seconds > 0 else {}
            self.client.set(
                self._key(namespace, key_payload),
                json.dumps(payload, sort_keys=True),
                **kwargs,
            )
        except Exception as exc:  # pragma: no cover - depends on Redis runtime
            logger.warning("Redis curation cache write failed: %s", exc)


_CACHE_BACKEND: Optional[CacheBackend] = None
_CACHE_BACKEND_CONFIG: Optional[tuple[str, str, str, str, str]] = None


def _cache_backend() -> CacheBackend:
    global _CACHE_BACKEND, _CACHE_BACKEND_CONFIG
    config = (
        os.environ.get("FUSION_GENE_CURATION_CACHE", "1"),
        os.environ.get("FUSION_GENE_CURATION_CACHE_BACKEND", "file"),
        os.environ.get("FUSION_GENE_CURATION_REDIS_URL")
        or os.environ.get("REDIS_URL")
        or "",
        os.environ.get(
            "FUSION_GENE_CURATION_CACHE_PREFIX",
            "fusion-annotation:gene-curation",
        ),
        os.environ.get("FUSION_GENE_CURATION_CACHE_DIR", ""),
    )
    if _CACHE_BACKEND is not None and _CACHE_BACKEND_CONFIG == config:
        return _CACHE_BACKEND

    if not _cache_enabled():
        _CACHE_BACKEND = NullCacheBackend()
        _CACHE_BACKEND_CONFIG = config
        return _CACHE_BACKEND
    backend = (
        os.environ.get("FUSION_GENE_CURATION_CACHE_BACKEND", "file")
        .strip()
        .lower()
    )
    if backend in {"0", "false", "no", "none", "disabled"}:
        _CACHE_BACKEND = NullCacheBackend()
        _CACHE_BACKEND_CONFIG = config
        return _CACHE_BACKEND
    if backend == "redis":
        redis_url = os.environ.get("FUSION_GENE_CURATION_REDIS_URL") or os.environ.get(
            "REDIS_URL"
        )
        if not redis_url:
            logger.warning(
                "FUSION_GENE_CURATION_CACHE_BACKEND=redis set without REDIS_URL; falling back to file cache."
            )
            _CACHE_BACKEND = FileCacheBackend(_cache_dir())
            _CACHE_BACKEND_CONFIG = config
            return _CACHE_BACKEND
        prefix = os.environ.get(
            "FUSION_GENE_CURATION_CACHE_PREFIX",
            "fusion-annotation:gene-curation",
        )
        try:
            _CACHE_BACKEND = RedisCacheBackend(redis_url, prefix)
        except RuntimeError as exc:
            logger.warning("%s Falling back to file cache.", exc)
            _CACHE_BACKEND = FileCacheBackend(_cache_dir())
        _CACHE_BACKEND_CONFIG = config
        return _CACHE_BACKEND
    if backend != "file":
        logger.warning(
            "Unknown curation cache backend %r; falling back to file cache.",
            backend,
        )
    _CACHE_BACKEND = FileCacheBackend(_cache_dir())
    _CACHE_BACKEND_CONFIG = config
    return _CACHE_BACKEND


def _read_cache(
    namespace: str,
    key_payload: dict[str, Any],
    ttl_seconds: int = 0,
) -> Optional[Any]:
    return _cache_backend().read(namespace, key_payload, ttl_seconds)


def _write_cache(
    namespace: str,
    key_payload: dict[str, Any],
    value: Any,
    ttl_seconds: int = 0,
) -> None:
    _cache_backend().write(namespace, key_payload, value, ttl_seconds)


def _cache_path(namespace: str, key_payload: dict[str, Any]) -> Path:
    return FileCacheBackend(_cache_dir())._path(namespace, key_payload)


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


def _fusion_label(fusion: object) -> str:
    return f"{getattr(fusion, 'five_gene')}::{getattr(fusion, 'three_gene')}"


def _optional_str(value: object) -> Optional[str]:
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def _input_exon(fusion: object, attr: str) -> Optional[str]:
    return _optional_str(getattr(fusion, attr, None))


def _breakpoint_genomic(partner: dict) -> Optional[str]:
    breakpoint = partner.get("breakpoint") or {}
    position = breakpoint.get("genomic_position")
    structure = partner.get("structure") or {}
    chrom = structure.get("chrom")
    if position is None or not chrom:
        return None
    return f"{chrom}:{position}"


def _resolved_exon(partner: dict, fallback: Optional[str]) -> Optional[str]:
    context = (partner.get("breakpoint") or {}).get("context") or {}
    return _optional_str(context.get("exon_rank")) or fallback


def _domain_label(domain: dict) -> str:
    name = str(domain.get("name") or domain.get("accession") or "domain")
    start = domain.get("start")
    end = domain.get("end")
    if start is not None and end is not None:
        return f"{name} ({start}-{end})"
    return name


def _domains_for_gene(result: Optional[dict], gene: str, status: str) -> tuple[str, ...]:
    if not result:
        return ()
    domains = result.get("interface", {}).get("domains") or []
    labels = []
    for domain in domains:
        if str(domain.get("gene", "")).upper() != gene.upper():
            continue
        if str(domain.get("status", "")).upper() != status:
            continue
        labels.append(_domain_label(domain))
    return tuple(dict.fromkeys(labels))


def _is_kinase_domain(domain: dict) -> bool:
    text = " ".join(
        str(domain.get(key, ""))
        for key in ("name", "type", "accession")
    ).lower()
    return "kinase" in text


def _kinase_signal(
    result: Optional[dict],
    five_gene: str,
    three_gene: str,
) -> tuple[Optional[str], Optional[str], Optional[str]]:
    if not result:
        return None, None, None
    domains = result.get("interface", {}).get("domains") or []
    kinase_domains = [domain for domain in domains if _is_kinase_domain(domain)]
    if not kinase_domains:
        return None, None, None

    status_rank = {"RETAINED": 3, "DISRUPTED": 2, "LOST": 1}
    kinase_domains.sort(
        key=lambda domain: status_rank.get(str(domain.get("status", "")).upper(), 0),
        reverse=True,
    )
    domain = kinase_domains[0]
    gene = str(domain.get("gene", "") or "")
    side = None
    if gene.upper() == str(five_gene).upper():
        side = "five_prime"
    elif gene.upper() == str(three_gene).upper():
        side = "three_prime"
    status = str(domain.get("status", "")).lower() or "unknown"
    if status not in {"retained", "lost", "disrupted"}:
        status = "unknown"
    return gene or None, side, status


def _context_for_gene(
    fusion: object,
    gene: str,
    *,
    result: Optional[dict] = None,
    error: Optional[str] = None,
) -> FusionCurationContext:
    five_gene = str(getattr(fusion, "five_gene"))
    three_gene = str(getattr(fusion, "three_gene"))
    side = "five_prime" if gene.upper() == five_gene.upper() else "three_prime"
    partner_gene = three_gene if side == "five_prime" else five_gene
    resolved = result.get("resolved", {}) if result else {}
    five = resolved.get("five") or {}
    three = resolved.get("three") or {}
    iface = result.get("interface", {}) if result else {}
    kinase_gene, kinase_side, kinase_status = _kinase_signal(result, five_gene, three_gene)
    return FusionCurationContext(
        gene=gene,
        fusion=_fusion_label(fusion),
        side=side,
        partner_gene=partner_gene,
        five_transcript=_optional_str(five.get("transcript")) or _optional_str(getattr(fusion, "five_transcript", None)),
        three_transcript=_optional_str(three.get("transcript")) or _optional_str(getattr(fusion, "three_transcript", None)),
        five_exon=_resolved_exon(five, _input_exon(fusion, "five_exon")),
        three_exon=_resolved_exon(three, _input_exon(fusion, "three_exon")),
        five_genomic=_breakpoint_genomic(five) or _optional_str(getattr(fusion, "five_genomic", None)),
        three_genomic=_breakpoint_genomic(three) or _optional_str(getattr(fusion, "three_genomic", None)),
        five_protein_breakpoint=f"p.{iface.get('five_last_aa')}" if iface.get("five_last_aa") is not None else None,
        three_protein_breakpoint=f"p.{iface.get('three_first_aa')}" if iface.get("three_first_aa") is not None else None,
        retained_domains=_domains_for_gene(result, gene, "RETAINED"),
        lost_domains=_domains_for_gene(result, gene, "LOST"),
        disrupted_domains=_domains_for_gene(result, gene, "DISRUPTED"),
        kinase_gene=kinase_gene,
        kinase_gene_side=kinase_side,
        kinase_domain_status=kinase_status,
        annotation_error=error,
    )


def fusion_contexts_by_gene(
    fusions: Iterable[object],
    annotation_results: Optional[Iterable[dict]] = None,
) -> dict[str, list[FusionCurationContext]]:
    result_items = list(annotation_results or [])
    contexts: dict[str, list[FusionCurationContext]] = {}
    for index, fusion in enumerate(fusions):
        item = result_items[index] if index < len(result_items) else {}
        result = item.get("result") if isinstance(item, dict) else None
        error = item.get("error") if isinstance(item, dict) else None
        for gene in (getattr(fusion, "five_gene"), getattr(fusion, "three_gene")):
            normalized = str(gene).strip().upper()
            if not normalized:
                continue
            contexts.setdefault(normalized, []).append(
                _context_for_gene(fusion, normalized, result=result, error=error)
            )
    return contexts


def _context_dicts(contexts: list[FusionCurationContext]) -> list[dict]:
    return [asdict(context) for context in contexts]


def _pubmed_queries(gene: str, fusion_contexts: list[FusionCurationContext]) -> list[str]:
    queries = [
        f'"{gene}"[Title/Abstract] AND (cancer OR tumor OR tumour OR carcinoma)',
    ]
    for context in fusion_contexts:
        fusion_hyphen = context.fusion.replace("::", "-")
        queries.extend([
            f'"{fusion_hyphen}" AND (cancer OR tumor OR tumour OR carcinoma)',
            f'"{gene}" AND "{context.partner_gene}" AND fusion',
        ])
        exon = context.five_exon if context.side == "five_prime" else context.three_exon
        if exon:
            queries.append(f'"{gene}" AND "exon {exon}" AND fusion')
        if context.kinase_gene and context.kinase_gene.upper() == gene.upper():
            queries.append(f'"{gene}" AND "kinase domain" AND fusion')
        if context.kinase_gene and context.kinase_domain_status == "retained":
            queries.append(f'"{context.kinase_gene}" AND "kinase domain retained"')
    return list(dict.fromkeys(queries))


def retrieve_pubmed_records(
    gene: str,
    *,
    ncbi_api_key: str = "",
    max_results: int = 8,
    fusion_contexts: Optional[list[FusionCurationContext]] = None,
) -> list[PubMedRecord]:
    fusion_contexts = fusion_contexts or []
    queries = _pubmed_queries(gene, fusion_contexts)
    cache_key = {
        "version": "pubmed-records-v1",
        "gene": gene.upper(),
        "max_results": max_results,
        "queries": queries,
    }
    ttl_seconds = max(0, int(os.environ.get("FUSION_GENE_CURATION_PUBMED_CACHE_TTL_SECONDS", "86400")))
    cached = _read_cache("pubmed-records", cache_key, ttl_seconds=ttl_seconds)
    if isinstance(cached, list):
        return [
            _record_from_cache(item)
            for item in cached
            if isinstance(item, dict) and item.get("pmid")
        ]

    pmids = []
    seen_pmids = set()
    for query in queries:
        params = {
            "db": "pubmed",
            "term": query,
            "retmode": "json",
            "retmax": str(max_results),
            "sort": "relevance",
        }
        if ncbi_api_key:
            params["api_key"] = ncbi_api_key
        search = requests.get(ESEARCH_URL, params=params, timeout=15)
        search.raise_for_status()
        for pmid in search.json().get("esearchresult", {}).get("idlist", []):
            if pmid not in seen_pmids:
                seen_pmids.add(pmid)
                pmids.append(pmid)
    if not pmids:
        _write_cache("pubmed-records", cache_key, [], ttl_seconds=ttl_seconds)
        return []

    fetch_params = {
        "db": "pubmed",
        "id": ",".join(pmids[: max_results * len(queries)]),
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
    _write_cache(
        "pubmed-records",
        cache_key,
        [_record_to_cache(record) for record in records],
        ttl_seconds=ttl_seconds,
    )
    return records


def _curation_cache_key(
    gene: str,
    records: list[PubMedRecord],
    model: str,
    fusion_contexts: list[FusionCurationContext],
) -> dict[str, Any]:
    return {
        "version": "gene-curation-v1",
        "gene": gene.upper(),
        "model": model,
        "fusion_contexts": _context_dicts(fusion_contexts),
        "records": [
            {
                "pmid": record.pmid,
                "title_sha256": _text_digest(record.title),
                "abstract_sha256": _text_digest(record.abstract),
            }
            for record in records
        ],
    }


def _no_pubmed_evidence_result(
    gene: str,
    fusion_contexts: Optional[list[FusionCurationContext]] = None,
) -> dict:
    return {
        "gene": gene,
        "cancer_associated": None,
        "rationale": "No PubMed abstracts were retrieved for this gene.",
        "supporting_pmids": [],
        "retrieved_pmids": [],
        "fusion_contexts": _context_dicts(fusion_contexts or []),
        "insufficient_evidence": True,
    }


def _format_domain_list(values: tuple[str, ...]) -> str:
    return ", ".join(values[:6]) if values else "none"


def _format_fusion_contexts_for_prompt(contexts: list[FusionCurationContext]) -> str:
    if not contexts:
        return "No Genome Nexus fusion-position context was available."
    lines = []
    for context in contexts:
        lines.extend([
            f"- Fusion: {context.fusion}",
            f"  Gene side: {context.side}; partner: {context.partner_gene}",
            (
                "  5' partner: "
                f"transcript={context.five_transcript or 'unknown'}; "
                f"exon={context.five_exon or 'unknown'}; "
                f"genomic={context.five_genomic or 'unknown'}; "
                f"protein={context.five_protein_breakpoint or 'unknown'}"
            ),
            (
                "  3' partner: "
                f"transcript={context.three_transcript or 'unknown'}; "
                f"exon={context.three_exon or 'unknown'}; "
                f"genomic={context.three_genomic or 'unknown'}; "
                f"protein={context.three_protein_breakpoint or 'unknown'}"
            ),
            (
                "  Domains for this gene: "
                f"retained={_format_domain_list(context.retained_domains)}; "
                f"lost={_format_domain_list(context.lost_domains)}; "
                f"disrupted={_format_domain_list(context.disrupted_domains)}"
            ),
            (
                "  Kinase context: "
                f"gene={context.kinase_gene or 'unknown'}; "
                f"side={context.kinase_gene_side or 'unknown'}; "
                f"domain_status={context.kinase_domain_status or 'unknown'}"
            ),
        ])
        if context.annotation_error:
            lines.append(f"  Fusion annotation limitation: {context.annotation_error}")
    return "\n".join(lines)


def curate_gene(
    gene: str,
    *,
    anthropic_api_key: str,
    ncbi_api_key: str = "",
    model: str = "claude-3-5-haiku-latest",
    max_results: int = 8,
    abstract_chars: int = 1200,
    fusion_contexts: Optional[list[FusionCurationContext]] = None,
) -> dict:
    if not anthropic_api_key:
        raise GeneCurationUnavailable("ANTHROPIC_API_KEY is not configured for server-side curation.")

    fusion_contexts = fusion_contexts or []
    records = retrieve_pubmed_records(
        gene,
        ncbi_api_key=ncbi_api_key,
        max_results=max_results,
        fusion_contexts=fusion_contexts,
    )
    if not records:
        return _no_pubmed_evidence_result(gene, fusion_contexts)

    cache_key = _curation_cache_key(gene, records, model, fusion_contexts)
    cache_key["abstract_chars"] = abstract_chars
    result_ttl_seconds = max(
        0,
        int(os.environ.get("FUSION_GENE_CURATION_RESULT_CACHE_TTL_SECONDS", "2592000")),
    )
    cached = _read_cache("gene-curation", cache_key, ttl_seconds=result_ttl_seconds)
    if isinstance(cached, dict):
        return cached

    context = "\n\n".join(
        f"PMID {record.pmid}\nTitle: {record.title}\nAbstract: {record.abstract[:abstract_chars]}"
        for record in records
    )
    prompt = f"""\
Gene: {gene}

Genome Nexus fusion-position context:
{_format_fusion_contexts_for_prompt(fusion_contexts)}

PubMed context:
{context}

Return one JSON object with:
- gene
- cancer_associated: true/false/null
- rationale: concise curator-facing scan text grounded only in the PubMed context.
  Write 1-2 short sentences, ideally 40-75 words total. Prioritize the
  classification-relevant conclusion, strongest mechanism/cancer context, and
  one caveat if needed. Do not enumerate every paper or make a literature-review paragraph.
- supporting_pmids: up to 4 PMIDs from the context
- retrieved_pmids: all PMIDs provided in the context
- fusion_contexts: echo the provided fusion context as compact JSON-compatible objects
- insufficient_evidence: true when the context is too sparse

Do not infer transcript, breakpoint, domain-retention, or kinase-domain status
beyond the Genome Nexus context. Do not mark a result cancer-associated solely
because a kinase domain is retained; supporting PubMed evidence is still required.
"""

    import anthropic

    client = anthropic.Anthropic(api_key=anthropic_api_key)
    response = client.messages.create(
        model=model,
        max_tokens=1024,
        system=(
            "You are a cancer genomics literature curator. "
            "Use only the provided PubMed context. Return valid JSON only. "
            "Keep rationale text concise and optimized for fast curator review."
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
            "fusion_contexts": _context_dicts(fusion_contexts),
            "insufficient_evidence": True,
        }
    payload.setdefault("gene", gene)
    payload.setdefault("retrieved_pmids", [record.pmid for record in records])
    payload["fusion_contexts"] = _context_dicts(fusion_contexts)
    _write_cache("gene-curation", cache_key, payload, ttl_seconds=result_ttl_seconds)
    return payload


def curate_fusion_genes(
    fusions: Iterable[object],
    annotation_results: Optional[Iterable[dict]] = None,
) -> dict:
    anthropic_api_key = os.environ.get("ANTHROPIC_API_KEY", "")
    if not anthropic_api_key:
        raise GeneCurationUnavailable("ANTHROPIC_API_KEY is not configured for server-side curation.")

    fusions = list(fusions)
    contexts_by_gene = fusion_contexts_by_gene(fusions, annotation_results)
    genes = list(contexts_by_gene) or unique_genes_from_fusions(fusions)
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
                fusion_contexts=contexts_by_gene.get(gene, []),
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
                    "fusion_contexts": _context_dicts(contexts_by_gene.get(gene, [])),
                })

    results.sort(key=lambda item: item.get("gene", ""))
    return {"genes": results}
