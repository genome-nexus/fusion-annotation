"""Server-side literature curation for Genome Nexus fusion batches."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
import logging
import os
from pathlib import Path
import tempfile
import threading
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from dataclasses import asdict
from dataclasses import field
from time import monotonic
from time import sleep
from time import time
from typing import Any
from typing import Iterable
from typing import Optional
from typing import Protocol

import requests


ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
DEFAULT_ONCOKB_API_BASE_URL = "https://www.oncokb.org/api/v1"
DEFAULT_CURATION_MODEL = "claude-haiku-4-5-20251001"
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
    breakpoint_context_available: bool = False
    fusion_specificity: str = "gene_pair_only"
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
    limitations: tuple[str, ...] = ()
    annotation_error: Optional[str] = None


class GeneCurationUnavailable(RuntimeError):
    """Raised when server-side curation is not configured."""


class PubMedRateLimitError(RuntimeError):
    """Raised when NCBI keeps rejecting PubMed retrieval for rate limiting."""


class OncoKBGeneLookupError(RuntimeError):
    """Raised when OncoKB curated-gene lookup fails."""


@dataclass
class NcbiClient:
    api_key: str = ""
    min_interval_seconds: Optional[float] = None
    max_retries: int = 3
    backoff_seconds: float = 1.0
    _lock: threading.Lock = field(default_factory=threading.Lock, init=False)
    _next_request_at: float = field(default=0.0, init=False)

    def __post_init__(self) -> None:
        if self.min_interval_seconds is None:
            configured = os.environ.get("FUSION_GENE_CURATION_NCBI_MIN_INTERVAL_SECONDS")
            if configured:
                self.min_interval_seconds = max(0.0, float(configured))
            else:
                self.min_interval_seconds = 0.12 if self.api_key else 0.4

    def get(self, url: str, *, params: dict[str, str], timeout: int) -> requests.Response:
        if self.api_key and "api_key" not in params:
            params = {**params, "api_key": self.api_key}

        response = None
        for attempt in range(self.max_retries + 1):
            self._wait_for_slot()
            response = requests.get(url, params=params, timeout=timeout)
            if response.status_code not in {429, 500, 502, 503, 504}:
                response.raise_for_status()
                return response
            if attempt >= self.max_retries:
                break
            retry_after = response.headers.get("Retry-After")
            try:
                delay = float(retry_after) if retry_after else self.backoff_seconds * (2 ** attempt)
            except ValueError:
                delay = self.backoff_seconds * (2 ** attempt)
            sleep(max(0.0, delay))

        if response is not None and response.status_code == 429:
            raise PubMedRateLimitError(
                "NCBI PubMed rate limit was reached while retrieving literature. "
                "Add NCBI_API_KEY or retry after the rate-limit window resets."
            )
        if response is not None:
            response.raise_for_status()
        raise PubMedRateLimitError("NCBI PubMed request failed before receiving a response.")

    def _wait_for_slot(self) -> None:
        with self._lock:
            now = monotonic()
            wait_seconds = max(0.0, self._next_request_at - now)
            if wait_seconds:
                sleep(wait_seconds)
                now = monotonic()
            self._next_request_at = now + (self.min_interval_seconds or 0.0)


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


def _fusion_value(fusion: object, attr: str) -> Any:
    if isinstance(fusion, dict):
        return fusion.get(attr)
    return getattr(fusion, attr, None)


def unique_genes_from_fusions(fusions: Iterable[object]) -> list[str]:
    genes: list[str] = []
    seen = set()
    for fusion in fusions:
        if fusion is None:
            continue
        for gene in (_fusion_value(fusion, "five_gene"), _fusion_value(fusion, "three_gene")):
            if gene is None:
                continue
            normalized = str(gene).strip().upper()
            if normalized and normalized not in seen:
                seen.add(normalized)
                genes.append(normalized)
    return genes


def _fusion_label(fusion: object) -> str:
    return f"{_fusion_value(fusion, 'five_gene')}::{_fusion_value(fusion, 'three_gene')}"


def _optional_str(value: object) -> Optional[str]:
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def _input_exon(fusion: object, attr: str) -> Optional[str]:
    return _optional_str(_fusion_value(fusion, attr))


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


def _fusion_specificity(result: Optional[dict]) -> str:
    if not result:
        return "gene_pair_only"

    iface = result.get("interface", {})
    if iface.get("frame_status") == "unknown":
        return "gene_pair_only"
    if iface.get("five_last_aa") is not None or iface.get("three_first_aa") is not None:
        return "protein_domain_level"
    domains = iface.get("domains") or []
    if any(str(domain.get("status", "")).upper() != "UNKNOWN" for domain in domains):
        return "protein_domain_level"

    resolved = result.get("resolved", {})
    five_breakpoint = (resolved.get("five") or {}).get("breakpoint") or {}
    three_breakpoint = (resolved.get("three") or {}).get("breakpoint") or {}
    if five_breakpoint.get("type") != "unknown" and three_breakpoint.get("type") != "unknown":
        return "exon_level"

    return "gene_pair_only"


def _context_limitations(
    *,
    result: Optional[dict],
    error: Optional[str],
) -> tuple[str, ...]:
    if result and _fusion_specificity(result) != "gene_pair_only":
        return ()

    limitations = [
        (
            "Exact Genome Nexus breakpoint context was unavailable, so this "
            "curation is based on the reported fusion gene pair and PubMed "
            "literature rather than the precise exon/protein breakpoint."
        )
    ]
    if error:
        limitations.append(f"Genome Nexus annotation did not complete: {error}")
    return tuple(limitations)


def _context_for_gene(
    fusion: object,
    gene: str,
    *,
    result: Optional[dict] = None,
    error: Optional[str] = None,
) -> FusionCurationContext:
    five_gene = str(_fusion_value(fusion, "five_gene"))
    three_gene = str(_fusion_value(fusion, "three_gene"))
    side = "five_prime" if gene.upper() == five_gene.upper() else "three_prime"
    partner_gene = three_gene if side == "five_prime" else five_gene
    resolved = result.get("resolved", {}) if result else {}
    five = resolved.get("five") or {}
    three = resolved.get("three") or {}
    iface = result.get("interface", {}) if result else {}
    kinase_gene, kinase_side, kinase_status = _kinase_signal(result, five_gene, three_gene)
    fusion_specificity = _fusion_specificity(result)
    breakpoint_context_available = fusion_specificity != "gene_pair_only"
    return FusionCurationContext(
        gene=gene,
        fusion=_fusion_label(fusion),
        side=side,
        partner_gene=partner_gene,
        breakpoint_context_available=breakpoint_context_available,
        fusion_specificity=fusion_specificity,
        five_transcript=_optional_str(five.get("transcript")) or _optional_str(_fusion_value(fusion, "five_transcript")),
        three_transcript=_optional_str(three.get("transcript")) or _optional_str(_fusion_value(fusion, "three_transcript")),
        five_exon=_resolved_exon(five, _input_exon(fusion, "five_exon")),
        three_exon=_resolved_exon(three, _input_exon(fusion, "three_exon")),
        five_genomic=_breakpoint_genomic(five) or _optional_str(_fusion_value(fusion, "five_genomic")),
        three_genomic=_breakpoint_genomic(three) or _optional_str(_fusion_value(fusion, "three_genomic")),
        five_protein_breakpoint=f"p.{iface.get('five_last_aa')}" if iface.get("five_last_aa") is not None else None,
        three_protein_breakpoint=f"p.{iface.get('three_first_aa')}" if iface.get("three_first_aa") is not None else None,
        retained_domains=_domains_for_gene(result, gene, "RETAINED"),
        lost_domains=_domains_for_gene(result, gene, "LOST"),
        disrupted_domains=_domains_for_gene(result, gene, "DISRUPTED"),
        kinase_gene=kinase_gene,
        kinase_gene_side=kinase_side,
        kinase_domain_status=kinase_status,
        limitations=_context_limitations(result=result, error=error),
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
        for gene in (_fusion_value(fusion, "five_gene"), _fusion_value(fusion, "three_gene")):
            if gene is None:
                continue
            normalized = str(gene).strip().upper()
            if not normalized:
                continue
            contexts.setdefault(normalized, []).append(
                _context_for_gene(fusion, normalized, result=result, error=error)
            )
    return contexts


def fusion_contexts_by_fusion(
    fusions: Iterable[object],
    annotation_results: Optional[Iterable[dict]] = None,
) -> dict[str, list[FusionCurationContext]]:
    result_items = list(annotation_results or [])
    contexts: dict[str, list[FusionCurationContext]] = {}
    for index, fusion in enumerate(fusions):
        item = result_items[index] if index < len(result_items) else {}
        result = item.get("result") if isinstance(item, dict) else None
        error = item.get("error") if isinstance(item, dict) else None
        label = _fusion_label(fusion)
        for gene in (_fusion_value(fusion, "five_gene"), _fusion_value(fusion, "three_gene")):
            if gene is None:
                continue
            normalized = str(gene).strip().upper()
            if not normalized:
                continue
            contexts.setdefault(label, []).append(
                _context_for_gene(fusion, normalized, result=result, error=error)
            )
    return contexts


def _context_dicts(contexts: list[FusionCurationContext]) -> list[dict]:
    return [asdict(context) for context in contexts]


def _strip_markdown_json_fence(text: str) -> str:
    cleaned = text.strip()
    if not cleaned.startswith("```"):
        return cleaned
    lines = cleaned.splitlines()
    if len(lines) >= 2 and lines[-1].strip().startswith("```"):
        return "\n".join(lines[1:-1]).strip()
    if len(lines) >= 2 and lines[0].strip().startswith("```"):
        return "\n".join(lines[1:]).strip()
    return cleaned


def _extract_json_object(text: str) -> dict:
    cleaned = _strip_markdown_json_fence(text)
    try:
        parsed = json.loads(cleaned)
    except json.JSONDecodeError:
        parsed = None
    if isinstance(parsed, dict):
        return parsed

    start = cleaned.find("{")
    if start == -1:
        raise json.JSONDecodeError("no JSON object found", cleaned, 0)

    in_string = False
    escaped = False
    depth = 0
    for index, char in enumerate(cleaned[start:], start=start):
        if escaped:
            escaped = False
            continue
        if char == "\\":
            escaped = True
            continue
        if char == '"':
            in_string = not in_string
            continue
        if in_string:
            continue
        if char == "{":
            depth += 1
        elif char == "}":
            depth -= 1
            if depth == 0:
                candidate = cleaned[start:index + 1]
                parsed = json.loads(candidate)
                if isinstance(parsed, dict):
                    return parsed
                break
    raise json.JSONDecodeError("no complete JSON object found", cleaned, start)


def _looks_like_json_blob(value: object) -> bool:
    text = str(value or "").strip()
    if not text:
        return False
    if text.startswith("{") or text.startswith("```"):
        return True
    return (
        '"rationale"' in text
        and (
            '"supporting_pmids"' in text
            or '"retrieved_pmids"' in text
            or '"fusion_literature_identified"' in text
        )
    )


def _jsonish_string_field(text: str, key: str) -> Optional[str]:
    marker = f'"{key}"'
    start = text.find(marker)
    if start == -1:
        return None
    colon = text.find(":", start + len(marker))
    if colon == -1:
        return None
    quote = text.find('"', colon + 1)
    if quote == -1:
        return None

    escaped = False
    for index in range(quote + 1, len(text)):
        char = text[index]
        if escaped:
            escaped = False
            continue
        if char == "\\":
            escaped = True
            continue
        if char == '"':
            try:
                value = json.loads(text[quote:index + 1])
            except json.JSONDecodeError:
                return None
            return value if isinstance(value, str) else None
    return None


def _jsonish_array_field(text: str, key: str) -> Optional[list]:
    marker = f'"{key}"'
    start = text.find(marker)
    if start == -1:
        return None
    colon = text.find(":", start + len(marker))
    if colon == -1:
        return None
    bracket = text.find("[", colon + 1)
    if bracket == -1:
        return None

    in_string = False
    escaped = False
    depth = 0
    for index, char in enumerate(text[bracket:], start=bracket):
        if escaped:
            escaped = False
            continue
        if char == "\\":
            escaped = True
            continue
        if char == '"':
            in_string = not in_string
            continue
        if in_string:
            continue
        if char == "[":
            depth += 1
        elif char == "]":
            depth -= 1
            if depth == 0:
                try:
                    value = json.loads(text[bracket:index + 1])
                except json.JSONDecodeError:
                    return None
                return value if isinstance(value, list) else None
    return None


def _jsonish_bool_or_null_field(text: str, key: str) -> Optional[bool]:
    marker = f'"{key}"'
    start = text.find(marker)
    if start == -1:
        return None
    colon = text.find(":", start + len(marker))
    if colon == -1:
        return None
    value = text[colon + 1:].lstrip()
    if value.startswith("true"):
        return True
    if value.startswith("false"):
        return False
    return None


def _repair_payload_from_jsonish_text(payload: dict, text: str) -> dict:
    repaired = dict(payload)
    cleaned = _strip_markdown_json_fence(text)
    for key in ("gene", "fusion", "rationale"):
        value = _jsonish_string_field(cleaned, key)
        if value:
            repaired[key] = value
    for key in ("supporting_pmids", "retrieved_pmids"):
        value = _jsonish_array_field(cleaned, key)
        if value is not None:
            repaired[key] = value
    for key in ("fusion_literature_identified", "cancer_associated", "insufficient_evidence"):
        value = _jsonish_bool_or_null_field(cleaned, key)
        if value is not None:
            repaired[key] = value
    return repaired


def _normal_model_parse_fallback(
    *,
    subject_label: str,
    subject: str,
    records: list[PubMedRecord],
    fusion_contexts: list[FusionCurationContext],
    fusion_result: bool,
    raw_text: str = "",
) -> dict:
    payload = {
        subject_label: subject,
        "cancer_associated": None,
        "rationale": (
            "The curation model returned a malformed structured response, so this "
            "result should be regenerated before review."
        ),
        "supporting_pmids": [],
        "retrieved_pmids": [record.pmid for record in records],
        "fusion_contexts": _context_dicts(fusion_contexts),
        "insufficient_evidence": True,
    }
    if fusion_result:
        payload["fusion_literature_identified"] = None
    if raw_text:
        payload = _repair_payload_from_jsonish_text(payload, raw_text)
    return payload


def _repair_json_rationale(payload: dict) -> dict:
    rationale = payload.get("rationale")
    if not _looks_like_json_blob(rationale):
        return payload
    try:
        parsed = _extract_json_object(str(rationale))
    except json.JSONDecodeError:
        repaired = _repair_payload_from_jsonish_text(payload, str(rationale))
        if repaired.get("rationale") != rationale and repaired.get("rationale"):
            return repaired
        payload["rationale"] = ""
        payload["insufficient_evidence"] = True
        return payload
    repaired = {**payload, **parsed}
    if _looks_like_json_blob(repaired.get("rationale")):
        repaired["rationale"] = (
            "The curation model returned malformed JSON in the rationale field; "
            "rerun curation to regenerate this result."
        )
        repaired["insufficient_evidence"] = True
    return repaired


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


def _fusion_pubmed_queries(fusion: str, fusion_contexts: list[FusionCurationContext]) -> list[str]:
    fusion_hyphen = fusion.replace("::", "-")
    genes = [part.strip() for part in fusion.replace("--", "::").split("::") if part.strip()]
    five_gene = fusion_contexts[0].gene if fusion_contexts and fusion_contexts[0].side == "five_prime" else None
    three_gene = fusion_contexts[0].partner_gene if five_gene else None
    if not five_gene and len(genes) == 2:
        five_gene, three_gene = genes
    if not three_gene and fusion_contexts:
        first = fusion_contexts[0]
        five_gene = first.gene if first.side == "five_prime" else first.partner_gene
        three_gene = first.partner_gene if first.side == "five_prime" else first.gene

    queries = [
        f'"{fusion_hyphen}" AND (cancer OR tumor OR tumour OR carcinoma OR sarcoma OR neoplasm)',
    ]
    if five_gene and three_gene:
        queries.extend([
            f'"{five_gene}" AND "{three_gene}" AND fusion',
            f'"{five_gene}" AND "{three_gene}" AND (oncogenic OR kinase OR inhibitor OR response)',
        ])
    kinase_genes = {
        context.kinase_gene
        for context in fusion_contexts
        if context.kinase_gene
    }
    for kinase_gene in sorted(kinase_genes):
        queries.append(f'"{fusion_hyphen}" AND "{kinase_gene}" AND "kinase domain"')
    return list(dict.fromkeys(queries))


def retrieve_pubmed_records_for_queries(
    subject: str,
    queries: list[str],
    *,
    ncbi_api_key: str = "",
    ncbi_client: Optional[NcbiClient] = None,
    max_results: int = 8,
) -> list[PubMedRecord]:
    ncbi_client = ncbi_client or NcbiClient(api_key=ncbi_api_key)
    cache_key = {
        "version": "pubmed-records-v2",
        "subject": subject.upper(),
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
        search = ncbi_client.get(ESEARCH_URL, params=params, timeout=15)
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
    fetch = ncbi_client.get(EFETCH_URL, params=fetch_params, timeout=30)

    records: list[PubMedRecord] = []
    root = ET.fromstring(fetch.text)
    for article in root.findall(".//PubmedArticle"):
        pmid_el = article.find(".//PMID")
        pmid = pmid_el.text if pmid_el is not None else None
        title_el = article.find(".//ArticleTitle")
        title = "".join(title_el.itertext()).strip() if title_el is not None else ""
        abstract_parts = article.findall(".//AbstractText")
        abstract = " ".join(
            "".join(part.itertext()).strip()
            for part in abstract_parts
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


def retrieve_pubmed_records(
    gene: str,
    *,
    ncbi_api_key: str = "",
    ncbi_client: Optional[NcbiClient] = None,
    max_results: int = 8,
    fusion_contexts: Optional[list[FusionCurationContext]] = None,
) -> list[PubMedRecord]:
    fusion_contexts = fusion_contexts or []
    return retrieve_pubmed_records_for_queries(
        gene,
        _pubmed_queries(gene, fusion_contexts),
        ncbi_api_key=ncbi_api_key,
        ncbi_client=ncbi_client,
        max_results=max_results,
    )


def _oncokb_api_token() -> str:
    return os.environ.get("ONCOKB_API_TOKEN", "") or os.environ.get("ONCOKB_TOKEN", "")


def fetch_oncokb_curated_gene_index(api_token: str = "") -> dict[str, dict]:
    if not api_token:
        return {}

    base_url = os.environ.get("ONCOKB_API_BASE_URL", DEFAULT_ONCOKB_API_BASE_URL).rstrip("/")
    params = {"includeEvidence": "true"}
    version = os.environ.get("ONCOKB_DATA_VERSION")
    if version:
        params["version"] = version

    cache_key = {
        "version": "oncokb-curated-genes-v1",
        "base_url": base_url,
        "data_version": version or "",
        "includeEvidence": True,
    }
    ttl_seconds = max(
        0,
        int(os.environ.get("FUSION_GENE_CURATION_ONCOKB_CACHE_TTL_SECONDS", "86400")),
    )
    cached = _read_cache("oncokb-curated-genes", cache_key, ttl_seconds=ttl_seconds)
    if isinstance(cached, dict):
        return {
            str(symbol).upper(): item
            for symbol, item in cached.items()
            if isinstance(item, dict)
        }

    response = requests.get(
        f"{base_url}/utils/allCuratedGenes",
        params=params,
        headers={
            "accept": "application/json",
            "Authorization": f"Bearer {api_token}",
        },
        timeout=30,
    )
    try:
        response.raise_for_status()
    except requests.HTTPError as exc:
        raise OncoKBGeneLookupError("OncoKB curated-gene lookup failed.") from exc

    data = response.json()
    if not isinstance(data, list):
        raise OncoKBGeneLookupError("OncoKB curated-gene lookup returned an unexpected response.")

    genes = {}
    for item in data:
        if not isinstance(item, dict):
            continue
        symbol = str(item.get("hugoSymbol") or "").strip().upper()
        if symbol:
            genes[symbol] = item
    _write_cache("oncokb-curated-genes", cache_key, genes, ttl_seconds=ttl_seconds)
    return genes


def _clean_oncokb_text(value: object) -> str:
    return " ".join(str(value or "").split())


def _truncate_at_sentence(text: str, max_chars: int = 420) -> str:
    if len(text) <= max_chars:
        return text
    clipped = text[:max_chars].rsplit(" ", 1)[0].strip()
    for marker in (". ", "; "):
        sentence = clipped.rsplit(marker, 1)[0].strip()
        if len(sentence) >= 80:
            return f"{sentence}."
    return f"{clipped.rstrip('.')}."


def _oncokb_gene_type_label(gene_type: str) -> str:
    labels = {
        "ONCOGENE_AND_TSG": "oncogene and tumor suppressor",
        "ONCOGENE": "oncogene",
        "TSG": "tumor suppressor",
        "INSUFFICIENT_EVIDENCE": "a gene with insufficient cancer-gene evidence",
        "NEITHER": "neither an oncogene nor tumor suppressor",
    }
    return labels.get(gene_type, gene_type.replace("_", " ").lower() or "a curated gene")


def _oncokb_cancer_association(gene_type: str) -> Optional[bool]:
    if gene_type in {"ONCOGENE_AND_TSG", "ONCOGENE", "TSG"}:
        return True
    if gene_type == "NEITHER":
        return False
    return None


def _oncokb_gene_result(
    gene: str,
    record: Optional[dict],
    fusion_contexts: list[FusionCurationContext],
) -> Optional[dict]:
    if not record:
        return None

    gene_type = str(record.get("geneType") or "").strip().upper()
    summary = _clean_oncokb_text(record.get("summary"))
    background = _clean_oncokb_text(record.get("background"))
    highest_sensitive = _clean_oncokb_text(record.get("highestSensitiveLevel"))
    highest_resistance = _clean_oncokb_text(
        record.get("highestResistanceLevel") or record.get("highestResistancLevel")
    )
    if not any([gene_type, summary, background, highest_sensitive, highest_resistance]):
        return None

    label = _oncokb_gene_type_label(gene_type)
    rationale_parts = [f"OncoKB curates {gene} as {label}."]
    if summary:
        rationale_parts.append(summary)
    elif background:
        rationale_parts.append(background)
    levels = []
    if highest_sensitive:
        levels.append(f"highest sensitive level {highest_sensitive}")
    if highest_resistance:
        levels.append(f"highest resistance level {highest_resistance}")
    if levels:
        rationale_parts.append(f"OncoKB reports {', '.join(levels)}.")

    return {
        "gene": gene,
        "cancer_associated": _oncokb_cancer_association(gene_type),
        "rationale": _truncate_at_sentence(" ".join(rationale_parts)),
        "supporting_pmids": [],
        "retrieved_pmids": [],
        "fusion_contexts": _context_dicts(fusion_contexts),
        "insufficient_evidence": gene_type == "INSUFFICIENT_EVIDENCE" and not summary and not background,
        "curation_source": "OncoKB",
        "oncokb_gene_type": gene_type,
        "oncokb_summary": summary,
        "oncokb_background": background,
        "oncokb_highest_sensitive_level": highest_sensitive,
        "oncokb_highest_resistance_level": highest_resistance,
        "oncokb_url": f"https://www.oncokb.org/gene/{gene}",
    }


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


def _fusion_curation_cache_key(
    fusion: str,
    records: list[PubMedRecord],
    model: str,
    fusion_contexts: list[FusionCurationContext],
) -> dict[str, Any]:
    return {
        "version": "fusion-curation-v1",
        "fusion": fusion.upper(),
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


def _no_fusion_pubmed_evidence_result(
    fusion: str,
    fusion_contexts: Optional[list[FusionCurationContext]] = None,
) -> dict:
    return {
        "fusion": fusion,
        "fusion_literature_identified": False,
        "cancer_associated": None,
        "rationale": "No PubMed abstracts were retrieved for this exact fusion.",
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
            (
                f"  Specificity: {context.fusion_specificity}; "
                f"Genome Nexus breakpoint context available: {context.breakpoint_context_available}"
            ),
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
        if context.limitations:
            lines.append(f"  Limitations: {'; '.join(context.limitations)}")
        if context.annotation_error:
            lines.append(f"  Fusion annotation limitation: {context.annotation_error}")
    return "\n".join(lines)


def curate_gene(
    gene: str,
    *,
    anthropic_api_key: str,
    ncbi_api_key: str = "",
    ncbi_client: Optional[NcbiClient] = None,
    model: str = DEFAULT_CURATION_MODEL,
    max_results: int = 8,
    abstract_chars: int = 1200,
    fusion_contexts: Optional[list[FusionCurationContext]] = None,
    oncokb_genes_by_symbol: Optional[dict[str, dict]] = None,
) -> dict:
    fusion_contexts = fusion_contexts or []
    oncokb_result = _oncokb_gene_result(
        gene,
        (oncokb_genes_by_symbol or {}).get(gene.upper()),
        fusion_contexts,
    )
    if oncokb_result:
        return oncokb_result

    if not anthropic_api_key:
        raise GeneCurationUnavailable("ANTHROPIC_API_KEY is not configured for server-side curation.")

    records = retrieve_pubmed_records(
        gene,
        ncbi_api_key=ncbi_api_key,
        ncbi_client=ncbi_client,
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
        cached = _repair_json_rationale(cached)
        cached.setdefault("gene", gene)
        cached.setdefault("retrieved_pmids", [record.pmid for record in records])
        cached["fusion_contexts"] = _context_dicts(fusion_contexts)
        cached.setdefault("curation_source", "PubMed + LLM")
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
When fusion_specificity is gene_pair_only, limit conclusions to gene-pair-level
literature evidence and explicitly avoid claims about the exact exon, protein
junction, retained/lost domains, or kinase-domain retention.
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
        payload = _extract_json_object(text)
    except json.JSONDecodeError:
        payload = _normal_model_parse_fallback(
            subject_label="gene",
            subject=gene,
            records=records,
            fusion_contexts=fusion_contexts,
            fusion_result=False,
            raw_text=text,
        )
    payload = _repair_json_rationale(payload)
    payload.setdefault("gene", gene)
    payload.setdefault("retrieved_pmids", [record.pmid for record in records])
    payload["fusion_contexts"] = _context_dicts(fusion_contexts)
    payload.setdefault("curation_source", "PubMed + LLM")
    _write_cache("gene-curation", cache_key, payload, ttl_seconds=result_ttl_seconds)
    return payload


def curate_fusion(
    fusion: str,
    *,
    anthropic_api_key: str,
    ncbi_api_key: str = "",
    ncbi_client: Optional[NcbiClient] = None,
    model: str = DEFAULT_CURATION_MODEL,
    max_results: int = 8,
    abstract_chars: int = 1200,
    fusion_contexts: Optional[list[FusionCurationContext]] = None,
) -> dict:
    if not anthropic_api_key:
        raise GeneCurationUnavailable("ANTHROPIC_API_KEY is not configured for server-side curation.")

    fusion_contexts = fusion_contexts or []
    records = retrieve_pubmed_records_for_queries(
        fusion,
        _fusion_pubmed_queries(fusion, fusion_contexts),
        ncbi_api_key=ncbi_api_key,
        ncbi_client=ncbi_client,
        max_results=max_results,
    )
    if not records:
        return _no_fusion_pubmed_evidence_result(fusion, fusion_contexts)

    cache_key = _fusion_curation_cache_key(fusion, records, model, fusion_contexts)
    cache_key["abstract_chars"] = abstract_chars
    result_ttl_seconds = max(
        0,
        int(os.environ.get("FUSION_GENE_CURATION_RESULT_CACHE_TTL_SECONDS", "2592000")),
    )
    cached = _read_cache("fusion-curation", cache_key, ttl_seconds=result_ttl_seconds)
    if isinstance(cached, dict):
        cached = _repair_json_rationale(cached)
        cached.setdefault("fusion", fusion)
        cached.setdefault("retrieved_pmids", [record.pmid for record in records])
        cached["fusion_contexts"] = _context_dicts(fusion_contexts)
        return cached

    context = "\n\n".join(
        f"PMID {record.pmid}\nTitle: {record.title}\nAbstract: {record.abstract[:abstract_chars]}"
        for record in records
    )
    prompt = f"""\
Fusion: {fusion}

Genome Nexus fusion-position context:
{_format_fusion_contexts_for_prompt(fusion_contexts)}

PubMed context:
{context}

Return one JSON object with:
- fusion
- fusion_literature_identified: true/false/null
- cancer_associated: true/false/null
- rationale: concise curator-facing scan text grounded only in the PubMed context.
  Write 1-2 short sentences, ideally 45-90 words total. State whether the exact
  fusion is described in the literature, the strongest cancer/mechanistic context,
  and one caveat if needed.
- supporting_pmids: up to 4 PMIDs from the context
- retrieved_pmids: all PMIDs provided in the context
- fusion_contexts: echo the provided fusion context as compact JSON-compatible objects
- insufficient_evidence: true when the exact fusion literature context is too sparse

Prioritize evidence about the exact fusion over separate evidence about either
partner gene. Do not infer transcript, breakpoint, domain-retention, or
kinase-domain status beyond the Genome Nexus context.
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
        payload = _extract_json_object(text)
    except json.JSONDecodeError:
        payload = _normal_model_parse_fallback(
            subject_label="fusion",
            subject=fusion,
            records=records,
            fusion_contexts=fusion_contexts,
            fusion_result=True,
            raw_text=text,
        )
    payload = _repair_json_rationale(payload)
    payload.setdefault("fusion", fusion)
    payload.setdefault("retrieved_pmids", [record.pmid for record in records])
    payload["fusion_contexts"] = _context_dicts(fusion_contexts)
    _write_cache("fusion-curation", cache_key, payload, ttl_seconds=result_ttl_seconds)
    return payload


def _fusion_curation_is_sufficient(result: dict) -> bool:
    if result.get("error") or result.get("insufficient_evidence"):
        return False
    if result.get("fusion_literature_identified") is False:
        return False
    return bool(result.get("supporting_pmids") or result.get("rationale"))


def curate_fusion_genes(
    fusions: Iterable[object],
    annotation_results: Optional[Iterable[dict]] = None,
    *,
    force_gene_curation: bool = False,
    requested_genes: Optional[Iterable[str]] = None,
) -> dict:
    anthropic_api_key = os.environ.get("ANTHROPIC_API_KEY", "")
    if not anthropic_api_key:
        raise GeneCurationUnavailable("ANTHROPIC_API_KEY is not configured for server-side curation.")

    fusions = list(fusions)
    contexts_by_gene = fusion_contexts_by_gene(fusions, annotation_results)
    contexts_by_fusion = fusion_contexts_by_fusion(fusions, annotation_results)
    ncbi_api_key = os.environ.get("NCBI_API_KEY", "")
    ncbi_client = NcbiClient(api_key=ncbi_api_key)
    model = os.environ.get("FUSION_GENE_CURATION_MODEL", DEFAULT_CURATION_MODEL)
    max_results = max(1, int(os.environ.get("FUSION_GENE_CURATION_MAX_RESULTS", "8")))
    abstract_chars = max(200, int(os.environ.get("FUSION_GENE_CURATION_ABSTRACT_CHARS", "1200")))
    max_workers = max(1, int(os.environ.get("FUSION_GENE_CURATION_WORKERS", "3")))
    fusion_results = []
    results = []
    requested_gene_set = {
        str(gene).strip().upper()
        for gene in (requested_genes or [])
        if str(gene).strip()
    }

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        fusion_futures = {
            executor.submit(
                curate_fusion,
                fusion,
                anthropic_api_key=anthropic_api_key,
                ncbi_api_key=ncbi_api_key,
                ncbi_client=ncbi_client,
                model=model,
                max_results=max_results,
                abstract_chars=abstract_chars,
                fusion_contexts=contexts,
            ): fusion
            for fusion, contexts in contexts_by_fusion.items()
        }
        for future in as_completed(fusion_futures):
            fusion = fusion_futures[future]
            try:
                fusion_results.append(future.result())
            except Exception as exc:
                fusion_results.append({
                    "fusion": fusion,
                    "error": str(exc),
                    "insufficient_evidence": True,
                    "supporting_pmids": [],
                    "retrieved_pmids": [],
                    "fusion_contexts": _context_dicts(contexts_by_fusion.get(fusion, [])),
                })

    sufficient_fusions = {
        result.get("fusion")
        for result in fusion_results
        if _fusion_curation_is_sufficient(result)
    }
    genes = []
    if force_gene_curation:
        genes = list(contexts_by_gene) or unique_genes_from_fusions(fusions)
    else:
        seen_genes = set()
        for fusion, contexts in contexts_by_fusion.items():
            if fusion in sufficient_fusions:
                continue
            for context in contexts:
                if context.gene not in seen_genes:
                    seen_genes.add(context.gene)
                    genes.append(context.gene)
    if requested_gene_set:
        genes = [gene for gene in genes if gene.upper() in requested_gene_set]

    oncokb_lookup_error = None
    oncokb_token = _oncokb_api_token()
    oncokb_genes_by_symbol = {}
    if genes:
        try:
            oncokb_genes_by_symbol = fetch_oncokb_curated_gene_index(oncokb_token)
        except OncoKBGeneLookupError as exc:
            if oncokb_token:
                oncokb_lookup_error = str(exc)

    if oncokb_lookup_error:
        for gene in genes:
            results.append({
                "gene": gene,
                "error": oncokb_lookup_error,
                "insufficient_evidence": True,
                "supporting_pmids": [],
                "retrieved_pmids": [],
                "fusion_contexts": _context_dicts(contexts_by_gene.get(gene, [])),
                "curation_source": "OncoKB",
            })
        results.sort(key=lambda item: item.get("gene", ""))
        fusion_results.sort(key=lambda item: item.get("fusion", ""))
        return {"fusions": fusion_results, "genes": results}

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                curate_gene,
                gene,
                anthropic_api_key=anthropic_api_key,
                ncbi_api_key=ncbi_api_key,
                ncbi_client=ncbi_client,
                model=model,
                max_results=max_results,
                abstract_chars=abstract_chars,
                fusion_contexts=contexts_by_gene.get(gene, []),
                oncokb_genes_by_symbol=oncokb_genes_by_symbol,
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
    fusion_results.sort(key=lambda item: item.get("fusion", ""))
    return {"fusions": fusion_results, "genes": results}
