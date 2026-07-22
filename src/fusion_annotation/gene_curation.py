"""Server-side literature curation for Genome Nexus fusion batches."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import os
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from dataclasses import asdict
from typing import Any
from typing import Iterable
from typing import Optional

import requests


ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


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
    if (
        iface.get("five_last_aa") is not None
        or iface.get("three_first_aa") is not None
        or iface.get("domains")
    ):
        return "protein_domain_level"

    resolved = result.get("resolved", {})
    if resolved.get("five") or resolved.get("three"):
        return "exon_level"

    return "gene_pair_only"


def _context_limitations(
    *,
    result: Optional[dict],
    error: Optional[str],
) -> tuple[str, ...]:
    if result:
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
    breakpoint_context_available = result is not None
    return FusionCurationContext(
        gene=gene,
        fusion=_fusion_label(fusion),
        side=side,
        partner_gene=partner_gene,
        breakpoint_context_available=breakpoint_context_available,
        fusion_specificity=_fusion_specificity(result),
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


def _context_dicts(contexts: list[FusionCurationContext]) -> list[dict]:
    return [asdict(context) for context in contexts]


def _strip_markdown_json_fence(text: str) -> str:
    cleaned = text.strip()
    if not cleaned.startswith("```"):
        return cleaned
    lines = cleaned.splitlines()
    if len(lines) >= 2 and lines[-1].strip().startswith("```"):
        return "\n".join(lines[1:-1]).strip()
    return cleaned


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
        title = "".join(title_el.itertext()).strip() if title_el is not None else ""
        abstract_parts = article.findall(".//AbstractText")
        abstract = " ".join(
            "".join(part.itertext()).strip()
            for part in abstract_parts
        ).strip()
        if pmid and abstract:
            records.append(PubMedRecord(pmid=pmid, title=title, abstract=abstract))
    return records


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
        payload = json.loads(_strip_markdown_json_fence(text))
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
