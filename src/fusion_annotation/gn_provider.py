"""DataProvider backed by Genome Nexus + UCSC + CIViC (no Ensembl REST).

Drop-in replacement for ``RestDataProvider`` that eliminates the slow Ensembl
REST dependency.  Per-annotation latency drops from ~50–95 s to ~1–2 s by:

  * Using Genome Nexus (cached Ensembl mirror) for transcript metadata + Pfam
    domains in two calls vs. five raw Ensembl calls.
  * Fetching the CDS nucleotide sequence from UCSC in one regional call and
    assembling it from the exon table, removing the two Ensembl /sequence calls.
  * Parallelising the two partner lookups with ``ThreadPoolExecutor``.
  * Resolving chromosome from HGNC/HUGO symbol metadata with NCBI esummary as
    a fallback.
  * Keeping the CIViC knowledge query unchanged.

Sources
-------
  - Genome Nexus (https://grch38.genomenexus.org / https://www.genomenexus.org)
  - UCSC Genome Browser REST API (https://api.genome.ucsc.edu)
  - HGNC REST API (https://rest.genenames.org) — primary symbol chrom lookup
  - NCBI eutils (https://eutils.ncbi.nlm.nih.gov)  — fallback chrom lookup only
  - CIViC GraphQL (https://civicdb.org/api/graphql)
  - InterPro API (https://www.ebi.ac.uk/interpro) — optional enrichment, on by default

Build mapping
-------------
  GRCh38  →  grch38.genomenexus.org  +  UCSC hg38
  GRCh37  →  www.genomenexus.org     +  UCSC hg19

See issue #13 for the full design rationale.
"""
from __future__ import annotations

import json
import os
import re
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, wait as futures_wait, FIRST_EXCEPTION
from typing import Optional

from .core import Transcript, translate, build_exon_cds_map, build_exon_genomic_map
from .rest_provider import (
    RestDataProvider, normalize_assembly, _request_with_retry, _UA,
    CIVIC_GRAPHQL, INTERPRO_BASE,
)

# ---------------------------------------------------------------------------
# Host / build routing
# ---------------------------------------------------------------------------
_GN_BASES = {
    "GRCh38": "https://grch38.genomenexus.org",
    "GRCh37": "https://www.genomenexus.org",
}
_UCSC_GENOMES = {"GRCh38": "hg38", "GRCh37": "hg19"}
_UCSC_SEQ_BASE = "https://api.genome.ucsc.edu/getData/sequence"
_HGNC_BASE = "https://rest.genenames.org"
_NCBI_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
_NCBI_ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

# Threshold above which we switch to per-exon UCSC fetches instead of one
# region call (avoids multi-MB payloads for very large genes).
_REGION_FETCH_MAX_SPAN_BP = 2_000_000

_HGNC_LOCATION_RE = re.compile(r"^(?:chr)?(?P<chrom>[1-9]|1[0-9]|2[0-2]|X|Y|MT|M)(?=[pq]|$)", re.IGNORECASE)


def _chrom_from_hgnc_location(location: str | None) -> str | None:
    """Convert an HGNC cytogenetic location like '17q21.31' to UCSC chrom."""
    if not location:
        return None
    match = _HGNC_LOCATION_RE.match(location.strip())
    if not match:
        return None
    chrom = match.group("chrom").upper()
    return _chrom_to_ucsc("MT" if chrom == "M" else chrom)


def _hgnc_docs(field: str, term: str) -> list[dict]:
    encoded = urllib.parse.quote(term, safe="")
    url = f"{_HGNC_BASE}/fetch/{field}/{encoded}"
    req = urllib.request.Request(url, headers={**dict(_UA), "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=10) as resp:
        payload = json.loads(resp.read())
    return payload.get("response", {}).get("docs", [])


def _best_hgnc_doc(docs: list[dict], symbol: str) -> dict | None:
    if not docs:
        return None
    sym = symbol.upper()
    for doc in docs:
        if doc.get("status") == "Approved" and doc.get("symbol", "").upper() == sym:
            return doc
    for doc in docs:
        if doc.get("status") == "Approved":
            return doc
    return docs[0]


def _hgnc_chrom(symbol: str) -> str | None:
    """Resolve approved or alias HGNC symbol to UCSC chromosome."""
    try:
        doc = _best_hgnc_doc(_hgnc_docs("symbol", symbol), symbol)
        if not doc:
            alias_docs = _hgnc_docs("alias_symbol", symbol) + _hgnc_docs("prev_symbol", symbol)
            doc = _best_hgnc_doc(alias_docs, symbol)
        return _chrom_from_hgnc_location(doc.get("location") if doc else None)
    except Exception:
        return None


def _ncbi_chrom(symbol: str) -> str | None:
    """Resolve gene symbol → chromosome via NCBI esearch + esummary (fallback)."""
    try:
        params = f"db=gene&term={symbol}[Gene+Name]+AND+human[Organism]+AND+alive[prop]&retmode=json&retmax=1"
        req = urllib.request.Request(f"{_NCBI_ESEARCH}?{params}", headers=dict(_UA))
        with urllib.request.urlopen(req, timeout=10) as resp:
            ids = json.loads(resp.read()).get("esearchresult", {}).get("idlist", [])
        if not ids:
            return None
        req2 = urllib.request.Request(f"{_NCBI_ESUMMARY}?db=gene&id={ids[0]}&retmode=json", headers=dict(_UA))
        with urllib.request.urlopen(req2, timeout=10) as resp2:
            doc = json.loads(resp2.read()).get("result", {}).get(ids[0], {})
        return doc.get("chromosome") or None
    except Exception:
        return None


def _chrom_to_ucsc(ncbi_chrom: str) -> str:
    """Convert NCBI chromosome string to UCSC format: '2' → 'chr2', 'MT' → 'chrM'."""
    if ncbi_chrom.upper() == "MT":
        return "chrM"
    return f"chr{ncbi_chrom}"


def _rc(seq: str) -> str:
    """Reverse complement a DNA sequence."""
    return seq.upper().translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _gn_exon_to_dict(exon: dict) -> dict:
    """Convert a GN exon record to the {start, end} shape build_exon_*_map expects."""
    return {"start": exon["exonStart"], "end": exon["exonEnd"]}


def _cds_bounds_from_utrs(utrs: list[dict], strand: int, exons: list[dict]) -> tuple[int, int]:
    """Derive (cds_g_start, cds_g_end) lo<hi from GN UTR records and strand.

    On the plus strand the 5'UTR is at lower coords; on the minus strand the
    5'UTR is at higher coords.  Either way, cds_g_start < cds_g_end.

    A transcript can be CDS-incomplete on one end (no annotated stop/start
    codon, e.g. NTRK1's canonical ENST00000524377 has no three_prime_UTR at
    all) — GN then simply omits that UTR record. When the UTR on a given
    side is missing there's nothing to trim, so that boundary falls back to
    the transcript's outermost exon coordinate on that side instead of
    raising.
    """
    if not utrs:
        raise ValueError("transcript has no UTR records; cannot derive CDS bounds.")
    five_utr = next((u for u in utrs if u["type"] == "five_prime_UTR"), None)
    three_utr = next((u for u in utrs if u["type"] == "three_prime_UTR"), None)
    # Whichever UTR sits at the lower/higher genomic coordinates depends on
    # strand: for +1, 5'UTR is lower and 3'UTR is higher; for -1, reversed.
    lower_utr, higher_utr = (five_utr, three_utr) if strand != -1 else (three_utr, five_utr)

    lo = (lower_utr["end"] + 1) if lower_utr else min(e["exonStart"] for e in exons)
    hi = (higher_utr["start"] - 1) if higher_utr else max(e["exonEnd"] for e in exons)
    return lo, hi


def _ucsc_region_seq(genome: str, chrom: str, lo: int, hi: int) -> str:
    """Fetch a genomic region from UCSC (0-based start, 1-based end = half-open)."""
    url = f"{_UCSC_SEQ_BASE}?genome={genome};chrom={chrom};start={lo - 1};end={hi}"
    req = urllib.request.Request(url, headers=dict(_UA))
    with urllib.request.urlopen(req, timeout=30) as resp:
        return json.loads(resp.read())["dna"]


def _assemble_cds(region_seq: str, region_lo: int, exons: list[dict],
                  cds_g_start: int, cds_g_end: int, strand: int) -> str:
    """Assemble CDS from a genomic region sequence and exon rank table.

    ``region_seq`` is the sequence of ``region_lo..region_hi`` (same as the
    span fetched from UCSC). Coding portions of each exon are sliced in rank
    order; minus-strand exons are reverse-complemented before concatenation.
    """
    exons_ranked = sorted(exons, key=lambda e: e["rank"])
    parts: list[str] = []
    for exon in exons_ranked:
        e_lo = max(exon["exonStart"], cds_g_start)
        e_hi = min(exon["exonEnd"], cds_g_end)
        if e_hi < e_lo:
            continue  # fully non-coding exon
        off_lo = e_lo - region_lo
        off_hi = e_hi - region_lo + 1
        seg = region_seq[off_lo:off_hi].upper()
        parts.append(_rc(seg) if strand == -1 else seg)
    return "".join(parts)


def _pfam_to_domain_dicts(pfam_domains: list[dict]) -> list[dict]:
    """Convert GN pfamDomains entries to the unified domain dict format."""
    out = []
    for d in pfam_domains:
        out.append({
            "accession": d["pfamDomainId"],
            "name": d.get("pfamDomainName") or d["pfamDomainId"],
            "type": "domain",
            "start": d["pfamDomainStart"],
            "end": d["pfamDomainEnd"],
        })
    return out


class GenomeNexusDataProvider:
    """Live provider backed by Genome Nexus + UCSC + CIViC (no Ensembl REST).

    Parameters
    ----------
    assembly : str
        Human genome build — 'GRCh38' (default) or 'GRCh37' (aliases hg38/hg19).
    interpro_enrichment : bool
        Whether to enrich Pfam domains with InterPro entries (default True).
        InterPro failures degrade to empty domains + a warning rather than
        aborting the annotation (issue #9).
    """

    def __init__(self, assembly: str = "GRCh38", interpro_enrichment: bool = True):
        self.assembly = normalize_assembly(assembly)
        self._gn_base = _GN_BASES[self.assembly]
        self._ucsc_genome = _UCSC_GENOMES[self.assembly]
        self._interpro_enrichment = interpro_enrichment
        # Per-instance caches (ephemeral, reset with each new provider instance).
        self._tx_cache: dict[str, dict] = {}     # GN raw transcript → cache
        self._chrom_cache: dict[str, str] = {}   # gene_symbol → NCBI chrom
        self._domain_warning: str | None = None  # set on InterPro degradation

    # ---- internal helpers --------------------------------------------------

    def _gn_get(self, path: str) -> dict:
        url = f"{self._gn_base}{path}"
        r = _request_with_retry("GET", url, headers={**_UA, "Accept": "application/json"})
        return r.json()

    def _resolve_chrom(self, gene_symbol: str) -> str:
        """Return UCSC-format chromosome for gene_symbol (e.g. 'chr2', 'chrX')."""
        sym = gene_symbol.upper()
        if sym in self._chrom_cache:
            return self._chrom_cache[sym]

        hgnc = _hgnc_chrom(gene_symbol)
        if hgnc:
            self._chrom_cache[sym] = hgnc
            return hgnc

        ncbi = _ncbi_chrom(gene_symbol)
        if ncbi:
            self._chrom_cache[sym] = _chrom_to_ucsc(ncbi)
            return self._chrom_cache[sym]
        raise ValueError(
            f"could not determine chromosome for gene {gene_symbol!r}; "
            "ensure HGNC or NCBI eutils is reachable.")

    def _fetch_gn_tx(self, gene_or_tx: str) -> tuple[dict, bool]:
        """Fetch raw GN transcript dict; return (data, is_canonical)."""
        if gene_or_tx in self._tx_cache:
            return self._tx_cache[gene_or_tx]
        user_pinned = gene_or_tx.upper().startswith("ENST")
        if user_pinned:
            data = self._gn_get(f"/ensembl/transcript/{gene_or_tx}")
            result = (data, None)  # is_canonical unknown for user-pinned tx
        else:
            data = self._gn_get(f"/ensembl/canonical-transcript/hgnc/{gene_or_tx}")
            result = (data, True)
        self._tx_cache[gene_or_tx] = result
        # Also cache by transcript ID so get_domains can reuse it
        tx_id = data.get("transcriptId", "")
        if tx_id and tx_id not in self._tx_cache:
            self._tx_cache[tx_id] = result
        return result

    def _fetch_cds_sequence(self, tx_data: dict, gene_symbol: str) -> str:
        """Fetch + assemble the CDS nucleotide sequence for a GN transcript."""
        strand = tx_data["exons"][0]["strand"]
        utrs = tx_data.get("utrs", [])
        cds_g_start, cds_g_end = _cds_bounds_from_utrs(utrs, strand, tx_data["exons"])
        span = cds_g_end - cds_g_start + 1

        chrom = self._resolve_chrom(gene_symbol)

        if span <= _REGION_FETCH_MAX_SPAN_BP:
            # One UCSC call for the whole CDS-spanning region (most genes)
            region_seq = _ucsc_region_seq(self._ucsc_genome, chrom, cds_g_start, cds_g_end)
            cds = _assemble_cds(region_seq, cds_g_start, tx_data["exons"],
                                cds_g_start, cds_g_end, strand)
        else:
            # Per-exon fallback for very large genes (avoids multi-MB payloads)
            exons_ranked = sorted(tx_data["exons"], key=lambda e: e["rank"])
            parts: list[str] = []
            for exon in exons_ranked:
                e_lo = max(exon["exonStart"], cds_g_start)
                e_hi = min(exon["exonEnd"], cds_g_end)
                if e_hi < e_lo:
                    continue
                seg = _ucsc_region_seq(self._ucsc_genome, chrom, e_lo, e_hi)
                parts.append(_rc(seg.upper()) if strand == -1 else seg.upper())
            cds = "".join(parts)

        return cds

    def _rest_fallback_transcript(self, tx_id: str, gene_symbol: str,
                                  is_canonical: Optional[bool]) -> Optional[Transcript]:
        """Fallback for rare GN transcripts whose assembled CDS is unusable.

        Some GN transcripts expose enough structure to resolve breakpoints and
        domains, but the assembled UCSC CDS translates to a protein that is far
        shorter than GN's reported protein length (observed: RAD51, CDH7 on
        GRCh37). In that narrow case, fall back to Ensembl REST for the same
        transcript id so downstream annotation can proceed.
        """
        try:
            fallback = RestDataProvider(assembly=self.assembly).get_transcript(tx_id)
        except Exception:
            return None
        fallback.gene_symbol = gene_symbol or fallback.gene_symbol
        fallback.is_canonical = is_canonical
        return fallback

    # ---- DataProvider protocol --------------------------------------------

    def get_transcript(self, gene_or_tx: str) -> Transcript:
        tx_data, is_canonical = self._fetch_gn_tx(gene_or_tx)

        strand = tx_data["exons"][0]["strand"]  # per-exon strand; top-level is null
        utrs = tx_data.get("utrs", [])
        cds_g_start, cds_g_end = _cds_bounds_from_utrs(utrs, strand, tx_data["exons"])

        # Derive gene symbol (needed for chromosome lookup + Transcript.gene_symbol)
        if gene_or_tx.upper().startswith("ENST"):
            hugo_list = tx_data.get("hugoSymbols") or []
            gene_symbol = hugo_list[0] if hugo_list else gene_or_tx
        else:
            gene_symbol = gene_or_tx

        cds = self._fetch_cds_sequence(tx_data, gene_symbol)
        protein = translate(cds).split("*")[0]

        # Validate length against GN's reported proteinLength
        reported_len = tx_data.get("proteinLength")
        gross_mismatch = (
            bool(reported_len) and
            abs(len(protein) - reported_len) > max(30, reported_len // 2)
        )
        if gross_mismatch:
            fallback = self._rest_fallback_transcript(
                tx_data.get("transcriptId", gene_or_tx), gene_symbol, is_canonical)
            if fallback is not None:
                return fallback
        if reported_len and abs(len(protein) - reported_len) > 1:
            # >1 aa difference is unusual — warn but don't abort
            self._domain_warning = (  # reuse the warning slot; callers check it
                f"translated protein length {len(protein)} differs from GN "
                f"proteinLength {reported_len} for {gene_symbol} "
                f"({tx_data.get('transcriptId')}); check for RNA editing or "
                "selenocysteine recoding")

        # Convert GN exon format to {start, end} for the exon-map builders
        exons_std = [_gn_exon_to_dict(e) for e in tx_data["exons"]]
        exon_cds = build_exon_cds_map(strand, exons_std, cds_g_start, cds_g_end)
        exon_genomic = build_exon_genomic_map(strand, exons_std)

        return Transcript(
            gene_symbol=gene_symbol,
            gene_id=tx_data.get("geneId", ""),
            transcript_id=tx_data.get("transcriptId", gene_or_tx),
            strand=strand,
            cds=cds,
            protein=protein,
            uniprot=tx_data.get("uniprotId"),
            exon_cds=exon_cds,
            exon_genomic=exon_genomic,
            cds_g_start=cds_g_start,
            cds_g_end=cds_g_end,
            is_canonical=is_canonical,
            chrom=self._resolve_chrom(gene_symbol),
        )

    def get_domains(self, uniprot: str) -> list[dict]:
        if not uniprot:
            return []

        # Pull Pfam domains from the cached GN transcript response
        pfam: list[dict] = []
        for cached_val in self._tx_cache.values():
            tx_data, _ = cached_val
            if tx_data.get("uniprotId") == uniprot:
                pfam = _pfam_to_domain_dicts(tx_data.get("pfamDomains") or [])
                break

        if not self._interpro_enrichment:
            return pfam

        # Optional InterPro enrichment (degradable — issue #9)
        try:
            out = list(pfam)
            url: str | None = f"{INTERPRO_BASE}/entry/interpro/protein/uniprot/{uniprot}"
            params: dict | None = {"page_size": 100}
            while url:
                r = _request_with_retry("GET", url, params=params, headers=_UA, timeout=15)
                payload = r.json()
                for res in payload.get("results", []):
                    meta = res["metadata"]
                    for prot_entry in res.get("proteins", []):
                        for loc in prot_entry.get("entry_protein_locations", []):
                            for frag in loc.get("fragments", []):
                                out.append({
                                    "accession": meta["accession"],
                                    "name": meta["name"],
                                    "type": meta["type"],
                                    "start": frag["start"],
                                    "end": frag["end"],
                                })
                url, params = payload.get("next"), None
            return out
        except Exception:
            if pfam:
                self._domain_warning = (
                    f"InterPro enrichment unavailable for {uniprot}; "
                    "Pfam domains only")
            else:
                self._domain_warning = (
                    f"domain annotation degraded (GN Pfam and InterPro both "
                    f"unavailable for {uniprot})")
            return pfam

    def get_fusion_knowledge(self, categorical_key: str) -> dict:
        """CIViC knowledge query — identical to RestDataProvider."""
        five, three = categorical_key.split("::")
        out: dict = {"oncogenic": None, "therapies": [], "evidence": [],
                     "diseases": [], "sources": []}
        headers = {**_UA, "Content-Type": "application/json"}
        api_key = os.environ.get("CIVIC_API_KEY")
        if api_key:
            headers["Authorization"] = f"Bearer {api_key}"
        try:
            mp_query = """
            query($name: String!) {
              molecularProfiles(name: $name, first: 10) {
                nodes { id name }
              }
            }"""
            r = _request_with_retry(
                "POST", CIVIC_GRAPHQL, headers=headers,
                json={"query": mp_query, "variables": {"name": f"{five}::{three} Fusion"}})
            nodes = r.json().get("data", {}).get("molecularProfiles", {}).get("nodes", [])
            mp = next((m for m in nodes if m.get("name", "").startswith(f"{five}::{three}")), None)
            if mp:
                out["sources"].append(f"CIViC MP {mp['id']}")
                ev_query = """
                query($mpId: Int!) {
                  evidenceItems(molecularProfileId: $mpId, first: 100) {
                    nodes {
                      id evidenceType evidenceLevel evidenceDirection
                      disease { name }
                      therapies { name }
                    }
                  }
                }"""
                r = _request_with_retry(
                    "POST", CIVIC_GRAPHQL, headers=headers,
                    json={"query": ev_query, "variables": {"mpId": mp["id"]}})
                ev = r.json().get("data", {}).get("evidenceItems", {}).get("nodes", [])
                thx, dis = set(), set()
                for rec in ev:
                    for t in (rec.get("therapies") or []):
                        thx.add(t["name"])
                    d = rec.get("disease")
                    if isinstance(d, dict) and d.get("name"):
                        dis.add(d["name"])
                    out["evidence"].append({
                        "id": rec.get("id"), "type": rec.get("evidenceType"),
                        "level": rec.get("evidenceLevel"),
                        "disease": d.get("name") if isinstance(d, dict) else d,
                        "therapies": [t["name"] for t in (rec.get("therapies") or [])],
                    })
                out["therapies"] = sorted(thx)
                out["diseases"] = sorted(dis)
        except Exception as e:
            out["sources"].append(f"CIViC error: {e}")
        return out
