"""A standalone ``DataProvider`` that talks directly to public REST/GraphQL
annotation sources over HTTPS — no Claude Science ``host.mcp`` dependency.

Use this in any deployment that runs outside a Claude Science session (e.g. the
public MCP server under ``server/app.py``, a cron job, a plain script). Inside a
Claude Science ``repl`` cell, prefer ``MCPDataProvider`` — it goes through the
connected, rate-limit-aware MCP servers instead of hitting upstream APIs raw.

Sources:
  - Ensembl REST   (https://rest.ensembl.org)      — transcript structure + CDS/protein sequence
  - InterPro API   (https://www.ebi.ac.uk/interpro) — protein domain architecture
  - CIViC GraphQL  (https://civicdb.org/api/graphql) — curated fusion knowledge

CIViC's public GraphQL endpoint answers unauthenticated reads for the queries
used here (molecular-profile search, evidence-item listing); set the
``CIVIC_API_KEY`` environment variable to send a bearer token if you have one
(higher rate limit / access to unreleased curation). No key is required to run.
"""
from __future__ import annotations
import os
import time
from typing import Optional

import requests

from .core import Transcript, build_exon_cds_map, build_exon_genomic_map

# Ensembl serves each human genome build from a different host. The default
# rest.ensembl.org is GRCh38; GRCh37/hg19 has its own archive endpoint. A
# genomic breakpoint is only meaningful against the build it was called on, so
# the assembly must be selectable (see the `assembly` arg of RestDataProvider).
ENSEMBL_BASE_GRCH38 = "https://rest.ensembl.org"
ENSEMBL_BASE_GRCH37 = "https://grch37.rest.ensembl.org"
ENSEMBL_BASE = ENSEMBL_BASE_GRCH38  # back-compat default

_ASSEMBLY_BASES = {"GRCh38": ENSEMBL_BASE_GRCH38, "GRCh37": ENSEMBL_BASE_GRCH37}
# Common spellings callers use for each build, normalized to our canonical key.
_ASSEMBLY_ALIASES = {
    "GRCH38": "GRCh38", "HG38": "GRCh38", "38": "GRCh38",
    "GRCH37": "GRCh37", "HG19": "GRCh37", "37": "GRCh37",
}

INTERPRO_BASE = "https://www.ebi.ac.uk/interpro/api"
CIVIC_GRAPHQL = "https://civicdb.org/api/graphql"

_UA = {"User-Agent": "fusion-annotation/0.1 (+https://github.com/genome-nexus/fusion-annotation)"}


def normalize_assembly(name: str | None) -> str:
    """Normalize a genome-build name to a canonical 'GRCh38' / 'GRCh37'.

    Accepts common aliases (hg38, hg19, "37", "38", any case). Returns 'GRCh38'
    for None or empty input. Raises ValueError on anything unrecognized so a typo
    can't silently fall back to the wrong build.
    """
    if not name:
        return "GRCh38"
    key = str(name).strip().upper()
    if key not in _ASSEMBLY_ALIASES:
        raise ValueError(
            f"unknown genome build {name!r}; use 'GRCh38' (aliases: hg38) "
            "or 'GRCh37' (aliases: hg19)")
    return _ASSEMBLY_ALIASES[key]

# Ensembl's REST API is noticeably less reliable from shared cloud NAT egress
# IPs (Cloud Run, etc.) than from a residential/office IP — transient 500s,
# 503s, and read timeouts are common even though the service is otherwise up.
# Retry transient failures a few times with exponential backoff before giving
# up; a 4xx (bad gene symbol, etc.) is a real error and is not retried.
_RETRYABLE_STATUS = {429, 500, 502, 503, 504}
_MAX_ATTEMPTS = 4
_BACKOFF_BASE = 1.5  # seconds; attempt i waits _BACKOFF_BASE * 2**i


def _request_with_retry(method: str, url: str, **kwargs) -> requests.Response:
    kwargs.setdefault("timeout", 30)
    last_exc: Exception | None = None
    for attempt in range(_MAX_ATTEMPTS):
        try:
            r = requests.request(method, url, **kwargs)
            if r.status_code in _RETRYABLE_STATUS and attempt < _MAX_ATTEMPTS - 1:
                time.sleep(_BACKOFF_BASE * (2 ** attempt))
                continue
            r.raise_for_status()
            return r
        except (requests.exceptions.Timeout, requests.exceptions.ConnectionError) as e:
            last_exc = e
            if attempt < _MAX_ATTEMPTS - 1:
                time.sleep(_BACKOFF_BASE * (2 ** attempt))
                continue
    raise last_exc if last_exc else RuntimeError(f"request to {url} failed after {_MAX_ATTEMPTS} attempts")


def _ensembl_get(path: str, base: str = ENSEMBL_BASE, **params) -> dict:
    params.setdefault("content-type", "application/json")
    r = _request_with_retry("GET", f"{base}{path}", params=params, headers=_UA)
    return r.json()


class RestDataProvider:
    """Live provider backed directly by public REST/GraphQL APIs (no MCP).

    Parameters
    ----------
    species : str
        Ensembl species for symbol lookups.
    assembly : str
        Human genome build — 'GRCh38' (default) or 'GRCh37' (aliases hg38/hg19).
        Selects which Ensembl REST host transcript structure and coordinates come
        from; a genomic breakpoint is interpreted against THIS build, so it must
        match the build the breakpoint was called on.
    """

    def __init__(self, species: str = "homo_sapiens", assembly: str = "GRCh38"):
        self.species = species
        self.assembly = normalize_assembly(assembly)
        self.ensembl_base = _ASSEMBLY_BASES[self.assembly]

    # ---- Layer 1 inputs: transcript structure + sequences -----------------
    def get_transcript(self, gene_or_tx: str) -> Transcript:
        def eget(path: str, **params) -> dict:
            return _ensembl_get(path, self.ensembl_base, **params)

        user_pinned_tx = gene_or_tx.upper().startswith("ENS")
        rec = eget(f"/lookup/id/{gene_or_tx}", expand=1) if user_pinned_tx \
            else eget(f"/lookup/symbol/{self.species}/{gene_or_tx}", expand=1)

        if rec.get("object_type") == "Gene":
            tx_id = rec["canonical_transcript"].split(".")[0]
            gene_id = rec["id"]
            gene_symbol = rec.get("display_name", gene_or_tx)
            rec = eget(f"/lookup/id/{tx_id}", expand=1)
            is_canonical = True
        else:
            # Reached only when a transcript id was passed directly (user-pinned);
            # a symbol resolves to a Gene and is handled above.
            tx_id = rec["id"]
            gene_id = rec.get("Parent", "")
            gene_symbol = gene_or_tx
            is_canonical = None

        tr = rec["Translation"]
        cds = eget(f"/sequence/id/{rec['id']}", type="cds")["seq"]
        prot = eget(f"/sequence/id/{rec['id']}", type="protein")["seq"]
        exon_cds = build_exon_cds_map(rec["strand"], rec["Exon"], tr["start"], tr["end"])
        exon_genomic = build_exon_genomic_map(rec["strand"], rec["Exon"])

        uniprot = None
        try:
            xr = eget(f"/xrefs/id/{gene_id or rec['id']}", external_db="Uniprot_gn")
            ids = [x.get("primary_id") for x in xr if x.get("primary_id")]
            uniprot = next((i for i in ids if i and i[0] in "OPQ" and len(i) == 6), ids[0] if ids else None)
        except Exception:
            pass

        return Transcript(
            gene_symbol=gene_symbol, gene_id=gene_id, transcript_id=rec["id"],
            strand=rec["strand"], cds=cds, protein=prot, uniprot=uniprot,
            exon_cds=exon_cds, exon_genomic=exon_genomic,
            cds_g_start=tr["start"], cds_g_end=tr["end"], is_canonical=is_canonical)

    # ---- Layer 1: domains --------------------------------------------------
    def get_domains(self, uniprot: str) -> list[dict]:
        if not uniprot:
            return []
        out, url = [], f"{INTERPRO_BASE}/entry/interpro/protein/uniprot/{uniprot}"
        params = {"page_size": 100}
        while url:
            r = _request_with_retry("GET", url, params=params, headers=_UA)
            payload = r.json()
            for res in payload.get("results", []):
                meta = res["metadata"]
                for prot_entry in res.get("proteins", []):
                    for loc in prot_entry.get("entry_protein_locations", []):
                        for frag in loc.get("fragments", []):
                            out.append({"accession": meta["accession"], "name": meta["name"],
                                        "type": meta["type"], "start": frag["start"], "end": frag["end"]})
            url, params = payload.get("next"), None  # `next` is a full URL already
        return out

    # ---- Layer 3: knowledge (CIViC) ----------------------------------------
    def get_fusion_knowledge(self, categorical_key: str) -> dict:
        five, three = categorical_key.split("::")
        out = {"oncogenic": None, "therapies": [], "evidence": [], "diseases": [], "sources": []}
        headers = dict(_UA)
        headers["Content-Type"] = "application/json"
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
            r = _request_with_retry("POST", CIVIC_GRAPHQL, headers=headers,
                                     json={"query": mp_query, "variables": {"name": f"{five}::{three} Fusion"}})
            nodes = r.json().get("data", {}).get("molecularProfiles", {}).get("nodes", [])
            mp = next((m for m in nodes if m.get("name", "").startswith(f"{five}::{three}")), None)
            if mp:
                out["sources"].append(f"CIViC MP {mp['id']}")
                ev_query = """
                query($mpId: Int!) {
                  evidenceItems(molecularProfileId: $mpId, first: 100) {
                    nodes {
                      id evidenceType evidenceLevel evidenceDirection significance: evidenceDirection
                      disease { name }
                      therapies { name }
                    }
                  }
                }"""
                r = _request_with_retry("POST", CIVIC_GRAPHQL, headers=headers,
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
                        "therapies": [t["name"] for t in (rec.get("therapies") or [])]})
                out["therapies"] = sorted(thx)
                out["diseases"] = sorted(dis)
        except Exception as e:
            out["sources"].append(f"CIViC error: {e}")
        return out
