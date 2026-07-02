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

from .core import Transcript, build_exon_cds_map

ENSEMBL_BASE = "https://rest.ensembl.org"
INTERPRO_BASE = "https://www.ebi.ac.uk/interpro/api"
CIVIC_GRAPHQL = "https://civicdb.org/api/graphql"

_UA = {"User-Agent": "fusion-annotation/0.1 (+https://github.com/genome-nexus/fusion-annotation)"}


def _ensembl_get(path: str, **params) -> dict:
    params.setdefault("content-type", "application/json")
    r = requests.get(f"{ENSEMBL_BASE}{path}", params=params, headers=_UA, timeout=20)
    r.raise_for_status()
    return r.json()


class RestDataProvider:
    """Live provider backed directly by public REST/GraphQL APIs (no MCP)."""

    def __init__(self, species: str = "homo_sapiens"):
        self.species = species

    # ---- Layer 1 inputs: transcript structure + sequences -----------------
    def get_transcript(self, gene_or_tx: str) -> Transcript:
        rec = _ensembl_get(f"/lookup/symbol/{self.species}/{gene_or_tx}", expand=1) \
            if not gene_or_tx.upper().startswith("ENS") \
            else _ensembl_get(f"/lookup/id/{gene_or_tx}", expand=1)

        if rec.get("object_type") == "Gene":
            tx_id = rec["canonical_transcript"].split(".")[0]
            gene_id = rec["id"]
            gene_symbol = rec.get("display_name", gene_or_tx)
            rec = _ensembl_get(f"/lookup/id/{tx_id}", expand=1)
        else:
            tx_id = rec["id"]
            gene_id = rec.get("Parent", "")
            gene_symbol = gene_or_tx

        tr = rec["Translation"]
        cds = _ensembl_get(f"/sequence/id/{rec['id']}", type="cds")["seq"]
        prot = _ensembl_get(f"/sequence/id/{rec['id']}", type="protein")["seq"]
        exon_cds = build_exon_cds_map(rec["strand"], rec["Exon"], tr["start"], tr["end"])

        uniprot = None
        try:
            xr = _ensembl_get(f"/xrefs/id/{gene_id or rec['id']}", external_db="Uniprot_gn")
            ids = [x.get("primary_id") for x in xr if x.get("primary_id")]
            uniprot = next((i for i in ids if i and i[0] in "OPQ" and len(i) == 6), ids[0] if ids else None)
        except Exception:
            pass

        return Transcript(
            gene_symbol=gene_symbol, gene_id=gene_id, transcript_id=rec["id"],
            strand=rec["strand"], cds=cds, protein=prot, uniprot=uniprot, exon_cds=exon_cds)

    # ---- Layer 1: domains --------------------------------------------------
    def get_domains(self, uniprot: str) -> list[dict]:
        if not uniprot:
            return []
        out, url = [], f"{INTERPRO_BASE}/entry/interpro/protein/uniprot/{uniprot}"
        params = {"page_size": 100}
        while url:
            r = requests.get(url, params=params, headers=_UA, timeout=20)
            r.raise_for_status()
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
            r = requests.post(CIVIC_GRAPHQL, headers=headers, timeout=20,
                               json={"query": mp_query, "variables": {"name": f"{five}::{three} Fusion"}})
            r.raise_for_status()
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
                r = requests.post(CIVIC_GRAPHQL, headers=headers, timeout=20,
                                   json={"query": ev_query, "variables": {"mpId": mp["id"]}})
                r.raise_for_status()
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
