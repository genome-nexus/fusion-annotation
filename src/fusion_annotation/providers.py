"""Data providers that feed the annotation engines.

Two implementations:

  * ``StaticProvider``  — offline, seeded from a dict. Used by the test suite and
    for reproducible examples. No network.

  * ``MCPDataProvider`` — wires the ``DataProvider`` protocol to live annotation
    sources through a Model Context Protocol (MCP) ``host.mcp`` callable:
        - Ensembl REST   (genomes server):   transcript structure + CDS/protein
        - InterPro       (protein-annotation server): domain architecture
        - CIViC + Open Targets (clinical-genomics server): fusion knowledge

    NOTE: ``host.mcp`` is only available inside the Claude Science ``repl`` tool,
    not the ``python`` kernel. Construct ``MCPDataProvider(host.mcp)`` in a repl
    cell (or pass any callable with signature ``mcp(server, method, **kwargs)``).
"""
from __future__ import annotations
from typing import Callable, Optional
from .core import Transcript, build_exon_cds_map


class StaticProvider:
    """Offline provider seeded from plain dicts — deterministic, no network."""

    def __init__(self, transcripts: dict[str, Transcript],
                 domains: Optional[dict[str, list[dict]]] = None,
                 knowledge: Optional[dict[str, dict]] = None):
        self._tx = transcripts
        self._dom = domains or {}
        self._kn = knowledge or {}

    def get_transcript(self, gene_or_tx: str) -> Transcript:
        if gene_or_tx not in self._tx:
            raise KeyError(f"no transcript seeded for {gene_or_tx!r}")
        return self._tx[gene_or_tx]

    def get_domains(self, uniprot: str) -> list[dict]:
        return self._dom.get(uniprot, [])

    def get_fusion_knowledge(self, categorical_key: str) -> dict:
        return self._kn.get(categorical_key, {})


class MCPDataProvider:
    """Live provider backed by MCP connectors.

    Parameters
    ----------
    mcp : callable
        ``host.mcp`` (or equivalent) with signature ``mcp(server, method, **kwargs)``.
    species : str
        Ensembl species for symbol lookups.
    """

    def __init__(self, mcp: Callable, species: str = "homo_sapiens"):
        self.mcp = mcp
        self.species = species

    # ---- Layer 1 inputs: transcript structure + sequences -----------------
    def get_transcript(self, gene_or_tx: str) -> Transcript:
        # Resolve a symbol to its canonical transcript, or accept a transcript id.
        rec = self.mcp("genomes", "ensembl_lookup", query=gene_or_tx,
                       species=self.species, expand=True)["record"]
        if rec.get("object_type") == "Gene" or "Transcript" in rec:
            # a gene was returned; descend to the canonical transcript
            tx_id = rec.get("canonical_transcript", "").split(".")[0]
            gene_id = rec["id"]
            gene_symbol = rec.get("display_name", gene_or_tx)
            rec = self.mcp("genomes", "ensembl_lookup", query=tx_id, expand=True)["record"]
        else:
            tx_id = rec["id"]
            gene_id = rec.get("Parent", "")
            gene_symbol = gene_or_tx

        tr = rec["Translation"]
        cds = self.mcp("genomes", "ensembl_sequence", stable_id=rec["id"], seq_type="cds")["seq"]
        prot = self.mcp("genomes", "ensembl_sequence", stable_id=rec["id"], seq_type="protein")["seq"]
        exon_cds = build_exon_cds_map(rec["strand"], rec["Exon"], tr["start"], tr["end"])

        # best-effort UniProt xref for the gene (for InterPro domains)
        uniprot = None
        try:
            xr = self.mcp("genomes", "ensembl_xrefs", stable_id=gene_id or rec["id"],
                          external_db="Uniprot_gn").get("xrefs", [])
            # prefer a reviewed-looking primary id (SwissProt style: [OPQ]….)
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
        res = self.mcp("protein-annotation", "get_domain_architecture", accessions=[uniprot])
        summ = res.get("summaries", {}).get(uniprot, {})
        out = []
        for e in summ.get("entries", []):
            for loc in e.get("locations", []):
                for f in loc.get("fragments", []):
                    out.append({"accession": e["accession"], "name": e["name"],
                                "type": e["type"], "start": f["start"], "end": f["end"]})
        return out

    # ---- Layer 3: knowledge (CIViC + Open Targets) -------------------------
    def get_fusion_knowledge(self, categorical_key: str) -> dict:
        five, three = categorical_key.split("::")
        out = {"oncogenic": None, "therapies": [], "evidence": [],
               "diseases": [], "sources": []}
        # Find the categorical molecular profile in CIViC (e.g. "EML4::ALK Fusion").
        try:
            mps = self.mcp("clinical-genomics", "civic_search_molecular_profiles",
                           name=f"{five}::{three} Fusion").get("records", [])
            mp = next((m for m in mps if m.get("name", "").startswith(f"{five}::{three}")), None)
            if mp:
                out["sources"].append(f"CIViC MP {mp['id']}")
                ev = self.mcp("clinical-genomics", "civic_search_evidence",
                              molecular_profile_id=mp["id"]).get("records", [])
                thx, dis = set(), set()
                for r in ev:
                    for t in (r.get("therapies") or []):
                        thx.add(t["name"])
                    d = r.get("disease")
                    if isinstance(d, dict) and d.get("name"):
                        dis.add(d["name"])
                    out["evidence"].append({
                        "id": r.get("id"), "type": r.get("evidenceType"),
                        "level": r.get("evidenceLevel"), "significance": r.get("significance"),
                        "disease": d.get("name") if isinstance(d, dict) else d,
                        "therapies": [t["name"] for t in (r.get("therapies") or [])]})
                out["therapies"] = sorted(thx)
                out["diseases"] = sorted(dis)
        except Exception as e:
            out["sources"].append(f"CIViC error: {e}")
        return out
