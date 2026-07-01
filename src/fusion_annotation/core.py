"""
fusion_annotator.py — a reference implementation of a three-layer gene-fusion
annotation system.

Design (mirrors VEP + OncoKB with an HGVS.p-like interface between them):

    Layer 1  EFFECT      "VEP-like": given an assayed fusion (two transcripts +
                          exon/coordinate breakpoints), reconstruct the chimeric
                          CDS, translate it, and determine reading-frame status,
                          the junction residue, and domain retention.

    Layer 2  INTERFACE   "HGVS.p-like": a normalized, computable representation of
                          the fusion at the protein level — the FusionProtein
                          dataclass and its .to_hgvsp()/.to_dict() serializers.
                          This is the contract that Layer 1 emits and Layer 3
                          consumes. Uses the HGVS/VICC "::" adjoined-junction
                          operator.

    Layer 3  KNOWLEDGE    "OncoKB-like": given a *categorical* fusion key
                          (e.g. "EML4::ALK"), attach curated clinical knowledge
                          (therapies, evidence, oncogenicity). Keyed on the
                          categorical gene pair, not the exact breakpoint.

The layers talk to real annotation sources through a pluggable `DataProvider`
protocol. `MCPDataProvider` wires them to the connected Ensembl / InterPro /
CIViC / Open Targets MCP servers (call from the `repl` tool). This module itself
is pure-Python and has no MCP dependency, so it is unit-testable offline.

Standards references:
  - HGVS Nomenclature: chimeric proteins as delins; "::" junction operator.
  - VICC Gene Fusion Specification (fusions.cancervariants.org): assayed vs.
    categorical fusions; minimal information model.
  - HGNC gene-fusion designation: GENE1::GENE2 with "::" separator.
"""
from __future__ import annotations
from dataclasses import dataclass, field, asdict
from typing import Optional, Protocol, Literal
import json

# ----------------------------------------------------------------------------
# Genetic code (standard) — used to translate the reconstructed chimeric CDS.
# ----------------------------------------------------------------------------
_CODON = {
 'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
 'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
 'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
 'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
 'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
 'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
 'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
 'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}

_AA3 = {'A':'Ala','R':'Arg','N':'Asn','D':'Asp','C':'Cys','Q':'Gln','E':'Glu',
 'G':'Gly','H':'His','I':'Ile','L':'Leu','K':'Lys','M':'Met','F':'Phe','P':'Pro',
 'S':'Ser','T':'Thr','W':'Trp','Y':'Tyr','V':'Val','*':'Ter','X':'Xaa'}

def translate(cds: str) -> str:
    cds = cds.upper()
    return ''.join(_CODON.get(cds[i:i+3], 'X') for i in range(0, len(cds) - len(cds) % 3, 3))

def aa3(one: str) -> str:
    return _AA3.get(one, 'Xaa')


# ----------------------------------------------------------------------------
# Layer-1 inputs: a transcript with enough structure to map exon -> CDS coords.
# ----------------------------------------------------------------------------
@dataclass
class Transcript:
    gene_symbol: str
    gene_id: str
    transcript_id: str
    strand: int
    cds: str                 # coding sequence (ATG..stop), 5'->3' in transcript orientation
    protein: str             # translated protein (no trailing stop)
    uniprot: Optional[str] = None
    # exon_cds[i] = (cds_start, cds_end) 1-based inclusive for coding exon i (rank order);
    # non-coding exons omitted. Built by build_exon_cds_map().
    exon_cds: list[tuple[int, int]] = field(default_factory=list)

    def cds_len(self) -> int:
        return len(self.cds)


def build_exon_cds_map(strand: int, exons: list[dict], cds_g_start: int, cds_g_end: int
                       ) -> list[tuple[int, int]]:
    """Map genomic exons to CDS-relative (start,end) coords, in transcription order.

    exons: list of {"start": g_lo, "end": g_hi} (genomic, lo<hi as Ensembl returns).
    cds_g_start/cds_g_end: genomic bounds of the CDS (lo<hi).
    Returns 1-based inclusive CDS coords for each *coding* exon, in rank order;
    non-coding exons contribute an empty (0,0) placeholder so the index == rank-1.
    """
    ordered = sorted(exons, key=lambda e: e["start"], reverse=(strand == -1))
    out, pos = [], 0
    for e in ordered:
        o_lo, o_hi = max(e["start"], cds_g_start), min(e["end"], cds_g_end)
        clen = max(0, o_hi - o_lo + 1)
        if clen:
            out.append((pos + 1, pos + clen)); pos += clen
        else:
            out.append((0, 0))
    return out


def cds_coord_at_exon_boundary(exon_cds: list[tuple[int, int]], exon_rank: int,
                               side: Literal["end", "start"]) -> int:
    """CDS coordinate at the 3' end (side='end') or 5' start (side='start') of a 1-based exon rank."""
    s, e = exon_cds[exon_rank - 1]
    if (s, e) == (0, 0):
        raise ValueError(f"exon {exon_rank} is non-coding")
    return e if side == "end" else s


# ----------------------------------------------------------------------------
# Layer 2: the HGVS.p-like interface object.
# ----------------------------------------------------------------------------
@dataclass
class DomainCall:
    accession: str
    name: str
    type: str
    start: int
    end: int
    status: Literal["RETAINED", "LOST", "DISRUPTED"]

@dataclass
class FusionProtein:
    """Normalized protein-level representation of a fusion — the interface layer.

    This is the 'hgvs.p-like' contract: it fully describes the chimeric protein
    product and is what the knowledge layer keys against (via .categorical_key()).
    """
    five_gene: str
    three_gene: str
    five_transcript: str
    three_transcript: str
    # protein breakpoints (1-based): last aa fully from 5' partner, first aa from 3' partner
    five_last_aa: int
    three_first_aa: int
    five_last_aa_res: str           # residue (1-letter) at five_last_aa in the 5' protein
    three_first_aa_res: str         # residue at three_first_aa in the 3' protein
    in_frame: bool
    hybrid_codon: bool              # True if a codon straddles the junction
    junction_residue: Optional[str] # residue encoded by a hybrid codon, else None
    fusion_length: int
    internal_stops: int
    fusion_protein_seq: str
    frame_status: Literal["in-frame", "out-of-frame", "frameshift-truncating"]
    domains: list[DomainCall] = field(default_factory=list)
    breakpoint_label: Optional[str] = None   # e.g. "E13;A20"

    def categorical_key(self) -> str:
        """The OncoKB/CIViC-style categorical fusion key, gene-pair only."""
        return f"{self.five_gene}::{self.three_gene}"

    def to_hgvsp(self) -> str:
        """HGVS.p-like junction string using the '::' adjoined-protein operator.

        Format:  FIVE:p.Met1_<Xaa><n>::THREE:p.<Xaa><m>_<Xaa><end>
        A hybrid junction codon is annotated in a trailing comment.
        """
        five = f"{self.five_gene}:p.Met1_{aa3(self.five_last_aa_res)}{self.five_last_aa}"
        end_res = self.fusion_protein_seq[-1] if self.fusion_protein_seq else 'X'
        three = (f"{self.three_gene}:p.{aa3(self.three_first_aa_res)}{self.three_first_aa}_"
                 f"{aa3(end_res)}{self._three_last_aa()}")
        s = f"{five}::{three}"
        if self.hybrid_codon and self.junction_residue:
            s += f"  (junction hybrid codon -> {aa3(self.junction_residue)})"
        return s

    def _three_last_aa(self) -> int:
        # length of 3' contribution = fusion_length - five_last_aa (- hybrid residue if any)
        c = self.fusion_length - self.five_last_aa - (1 if self.hybrid_codon else 0)
        return self.three_first_aa + c - 1

    def to_dict(self) -> dict:
        d = asdict(self)
        d["categorical_key"] = self.categorical_key()
        d["hgvsp_like"] = self.to_hgvsp()
        return d


# ----------------------------------------------------------------------------
# Data provider protocol (Layers pull real annotation through this).
# ----------------------------------------------------------------------------
class DataProvider(Protocol):
    def get_transcript(self, gene_or_tx: str) -> Transcript: ...
    def get_domains(self, uniprot: str) -> list[dict]: ...            # [{accession,name,type,start,end}]
    def get_fusion_knowledge(self, categorical_key: str) -> dict: ...  # OncoKB-like


# ----------------------------------------------------------------------------
# Layer 1: effect predictor (VEP-like).
# ----------------------------------------------------------------------------
def annotate_effect(five: Transcript, three: Transcript,
                    five_cds_end: int, three_cds_start: int,
                    breakpoint_label: Optional[str] = None,
                    domains_five: Optional[list[dict]] = None,
                    domains_three: Optional[list[dict]] = None) -> FusionProtein:
    """Reconstruct and characterize the chimeric protein.

    five_cds_end     : 1-based CDS coord (inclusive) of the last 5'-partner base kept.
    three_cds_start  : 1-based CDS coord of the first 3'-partner base kept.
    """
    fusion_cds = five.cds[:five_cds_end] + three.cds[three_cds_start - 1:]
    prot_full = translate(fusion_cds)
    core = prot_full.split('*')[0]                 # up to first stop
    internal_stops = core.count('*')
    in_frame = (len(fusion_cds) % 3 == 0)

    five_complete_codons = five_cds_end // 3
    remainder = five_cds_end % 3
    hybrid = remainder != 0
    five_last_aa = five_complete_codons
    junction_res = core[five_last_aa] if hybrid and len(core) > five_last_aa else None
    # Exact 3'-first-residue via phase: a hybrid codon consumes (3-remainder) 3'
    # nucleotides, so the first *fully-retained* 3' residue starts after them.
    alk_nt_in_hybrid = (3 - remainder) if hybrid else 0
    three_first_full_cds = three_cds_start + alk_nt_in_hybrid
    three_first_aa = (three_first_full_cds - 1) // 3 + 1  # first full 3'-partner residue

    # frame status
    downstream_len = len(prot_full) - five_complete_codons
    if not in_frame:
        status = "frameshift-truncating"
    elif internal_stops > 0:
        status = "out-of-frame"
    else:
        status = "in-frame"

    # domain retention
    dcalls: list[DomainCall] = []
    for d in (domains_five or []):
        if d["end"] <= five_last_aa:              st = "RETAINED"
        elif d["start"] > five_last_aa:           st = "LOST"
        else:                                     st = "DISRUPTED"
        dcalls.append(DomainCall(d["accession"], d["name"], d["type"], d["start"], d["end"], st))
    for d in (domains_three or []):
        if d["start"] >= three_first_aa:          st = "RETAINED"
        elif d["end"] < three_first_aa:           st = "LOST"
        else:                                     st = "DISRUPTED"
        dcalls.append(DomainCall(d["accession"], d["name"], d["type"], d["start"], d["end"], st))

    return FusionProtein(
        five_gene=five.gene_symbol, three_gene=three.gene_symbol,
        five_transcript=five.transcript_id, three_transcript=three.transcript_id,
        five_last_aa=five_last_aa, three_first_aa=three_first_aa,
        five_last_aa_res=five.protein[five_last_aa - 1] if five_last_aa <= len(five.protein) else 'X',
        three_first_aa_res=three.protein[three_first_aa - 1] if three_first_aa <= len(three.protein) else 'X',
        in_frame=in_frame, hybrid_codon=hybrid, junction_residue=junction_res,
        fusion_length=len(core), internal_stops=internal_stops,
        fusion_protein_seq=core, frame_status=status,
        domains=dcalls, breakpoint_label=breakpoint_label)


# ----------------------------------------------------------------------------
# Layer 3: knowledge annotator (OncoKB-like), keyed on categorical fusion.
# ----------------------------------------------------------------------------
@dataclass
class FusionKnowledge:
    categorical_key: str
    oncogenic: Optional[str] = None
    therapies: list[str] = field(default_factory=list)
    evidence: list[dict] = field(default_factory=list)
    diseases: list[str] = field(default_factory=list)
    sources: list[str] = field(default_factory=list)

def annotate_knowledge(fp: FusionProtein, provider: DataProvider) -> FusionKnowledge:
    raw = provider.get_fusion_knowledge(fp.categorical_key())
    return FusionKnowledge(categorical_key=fp.categorical_key(), **raw)


# ----------------------------------------------------------------------------
# Orchestration: the full three-layer annotate() call.
# ----------------------------------------------------------------------------
def annotate_fusion(provider: DataProvider,
                    five_gene: str, three_gene: str,
                    five_exon: int, three_exon: int,
                    five_tx: Optional[str] = None, three_tx: Optional[str] = None
                    ) -> dict:
    """End-to-end: resolve transcripts -> effect -> interface -> knowledge."""
    five = provider.get_transcript(five_tx or five_gene)
    three = provider.get_transcript(three_tx or three_gene)
    five_cds_end = cds_coord_at_exon_boundary(five.exon_cds, five_exon, "end")
    three_cds_start = cds_coord_at_exon_boundary(three.exon_cds, three_exon, "start")
    dfive = provider.get_domains(five.uniprot) if five.uniprot else []
    dthree = provider.get_domains(three.uniprot) if three.uniprot else []
    fp = annotate_effect(five, three, five_cds_end, three_cds_start,
                         breakpoint_label=f"E{five_exon};A{three_exon}",
                         domains_five=dfive, domains_three=dthree)
    kn = annotate_knowledge(fp, provider)
    return {"interface": fp.to_dict(), "knowledge": asdict(kn)}
