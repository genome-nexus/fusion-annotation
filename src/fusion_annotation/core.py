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
from concurrent.futures import ThreadPoolExecutor
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
    # Genomic structure, needed to map a genomic breakpoint -> CDS coord. Optional so
    # exon-only callers/fixtures keep working; populated by build_exon_genomic_map().
    # exon_genomic[i] = (g_lo, g_hi) for exon i, same rank order/index as exon_cds.
    exon_genomic: list[tuple[int, int]] = field(default_factory=list)
    cds_g_start: int = 0     # genomic CDS bounds (lo < hi), regardless of strand
    cds_g_end: int = 0
    is_canonical: Optional[bool] = None   # True if this is the gene's canonical/MANE transcript

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


def build_exon_genomic_map(strand: int, exons: list[dict]) -> list[tuple[int, int]]:
    """Genomic (g_lo, g_hi) bounds per exon, in transcription order.

    Same ordering and indexing (index == rank-1) as build_exon_cds_map(), so the
    two lists line up. Kept separate from exon_cds so a genomic breakpoint can be
    resolved to a CDS coordinate without re-fetching the transcript structure.
    """
    ordered = sorted(exons, key=lambda e: e["start"], reverse=(strand == -1))
    return [(e["start"], e["end"]) for e in ordered]


def cds_coord_at_genomic(strand: int, exon_genomic: list[tuple[int, int]],
                         cds_g_start: int, cds_g_end: int,
                         g_pos: int, side: Literal["end", "start"]) -> int:
    """Map a genomic breakpoint to a 1-based CDS coordinate on this transcript.

    Unlike an exon number, a genomic coordinate pins the isoform: it only resolves
    against the transcript whose exon table actually spans it, which is exactly what
    disambiguates overlapping isoforms with divergent exon numbering (see issue #3).

    Parameters
    ----------
    exon_genomic : (g_lo, g_hi) per exon in transcription order (build_exon_genomic_map()).
    cds_g_start / cds_g_end : genomic CDS bounds (lo < hi); UTR-only exon portions are ignored.
    g_pos : genomic position of the breakpoint.
    side : 'end'   -> last CDS base the 5' partner retains at/upstream of g_pos.
           'start' -> first CDS base the 3' partner retains at/downstream of g_pos.

    A breakpoint inside a coding exon maps to that exact base; a breakpoint in an
    intron/UTR snaps to the nearest flanking coding-exon boundary in the relevant
    direction. Raises ValueError if g_pos maps outside the coding region entirely.
    """
    if not exon_genomic:
        raise ValueError("transcript has no genomic exon map; cannot resolve a genomic breakpoint "
                         "(populate Transcript.exon_genomic via build_exon_genomic_map)")

    def tc(g: int) -> int:                      # transcription coordinate (increases 5'->3')
        return g if strand == 1 else -g

    tp = tc(g_pos)
    pos = 0
    upstream_cds_end: Optional[int] = None      # cds_end of the last coding exon fully 5' of g_pos
    for g_lo, g_hi in exon_genomic:
        o_lo, o_hi = max(g_lo, cds_g_start), min(g_hi, cds_g_end)
        clen = o_hi - o_lo + 1
        if clen <= 0:                           # non-coding / UTR-only exon
            continue
        c_start, c_end = pos + 1, pos + clen
        pos += clen
        t_lo, t_hi = sorted((tc(o_lo), tc(o_hi)))
        if t_lo <= tp <= t_hi:                  # breakpoint inside this coding exon
            return c_start + (tp - t_lo)        # exact base (c_start is the 5' end of the exon)
        if tp > t_hi:                           # exon lies entirely upstream (5') of g_pos
            upstream_cds_end = c_end
        elif side == "start":                   # first exon entirely downstream (3') of g_pos
            return c_start
    if side == "end" and upstream_cds_end is not None:
        return upstream_cds_end
    raise ValueError(f"genomic position {g_pos} does not map into the CDS of this transcript")


def parse_genomic_breakpoint(value) -> int:
    """Parse a genomic breakpoint into an integer position.

    Accepts a plain int, a bare position string ("117324415"), a "chrom:pos" form
    ("chr6:117324415" / "6:117324415"), or an HGVS genomic term ("g.117324415",
    "chr6:g.117324415"). Only the 1-based genomic position is used — the transcript's
    own strand/exon table supplies chromosome context — so the chromosome token, if
    present, is accepted but not validated here.
    """
    if isinstance(value, int):
        return value
    s = str(value).strip()
    if ":" in s:
        s = s.rsplit(":", 1)[1].strip()
    if s.lower().startswith("g."):
        s = s[2:].strip()
    s = s.replace(",", "").replace("_", "")
    try:
        return int(s)
    except ValueError as exc:
        raise ValueError(f"could not parse genomic breakpoint {value!r}") from exc


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
    # Which fusion partner this domain call belongs to (the 5' or 3' gene
    # symbol), and that partner's full-length protein size — both needed by
    # UI consumers (e.g. web/) to lay out a two-track domain-retention
    # diagram without re-deriving partner attribution from list order.
    gene: str = ""
    partner_protein_length: int = 0

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
        dcalls.append(DomainCall(d["accession"], d["name"], d["type"], d["start"], d["end"], st,
                                 gene=five.gene_symbol, partner_protein_length=len(five.protein)))
    for d in (domains_three or []):
        if d["start"] >= three_first_aa:          st = "RETAINED"
        elif d["end"] < three_first_aa:           st = "LOST"
        else:                                     st = "DISRUPTED"
        dcalls.append(DomainCall(d["accession"], d["name"], d["type"], d["start"], d["end"], st,
                                 gene=three.gene_symbol, partner_protein_length=len(three.protein)))

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


# Curated gene pairs that are established oncogenic fusion drivers. Used only to
# raise a sanity flag when such a pair reconstructs out-of-frame — an in-frame
# driver that comes back frameshifted almost always means a wrong exon number or
# transcript isoform was supplied (issue #3), not a real frameshift. Symbols are
# order-independent (either partner may be 5' or 3').
KNOWN_ONCOGENIC_PAIRS: frozenset[frozenset[str]] = frozenset({
    frozenset({"EML4", "ALK"}), frozenset({"NPM1", "ALK"}),
    frozenset({"LMNA", "NTRK1"}), frozenset({"CD74", "ROS1"}),
    frozenset({"TPM3", "NTRK1"}), frozenset({"ETV6", "NTRK3"}),
    frozenset({"KIF5B", "RET"}), frozenset({"CCDC6", "RET"}),
    frozenset({"SLC34A2", "ROS1"}), frozenset({"BCR", "ABL1"}),
    frozenset({"FGFR3", "TACC3"}), frozenset({"PAX8", "PPARG"}),
})


def _resolve_breakpoint(tx: Transcript, side: Literal["end", "start"],
                        exon: Optional[int], genomic) -> tuple[int, dict]:
    """Resolve a partner's breakpoint to a CDS coordinate, preferring genomic input.

    Returns (cds_coord, provenance) where provenance records how it was derived so
    the caller can echo it back to the user.
    """
    if genomic is not None:
        g_pos = parse_genomic_breakpoint(genomic)
        cds = cds_coord_at_genomic(tx.strand, tx.exon_genomic, tx.cds_g_start,
                                   tx.cds_g_end, g_pos, side)
        return cds, {"type": "genomic", "genomic_position": g_pos, "cds_coord": cds}
    if exon is not None:
        cds = cds_coord_at_exon_boundary(tx.exon_cds, exon, side)
        return cds, {"type": "exon", "exon": exon, "cds_coord": cds}
    raise ValueError(f"no breakpoint given for the {'5' if side == 'end' else '3'}' partner: "
                     "supply an exon number or a genomic position")


# ----------------------------------------------------------------------------
# Orchestration: the full three-layer annotate() call.
# ----------------------------------------------------------------------------
def annotate_fusion(provider: DataProvider,
                    five_gene: str, three_gene: str,
                    five_exon: Optional[int] = None, three_exon: Optional[int] = None,
                    five_tx: Optional[str] = None, three_tx: Optional[str] = None,
                    five_genomic=None, three_genomic=None,
                    ) -> dict:
    """End-to-end: resolve transcripts -> effect -> interface -> knowledge.

    Each partner's breakpoint may be given as an exon number (``five_exon`` /
    ``three_exon``) or as a genomic position (``five_genomic`` / ``three_genomic``:
    int, "chr6:117324415", or HGVS "g." form). A genomic breakpoint pins the isoform
    and is preferred when both are supplied — this is what resolves the isoform
    ambiguity described in issue #3. The resolved transcript ids and how each
    breakpoint was interpreted are echoed back under ``resolved``; a ``warnings``
    list flags a known oncogenic pair that reconstructs out-of-frame.
    """
    # Parallelise the two independent partner lookups (each may involve several
    # network round-trips; doing them concurrently roughly halves this phase).
    with ThreadPoolExecutor(max_workers=2) as pool:
        f_five = pool.submit(provider.get_transcript, five_tx or five_gene)
        f_three = pool.submit(provider.get_transcript, three_tx or three_gene)
        five = f_five.result()
        three = f_three.result()

    five_cds_end, five_bp = _resolve_breakpoint(five, "end", five_exon, five_genomic)
    three_cds_start, three_bp = _resolve_breakpoint(three, "start", three_exon, three_genomic)

    dfive = provider.get_domains(five.uniprot) if five.uniprot else []
    dthree = provider.get_domains(three.uniprot) if three.uniprot else []

    label_5 = f"g.{five_bp['genomic_position']}" if five_bp["type"] == "genomic" else f"E{five_exon}"
    label_3 = f"g.{three_bp['genomic_position']}" if three_bp["type"] == "genomic" else f"A{three_exon}"
    fp = annotate_effect(five, three, five_cds_end, three_cds_start,
                         breakpoint_label=f"{label_5};{label_3}",
                         domains_five=dfive, domains_three=dthree)
    kn = annotate_knowledge(fp, provider)

    resolved = {
        # Echo the genome build so a caller can confirm their genomic breakpoints
        # were interpreted against the build they intended (defaults to GRCh38).
        "genome_build": getattr(provider, "assembly", "GRCh38"),
        "five": _echo_partner(five, five_tx, five_bp),
        "three": _echo_partner(three, three_tx, three_bp),
    }
    both_genomic = five_bp["type"] == "genomic" and three_bp["type"] == "genomic"
    warnings = _sanity_warnings(five_gene, three_gene, fp, kn, both_genomic)
    # Surface any domain-degradation warning set by the provider (e.g. InterPro
    # throttle, GN Pfam unavailable).
    domain_warn = getattr(provider, "_domain_warning", None)
    if domain_warn:
        warnings.append(domain_warn)

    return {"interface": fp.to_dict(), "knowledge": asdict(kn),
            "resolved": resolved, "warnings": warnings}


def _echo_partner(tx: Transcript, user_tx: Optional[str], breakpoint_prov: dict) -> dict:
    """Echo the transcript actually used and how it was selected, for caller transparency."""
    if user_tx:
        source = "user-specified"
    elif tx.is_canonical:
        source = "canonical"
    elif tx.is_canonical is False:
        source = "non-canonical"
    else:
        source = "provider-default"
    return {
        "gene": tx.gene_symbol,
        "transcript": tx.transcript_id,
        "transcript_source": source,
        "breakpoint": breakpoint_prov,
    }


def _sanity_warnings(five_gene: str, three_gene: str,
                     fp: FusionProtein, kn: "FusionKnowledge",
                     both_genomic: bool) -> list[str]:
    """Flag a known oncogenic partner pair that reconstructs out-of-frame (issue #3)."""
    warnings: list[str] = []
    pair = frozenset({five_gene.upper(), three_gene.upper()})
    known = pair in KNOWN_ONCOGENIC_PAIRS or bool(kn.oncogenic)
    if known and fp.frame_status != "in-frame":
        # Genomic breakpoints already pin the isoform, so don't suggest them again.
        remedy = ("Double-check the genomic breakpoints and the resolved transcripts "
                  "against the assayed fusion." if both_genomic else
                  "Re-check the exon numbers against the pinned transcripts, or supply "
                  "genomic breakpoints (five_genomic/three_genomic) to remove "
                  "exon-numbering ambiguity between isoforms.")
        warnings.append(
            f"{fp.categorical_key()} is a known oncogenic fusion pair but this breakpoint "
            f"reconstructs as '{fp.frame_status}'. This often means the breakpoints or "
            f"transcript isoforms do not match the assayed fusion. {remedy}")
    return warnings
