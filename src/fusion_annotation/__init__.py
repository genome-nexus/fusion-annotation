"""fusion_annotation — a standards-aligned gene-fusion annotation system.

Two engines joined by an HGVS.p-like protein-level interface:

  * EFFECT engine (VEP-like)    -> annotate_effect / annotate_fusion
  * INTERFACE (HGVS.p-like)     -> FusionProtein (::-junction, VICC GFS model)
  * KNOWLEDGE engine (OncoKB-like) -> annotate_knowledge

See docs/DESIGN.md for the architecture and the EML4::ALK worked example.
"""
from .core import (
    Transcript,
    FusionProtein,
    DomainCall,
    FusionKnowledge,
    DataProvider,
    build_exon_cds_map,
    build_exon_genomic_map,
    cds_coord_at_exon_boundary,
    cds_coord_at_genomic,
    parse_genomic_breakpoint,
    annotate_effect,
    annotate_knowledge,
    annotate_fusion,
    KNOWN_ONCOGENIC_PAIRS,
    translate,
    aa3,
)

__all__ = [
    "Transcript", "FusionProtein", "DomainCall", "FusionKnowledge",
    "DataProvider", "build_exon_cds_map", "build_exon_genomic_map",
    "cds_coord_at_exon_boundary", "cds_coord_at_genomic", "parse_genomic_breakpoint",
    "annotate_effect", "annotate_knowledge", "annotate_fusion",
    "KNOWN_ONCOGENIC_PAIRS", "translate", "aa3",
]
__version__ = "0.4.0"
