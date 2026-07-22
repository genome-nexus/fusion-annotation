/**
 * Compute equivalent breakpoint representations from an annotation result.
 *
 * After annotating with any input mode the engine resolves a full transcript
 * model for each partner. This module extracts all four representations —
 * genomic chr:pos, exon rank, CDS nucleotide offset, and protein amino-acid
 * position — so the UI can show them as read-only derived values and pre-fill
 * the form when the user switches input modes.
 */
import type { AnnotationResult, AnnotateParams, ResolvedPartner } from "./types";

export interface PartnerDerived {
  gene: string;
  transcript: string;
  /** Genomic chr:pos, e.g. "chr2:42295516" — from input or derived from exon boundary */
  genomic: string | null;
  /** 1-based exon rank (last exon for 5′ partner, first exon for 3′ partner) */
  exon: number | null;
  /** 1-based CDS nucleotide coordinate */
  cds: number | null;
  /** Amino-acid position (last for 5′, first for 3′) */
  protein: number | null;
}

export interface DerivedInputs {
  five: PartnerDerived;
  three: PartnerDerived;
  genomeBuild: string;
}

/** Derive the exon-boundary genomic position when the input was exon-based.
 *
 * For the 5′ partner we want the *end* of the last exon it contributes;
 * for the 3′ partner we want the *start* of the first exon it contributes.
 * Genomic start/end are always lo < hi regardless of strand; for a minus-
 * strand gene the "end" in transcript space is the lower genomic coordinate.
 */
function exonBoundaryGenomic(
  partner: ResolvedPartner,
  side: "five" | "three",
): string | null {
  const { structure, breakpoint } = partner;
  if (!structure || !structure.chrom) return null;
  const rank = breakpoint.context.exon_rank;
  if (rank == null) return null;
  const exon = structure.exons[rank - 1];
  if (!exon) return null;

  let pos: number;
  if (side === "five") {
    // Breakpoint is at the end of the last 5′ exon (in transcript direction)
    pos = structure.strand === 1 ? exon.genomic_end : exon.genomic_start;
  } else {
    // Breakpoint is at the start of the first 3′ exon (in transcript direction)
    pos = structure.strand === 1 ? exon.genomic_start : exon.genomic_end;
  }
  return `${structure.chrom}:${pos}`;
}

function derivePartner(
  partner: ResolvedPartner,
  proteinAa: number | null,
  side: "five" | "three",
): PartnerDerived {
  const { breakpoint, gene, transcript } = partner;

  // Genomic: use input value if provided, otherwise derive from exon boundary
  let genomic: string | null = null;
  if (breakpoint.genomic_position != null && partner.structure?.chrom) {
    genomic = `${partner.structure.chrom}:${breakpoint.genomic_position}`;
  } else {
    genomic = exonBoundaryGenomic(partner, side);
  }

  return {
    gene,
    transcript,
    genomic,
    exon: breakpoint.context.exon_rank ?? null,
    cds: breakpoint.cds_coord ?? null,
    protein: proteinAa,
  };
}

export function computeDerivedInputs(result: AnnotationResult): DerivedInputs {
  const { interface: iface, resolved } = result;
  return {
    five: derivePartner(resolved.five, iface.five_last_aa, "five"),
    three: derivePartner(resolved.three, iface.three_first_aa, "three"),
    genomeBuild: resolved.genome_build,
  };
}

/** Convert derived inputs back into AnnotateParams for a specific mode,
 * so the form can pre-fill when the user switches modes. */
export function derivedToParams(
  derived: DerivedInputs,
  mode: "exon" | "genomic",
  current: AnnotateParams,
): AnnotateParams {
  const base: AnnotateParams = {
    ...current,
    five_gene: derived.five.gene,
    three_gene: derived.three.gene,
    genome_build: derived.genomeBuild,
    input_mode: mode,
    // keep transcripts if already set
    five_transcript: current.five_transcript || derived.five.transcript,
    three_transcript: current.three_transcript || derived.three.transcript,
  };

  if (mode === "exon") {
    return {
      ...base,
      five_exon: derived.five.exon != null ? String(derived.five.exon) : "",
      three_exon: derived.three.exon != null ? String(derived.three.exon) : "",
      five_genomic: "",
      three_genomic: "",
    };
  }
  // genomic mode
  return {
    ...base,
    five_genomic: derived.five.genomic ?? "",
    three_genomic: derived.three.genomic ?? "",
    five_exon: "",
    three_exon: "",
  };
}
