// Pure domain-diagram logic shared by DomainDiagram.tsx: domain de-duplication,
// category-based coloring, structure-track layout, and label-collision
// avoidance. Split out from the component file so Vite's fast-refresh only
// sees component exports there.
import type { DomainCall, ResolvedPartner, TranscriptStructure } from "./types";

export type Status = "RETAINED" | "LOST" | "DISRUPTED" | "UNKNOWN";

export interface CanonDomain {
  name: string;
  start: number;
  end: number;
  status: Status;
}

export interface StructureLayoutExon {
  rank: number;
  start: number;
  end: number;
  width: number;
  length: number;
  segments: TranscriptStructure["exons"][number]["segments"];
}

export interface StructureLayout {
  promoterStart: number;
  promoterEnd: number;
  width: number;
  exons: StructureLayoutExon[];
}

// ---------------------------------------------------------------------------
// Domain -> color: a single lookup used for every track (5' partner, 3'
// partner, and the fusion protein), so the *same* domain always renders in
// the same color wherever it appears.
// ---------------------------------------------------------------------------
const CATEGORY_KEYWORDS: Array<[string, string]> = [
  ["kinase", "#e8590c"],
  ["wd40", "#2f9e44"],
  ["beta-propeller", "#2f9e44"],
  ["help", "#0c8599"],
  ["mam domain", "#1971c2"],
];
const FALLBACK_PALETTE = [
  "#495057", "#c2255c", "#5f3dc4", "#0b7285", "#e67700", "#1864ab", "#862e9c",
];
const PROMOTER_WIDTH = 34;
const EXON_GAP = 12;
const EXON_MIN_WIDTH = 18;
const EXON_MAX_WIDTH = 34;
const LABEL_HALF_WIDTH_SCALE = 0.007;
const EDGE_LABEL_CHAR_WIDTH = 4.2;

function hashString(s: string): number {
  let h = 0;
  for (let i = 0; i < s.length; i++) {
    h = (Math.imul(31, h) + s.charCodeAt(i)) | 0;
  }
  return Math.abs(h);
}

export function colorFor(name: string): string {
  const lower = name.toLowerCase();
  for (const [keyword, color] of CATEGORY_KEYWORDS) {
    if (lower.includes(keyword)) return color;
  }
  return FALLBACK_PALETTE[hashString(lower) % FALLBACK_PALETTE.length];
}

export function structureSegmentColor(kind: "utr5" | "coding" | "utr3"): string {
  if (kind === "coding") return "#4c6ef5";
  if (kind === "utr5") return "#a5d8ff";
  return "#d0ebff";
}

export function edgeAwareTextPlacement(
  x: number,
  text: string,
  minX: number,
  maxX: number,
): { x: number; anchor: "start" | "middle" | "end" } {
  const halfWidth = (text.length * EDGE_LABEL_CHAR_WIDTH) / 2;
  if (x - halfWidth < minX) return { x: minX + 2, anchor: "start" };
  if (x + halfWidth > maxX) return { x: maxX - 2, anchor: "end" };
  return { x, anchor: "middle" };
}

// ---------------------------------------------------------------------------
// Domain de-duplication: Genome Nexus reports many overlapping InterPro/Pfam
// records describing the same physical domain (different source databases,
// slightly different boundaries). Rendering every raw record as its own
// rectangle is illegible, so overlapping same-type records are merged into one
// representative block, preferring more granular types (repeat > domain >
// conserved_site) and dropping coarser records once a finer one already
// covers the same range. If a cluster has both curated labels and accession-
// style fallback names, prefer the curated label; otherwise keep the
// accession so Pfam-only regions still render. Same algorithm as
// docs/generate_domain_map.py.
// ---------------------------------------------------------------------------
const KEEP_TYPES = new Set(["domain", "repeat", "conserved_site"]);
const TYPE_PRIORITY = ["repeat", "domain", "conserved_site"];

interface Span {
  start: number;
  end: number;
}

function overlaps(a: Span, b: Span): boolean {
  return !(a.end < b.start || a.start > b.end);
}

function clusterByOverlap(items: DomainCall[]): DomainCall[][] {
  const clusters: DomainCall[][] = [];
  for (const it of items) {
    let merged = [it];
    const remaining: DomainCall[][] = [];
    for (const c of clusters) {
      if (c.some((m) => overlaps(it, m))) {
        merged = merged.concat(c);
      } else {
        remaining.push(c);
      }
    }
    remaining.push(merged);
    clusters.splice(0, clusters.length, ...remaining);
  }
  return clusters;
}

function spanOverlapFrac(span: Span, clusters: DomainCall[][]): number {
  const total = span.end - span.start + 1;
  let covered = 0;
  for (const c of clusters) {
    const cs = Math.min(...c.map((m) => m.start));
    const ce = Math.max(...c.map((m) => m.end));
    const os = Math.max(span.start, cs);
    const oe = Math.min(span.end, ce);
    if (oe >= os) covered += oe - os + 1;
  }
  return total > 0 ? covered / total : 0;
}

export function canonicalizeDomains(raw: DomainCall[], gene: string): CanonDomain[] {
  const items = raw.filter((d) => d.gene === gene && KEEP_TYPES.has(d.type));
  const kept: DomainCall[][] = [];
  for (const type of TYPE_PRIORITY) {
    const clusters = clusterByOverlap(items.filter((d) => d.type === type));
    for (const c of clusters) {
      const span = { start: Math.min(...c.map((m) => m.start)), end: Math.max(...c.map((m) => m.end)) };
      if (kept.length > 0 && spanOverlapFrac(span, kept) >= 0.5) continue;
      kept.push(c);
    }
  }
  const reps = kept.map((c) => {
    const start = Math.min(...c.map((m) => m.start));
    const end = Math.max(...c.map((m) => m.end));
    const preferred = c.filter((m) => m.name !== m.accession);
    const targets = preferred.length > 0 ? preferred : c;
    const name = targets.reduce(
      (shortest, m) => (m.name.length < shortest.length ? m.name : shortest),
      targets[0].name,
    );
    const statuses = new Set(c.map((m) => m.status));
    const status: Status =
      statuses.size === 1 && statuses.has("RETAINED")
        ? "RETAINED"
        : statuses.size === 1 && statuses.has("LOST")
          ? "LOST"
          : "DISRUPTED";
    return { name, start, end, status };
  });
  reps.sort((a, b) => a.start - b.start);
  return reps;
}

function exonWidth(length: number): number {
  return Math.max(EXON_MIN_WIDTH, Math.min(EXON_MAX_WIDTH, 10 + Math.sqrt(Math.max(length, 1)) * 0.8));
}

export function layoutTranscriptStructure(structure: TranscriptStructure): StructureLayout {
  const exons: StructureLayoutExon[] = [];
  let cursor = PROMOTER_WIDTH + EXON_GAP;
  for (const exon of structure.exons) {
    const width = exonWidth(exon.length);
    exons.push({
      rank: exon.rank,
      start: cursor,
      end: cursor + width,
      width,
      length: exon.length,
      segments: exon.segments,
    });
    cursor += width + EXON_GAP;
  }
  return {
    promoterStart: 0,
    promoterEnd: PROMOTER_WIDTH,
    width: exons.length > 0 ? exons[exons.length - 1].end + EXON_GAP / 2 : PROMOTER_WIDTH + EXON_GAP,
    exons,
  };
}

export function transcriptBreakpointPosition(partner: ResolvedPartner, layout: StructureLayout): number | null {
  const { context } = partner.breakpoint;
  if (context.region === "upstream") return (layout.promoterStart + layout.promoterEnd) / 2;
  if (context.region === "downstream") return layout.width - EXON_GAP / 4;
  if (context.intron_rank != null) {
    const left = layout.exons[context.intron_rank - 1];
    const right = layout.exons[context.intron_rank];
    if (!left || !right) return null;
    return (left.end + right.start) / 2;
  }
  if (context.exon_rank != null) {
    const exon = layout.exons[context.exon_rank - 1];
    if (!exon) return null;
    if (context.boundary === "before") return exon.start;
    if (context.boundary === "after") return exon.end;
    if (context.exon_offset != null && context.exon_length != null) {
      const frac = context.exon_length <= 1 ? 0.5 : (context.exon_offset - 1) / (context.exon_length - 1);
      return exon.start + frac * exon.width;
    }
    return (exon.start + exon.end) / 2;
  }
  return null;
}

export function transcriptBreakpointLabel(partner: ResolvedPartner): string {
  const { breakpoint } = partner;
  if (breakpoint.type === "genomic" && breakpoint.genomic_position != null) {
    return `g.${breakpoint.genomic_position} · ${breakpoint.context.label}`;
  }
  return breakpoint.context.label;
}

/** Group adjacent same-name domains under one shared label, then stagger
 * labels into extra rows when their estimated text width would collide with
 * a neighboring label. */
export function labelRows(domains: CanonDomain[], proteinLength: number) {
  const groups: { name: string; start: number; end: number }[] = [];
  const gapTolerance = proteinLength * 0.03;
  for (const d of domains) {
    const last = groups[groups.length - 1];
    if (last && last.name === d.name && d.start - last.end <= gapTolerance) {
      last.end = Math.max(last.end, d.end);
    } else {
      groups.push({ name: d.name, start: d.start, end: d.end });
    }
  }
  const rowFreeAt: number[] = [];
  return groups.map((g) => {
    const center = (g.start + g.end) / 2;
    const halfWidth = g.name.length * proteinLength * LABEL_HALF_WIDTH_SCALE;
    let row = 0;
    while (row < rowFreeAt.length && rowFreeAt[row] > center - halfWidth) row++;
    if (row === rowFreeAt.length) rowFreeAt.push(center + halfWidth);
    else rowFreeAt[row] = center + halfWidth;
    return { center, row, name: g.name };
  });
}
