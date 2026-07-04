// Pure domain-diagram logic shared by DomainDiagram.tsx: domain de-duplication,
// category-based coloring, and label-collision avoidance. Split out from the
// component file so Vite's fast-refresh only sees component exports there.
import type { DomainCall } from "./types";

export type Status = "RETAINED" | "LOST" | "DISRUPTED";

export interface CanonDomain {
  name: string;
  start: number;
  end: number;
  status: Status;
}

// ---------------------------------------------------------------------------
// Domain -> color: a single lookup used for every track (5' partner, 3'
// partner, and the fusion protein), so the *same* domain always renders in
// the same color wherever it appears. This mirrors docs/generate_domain_map.py,
// which regenerates docs/fusion_domain_map.png with the same scheme — fixing
// https://github.com/genome-nexus/fusion-annotation/issues/17, where the ALK
// kinase domain used to render orange in one track and blue in another.
// Curated colors cover the domain categories seen in common fusions; any
// other domain name still gets a stable (hash-based) color from the
// fallback palette, so it's consistent across renders without being random.
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

// ---------------------------------------------------------------------------
// Domain de-duplication: Genome Nexus reports many overlapping InterPro/Pfam
// records describing the same physical domain (different source databases,
// slightly different boundaries — e.g. 7+ records covering the ALK kinase
// region alone). Rendering every raw record as its own rectangle is
// illegible, so overlapping same-type records are merged into one
// representative block, preferring more granular types (repeat > domain >
// conserved_site) and dropping coarser records once a finer one already
// covers the same range. Same algorithm as docs/generate_domain_map.py.
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
  const items = raw.filter((d) => d.gene === gene && KEEP_TYPES.has(d.type) && d.name !== d.accession);
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
    const name = c.reduce((shortest, m) => (m.name.length < shortest.length ? m.name : shortest), c[0].name);
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
    const halfWidth = g.name.length * proteinLength * 0.0055;
    let row = 0;
    while (row < rowFreeAt.length && rowFreeAt[row] > center - halfWidth) row++;
    if (row === rowFreeAt.length) rowFreeAt.push(center + halfWidth);
    else rowFreeAt[row] = center + halfWidth;
    return { center, row, name: g.name };
  });
}
