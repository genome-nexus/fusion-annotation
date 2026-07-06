import { useState } from "react";
import {
  canonicalizeDomains,
  colorFor,
  labelRows,
  layoutTranscriptStructure,
  structureSegmentColor,
  transcriptBreakpointLabel,
  transcriptBreakpointPosition,
  type CanonDomain,
  type Status,
} from "../lib/domainDiagram";
import type { DomainCall, ResolvedPartner, TranscriptStructureSegment } from "../lib/types";

interface Props {
  domains: DomainCall[];
  fiveGene: string;
  threeGene: string;
  fivePartner: ResolvedPartner;
  threePartner: ResolvedPartner;
  fiveLastAa: number;
  threeFirstAa: number;
  hybridCodon: boolean;
  fusionLength: number;
}

const TRACK_WIDTH = 640;
const TRACK_HEIGHT = 24;
const LABEL_ROW_HEIGHT = 14;
const STRUCTURE_LABEL_ROW = 16;
const MARGIN = 16;
const TRACK_TOP_PAD = 22;
const TRACK_BOTTOM_PAD = 14;
const TITLE_SPACE = 20;

function segmentTitle(kind: TranscriptStructureSegment["kind"]): string {
  if (kind === "coding") return "coding";
  if (kind === "utr5") return "5' UTR";
  return "3' UTR";
}

function ProteinTrack({
  label,
  proteinLength,
  breakpointAa,
  breakpointLabel,
  domains,
  y,
  onHover,
}: {
  label: string;
  proteinLength: number;
  breakpointAa: number | null;
  breakpointLabel?: string;
  domains: CanonDomain[];
  y: number;
  onHover: (d: CanonDomain | null) => void;
}) {
  const scale = (aa: number) => MARGIN + (aa / Math.max(proteinLength, 1)) * TRACK_WIDTH;
  const rows = labelRows(domains, proteinLength);
  const nRows = Math.max(0, ...rows.map((r) => r.row + 1));
  const bodyY = TRACK_TOP_PAD + nRows * LABEL_ROW_HEIGHT;

  return (
    <g transform={`translate(0, ${y})`}>
      <text x={MARGIN} y={-6} className="track-label">
        {label} <tspan className="track-sublabel">({proteinLength} aa)</tspan>
      </text>
      {rows.map((r) => (
        <text
          key={`${r.name}-${r.center}`}
          x={scale(r.center)}
          y={bodyY - 4 - r.row * LABEL_ROW_HEIGHT}
          className="domain-label"
          textAnchor="middle"
        >
          {r.name}
        </text>
      ))}
      <rect x={MARGIN} y={bodyY} width={TRACK_WIDTH} height={TRACK_HEIGHT} rx={4}
            fill="#e9ecef" stroke="#ced4da" />
      {domains.map((d) => {
        const x1 = scale(d.start);
        const x2 = scale(d.end);
        const opacity = d.status === "RETAINED" ? 1 : d.status === "DISRUPTED" ? 0.85 : 0.35;
        const dashed = d.status === "DISRUPTED";
        return (
          <rect
            key={`${d.name}-${d.start}-${d.end}`}
            x={x1}
            y={bodyY}
            width={Math.max(x2 - x1, 1)}
            height={TRACK_HEIGHT}
            fill={colorFor(d.name)}
            fillOpacity={opacity}
            stroke="#1a1a1a"
            strokeWidth={0.75}
            strokeDasharray={dashed ? "3,2" : undefined}
            onMouseEnter={() => onHover(d)}
            onMouseLeave={() => onHover(null)}
          >
            <title>
              {d.name} ({d.start}-{d.end}) — {d.status}
            </title>
          </rect>
        );
      })}
      {breakpointAa != null && (
        <>
          <line x1={scale(breakpointAa)} x2={scale(breakpointAa)} y1={bodyY - 6} y2={bodyY + TRACK_HEIGHT + 6}
                stroke="#e03131" strokeWidth={2} strokeDasharray="4,3" />
          {breakpointLabel && (
            <text x={scale(breakpointAa)} y={bodyY + TRACK_HEIGHT + 20} className="breakpoint-label"
                  textAnchor="middle">
              {breakpointLabel}
            </text>
          )}
        </>
      )}
    </g>
  );
}

function TranscriptTrack({
  label,
  partner,
  y,
}: {
  label: string;
  partner: ResolvedPartner;
  y: number;
}) {
  if (!partner.structure) return null;

  const layout = layoutTranscriptStructure(partner.structure);
  const scale = (pos: number) => MARGIN + (pos / Math.max(layout.width, 1)) * TRACK_WIDTH;
  const bodyY = TRACK_TOP_PAD + STRUCTURE_LABEL_ROW;
  const breakpointX = transcriptBreakpointPosition(partner, layout);
  const breakpointLabel = transcriptBreakpointLabel(partner);
  const strandLabel = partner.structure.strand === 1 ? "+ strand" : "- strand";

  return (
    <g transform={`translate(0, ${y})`}>
      <text x={MARGIN} y={-6} className="track-label">
        {label} <tspan className="track-sublabel">({partner.transcript} · {strandLabel})</tspan>
      </text>
      <text x={scale((layout.promoterStart + layout.promoterEnd) / 2)} y={bodyY - 6}
            className="domain-label" textAnchor="middle">
        promoter
      </text>
      <rect
        x={scale(layout.promoterStart)}
        y={bodyY}
        width={Math.max(scale(layout.promoterEnd) - scale(layout.promoterStart), 1)}
        height={TRACK_HEIGHT}
        fill="#fff3bf"
        stroke="#c92a2a"
        strokeWidth={0.9}
        strokeDasharray="4,2"
        rx={4}
      />
      {layout.exons.map((exon, index) => (
        <g key={`exon-${exon.rank}`}>
          {index > 0 && (
            <line
              x1={scale(layout.exons[index - 1].end)}
              x2={scale(exon.start)}
              y1={bodyY + TRACK_HEIGHT / 2}
              y2={bodyY + TRACK_HEIGHT / 2}
              stroke="#868e96"
              strokeWidth={1.5}
            />
          )}
          {index === 0 && (
            <line
              x1={scale(layout.promoterEnd)}
              x2={scale(exon.start)}
              y1={bodyY + TRACK_HEIGHT / 2}
              y2={bodyY + TRACK_HEIGHT / 2}
              stroke="#868e96"
              strokeWidth={1.5}
            />
          )}
          <text x={scale((exon.start + exon.end) / 2)} y={bodyY - 6}
                className="domain-label" textAnchor="middle">
            {exon.rank}
          </text>
          <rect
            x={scale(exon.start)}
            y={bodyY}
            width={Math.max(scale(exon.end) - scale(exon.start), 1)}
            height={TRACK_HEIGHT}
            fill="#f1f3f5"
            stroke="#495057"
            strokeWidth={0.9}
            rx={3}
          />
          {exon.segments.map((segment) => {
            const segStart = exon.start + ((segment.start - 1) / exon.length) * exon.width;
            const segEnd = exon.start + (segment.end / exon.length) * exon.width;
            return (
              <rect
                key={`${exon.rank}-${segment.kind}-${segment.start}-${segment.end}`}
                x={scale(segStart)}
                y={bodyY}
                width={Math.max(scale(segEnd) - scale(segStart), 1)}
                height={TRACK_HEIGHT}
                fill={structureSegmentColor(segment.kind)}
              >
                <title>
                  Exon {exon.rank}: {segmentTitle(segment.kind)}
                </title>
              </rect>
            );
          })}
        </g>
      ))}
      {breakpointX != null && (
        <>
          <line
            x1={scale(breakpointX)}
            x2={scale(breakpointX)}
            y1={bodyY - 6}
            y2={bodyY + TRACK_HEIGHT + 6}
            stroke="#e03131"
            strokeWidth={2}
            strokeDasharray="4,3"
          />
          <text x={scale(breakpointX)} y={bodyY + TRACK_HEIGHT + 20} className="breakpoint-label"
                textAnchor="middle">
            {breakpointLabel}
          </text>
        </>
      )}
    </g>
  );
}

function trackHeight(domains: CanonDomain[], proteinLength: number, hasBreakpointLabel: boolean) {
  const rows = labelRows(domains, proteinLength);
  const nRows = Math.max(0, ...rows.map((r) => r.row + 1));
  return TRACK_TOP_PAD + nRows * LABEL_ROW_HEIGHT + TRACK_HEIGHT + TRACK_BOTTOM_PAD + (hasBreakpointLabel ? 6 : 0);
}

function transcriptTrackHeight(hasBreakpointLabel: boolean) {
  return TRACK_TOP_PAD + STRUCTURE_LABEL_ROW + TRACK_HEIGHT + TRACK_BOTTOM_PAD + (hasBreakpointLabel ? 12 : 0);
}

export function DomainDiagram({
  domains,
  fiveGene,
  threeGene,
  fivePartner,
  threePartner,
  fiveLastAa,
  threeFirstAa,
  hybridCodon,
  fusionLength,
}: Props) {
  const [hovered, setHovered] = useState<CanonDomain | null>(null);

  const fiveDomains = canonicalizeDomains(domains, fiveGene);
  const threeDomains = canonicalizeDomains(domains, threeGene);
  const fiveLength = fivePartner.protein_length;
  const threeLength = threePartner.protein_length;

  const threeOffset = fiveLastAa + (hybridCodon ? 1 : 0) - threeFirstAa + 1;
  const fusionDomains: CanonDomain[] = [
    ...fiveDomains
      .filter((d) => d.status !== "LOST")
      .map((d) => ({ ...d, end: Math.min(d.end, fiveLastAa) })),
    ...threeDomains
      .filter((d) => d.status !== "LOST")
      .map((d) => ({
        name: d.name,
        status: "RETAINED" as Status,
        start: Math.max(d.start, threeFirstAa) + threeOffset,
        end: d.end + threeOffset,
      })),
  ].sort((a, b) => a.start - b.start);

  const s1 = fivePartner.structure ? transcriptTrackHeight(true) : 0;
  const h1 = trackHeight(fiveDomains, fiveLength, true);
  const s2 = threePartner.structure ? transcriptTrackHeight(true) : 0;
  const h2 = trackHeight(threeDomains, threeLength, true);
  const h3 = trackHeight(fusionDomains, fusionLength, false);
  const y1 = TITLE_SPACE;
  const y1Protein = y1 + s1;
  const y2 = y1Protein + h1;
  const y2Protein = y2 + s2;
  const y3 = y2Protein + h2;
  const totalHeight = y3 + h3;

  return (
    <div className="domain-diagram">
      <svg viewBox={`0 0 ${TRACK_WIDTH + MARGIN * 2} ${totalHeight}`} width="100%" role="img"
           aria-label="Transcript structure and domain retention diagram">
        {fivePartner.structure && (
          <TranscriptTrack
            label={`${fiveGene} transcript structure`}
            partner={fivePartner}
            y={y1}
          />
        )}
        <ProteinTrack
          label={fiveGene}
          proteinLength={fiveLength}
          breakpointAa={fiveLastAa}
          breakpointLabel={`breakpoint aa ${fiveLastAa}`}
          domains={fiveDomains}
          y={y1Protein}
          onHover={setHovered}
        />
        {threePartner.structure && (
          <TranscriptTrack
            label={`${threeGene} transcript structure`}
            partner={threePartner}
            y={y2}
          />
        )}
        <ProteinTrack
          label={threeGene}
          proteinLength={threeLength}
          breakpointAa={threeFirstAa}
          breakpointLabel={`breakpoint aa ${threeFirstAa}`}
          domains={threeDomains}
          y={y2Protein}
          onHover={setHovered}
        />
        <ProteinTrack
          label={`${fiveGene}::${threeGene} fusion protein`}
          proteinLength={fusionLength}
          breakpointAa={fiveLastAa}
          domains={fusionDomains}
          y={y3}
          onHover={setHovered}
        />
      </svg>
      <div className="domain-legend">
        <span><i style={{ background: "#fff3bf", border: "1px dashed #c92a2a" }} /> Promoter / 5' upstream</span>
        <span><i style={{ background: "#a5d8ff" }} /> 5' UTR</span>
        <span><i style={{ background: "#4c6ef5" }} /> Coding exon</span>
        <span><i style={{ background: "#d0ebff" }} /> 3' UTR</span>
        <span><i style={{ background: "#2f9e44" }} /> WD40 / β-propeller</span>
        <span><i style={{ background: "#e8590c" }} /> Kinase</span>
        <span><i style={{ background: "#0c8599" }} /> HELP</span>
        <span><i style={{ background: "#1971c2" }} /> MAM domain</span>
        <span><i style={{ background: "#495057" }} /> Other (stable per-name color)</span>
        <span className="legend-note">Protein domains: retained = solid · disrupted = dashed border, 85% opacity · lost = 35% opacity</span>
        <span><i className="junction-swatch" /> Breakpoint</span>
      </div>
      {hovered && (
        <div className="domain-tooltip">
          <strong>{hovered.name}</strong> ({hovered.start}-{hovered.end}) — {hovered.status}
        </div>
      )}
    </div>
  );
}
