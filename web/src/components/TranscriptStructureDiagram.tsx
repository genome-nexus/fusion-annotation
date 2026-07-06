import {
  layoutTranscriptStructure,
  structureSegmentColor,
  transcriptBreakpointLabel,
  transcriptBreakpointPosition,
} from "../lib/domainDiagram";
import type { ResolvedPartner, TranscriptStructureSegment } from "../lib/types";

interface Props {
  fiveGene: string;
  threeGene: string;
  fivePartner: ResolvedPartner;
  threePartner: ResolvedPartner;
}

const TRACK_WIDTH = 640;
const TRACK_HEIGHT = 24;
const STRUCTURE_LABEL_ROW = 16;
const MARGIN = 16;
const TRACK_TOP_PAD = 22;
const TRACK_BOTTOM_PAD = 20;
const TITLE_SPACE = 20;

function segmentTitle(kind: TranscriptStructureSegment["kind"]): string {
  if (kind === "coding") return "coding";
  if (kind === "utr5") return "5' UTR";
  return "3' UTR";
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

function transcriptTrackHeight() {
  return TRACK_TOP_PAD + STRUCTURE_LABEL_ROW + TRACK_HEIGHT + TRACK_BOTTOM_PAD + 12;
}

export function TranscriptStructureDiagram({
  fiveGene,
  threeGene,
  fivePartner,
  threePartner,
}: Props) {
  const hasStructure = Boolean(fivePartner.structure || threePartner.structure);
  if (!hasStructure) return null;

  const h1 = fivePartner.structure ? transcriptTrackHeight() : 0;
  const h2 = threePartner.structure ? transcriptTrackHeight() : 0;
  const y1 = TITLE_SPACE;
  const y2 = y1 + h1;
  const totalHeight = y2 + h2;

  return (
    <div className="domain-diagram">
      <svg viewBox={`0 0 ${TRACK_WIDTH + MARGIN * 2} ${totalHeight}`} width="100%" role="img"
           aria-label="Transcript structure diagram">
        {fivePartner.structure && (
          <TranscriptTrack
            label={`${fiveGene} transcript structure`}
            partner={fivePartner}
            y={y1}
          />
        )}
        {threePartner.structure && (
          <TranscriptTrack
            label={`${threeGene} transcript structure`}
            partner={threePartner}
            y={y2}
          />
        )}
      </svg>
      <div className="domain-legend">
        <span><i style={{ background: "#fff3bf", border: "1px dashed #c92a2a" }} /> Promoter / 5' upstream</span>
        <span><i style={{ background: "#a5d8ff" }} /> 5' UTR</span>
        <span><i style={{ background: "#4c6ef5" }} /> Coding exon</span>
        <span><i style={{ background: "#d0ebff" }} /> 3' UTR</span>
        <span><i className="junction-swatch" /> Breakpoint</span>
      </div>
    </div>
  );
}
