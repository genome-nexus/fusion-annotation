import { useState } from "react";
import type { DomainCall } from "../lib/types";

interface Props {
  domains: DomainCall[];
  fiveGene: string;
  threeGene: string;
  fiveLastAa: number;
  threeFirstAa: number;
}

const STATUS_COLOR: Record<string, string> = {
  RETAINED: "#2f9e44",
  LOST: "#adb5bd",
  DISRUPTED: "#f08c00",
};

const TRACK_WIDTH = 640;
const TRACK_HEIGHT = 26;
const TRACK_GAP = 56;
const MARGIN = 16;

/** One partner's protein track: a grey full-length backbone, colored domain
 * rectangles on top, and a junction marker at the breakpoint residue. */
function ProteinTrack({
  label,
  proteinLength,
  breakpointAa,
  domains,
  y,
  onHover,
}: {
  label: string;
  proteinLength: number;
  breakpointAa: number;
  domains: DomainCall[];
  y: number;
  onHover: (d: DomainCall | null) => void;
}) {
  const scale = (aa: number) => MARGIN + (aa / Math.max(proteinLength, 1)) * TRACK_WIDTH;

  return (
    <g transform={`translate(0, ${y})`}>
      <text x={MARGIN} y={-8} className="track-label">
        {label} <tspan className="track-sublabel">({proteinLength} aa)</tspan>
      </text>
      {/* backbone */}
      <rect
        x={MARGIN}
        y={0}
        width={TRACK_WIDTH}
        height={TRACK_HEIGHT}
        rx={4}
        fill="#e9ecef"
        stroke="#ced4da"
      />
      {domains.map((d) => {
        const x1 = scale(d.start);
        const x2 = scale(d.end);
        return (
          <rect
            key={`${d.accession}-${d.start}-${d.end}`}
            x={x1}
            y={0}
            width={Math.max(x2 - x1, 1)}
            height={TRACK_HEIGHT}
            fill={STATUS_COLOR[d.status] ?? "#495057"}
            opacity={d.status === "LOST" ? 0.5 : 0.9}
            onMouseEnter={() => onHover(d)}
            onMouseLeave={() => onHover(null)}
          >
            <title>
              {d.name} ({d.start}-{d.end}) — {d.status}
            </title>
          </rect>
        );
      })}
      {/* junction marker */}
      <line
        x1={scale(breakpointAa)}
        x2={scale(breakpointAa)}
        y1={-4}
        y2={TRACK_HEIGHT + 4}
        stroke="#e03131"
        strokeWidth={2}
      />
    </g>
  );
}

/** Interactive SVG domain-retention diagram — the in-browser equivalent of
 * docs/fusion_domain_map.png: one track per partner, domains colored by
 * retained/lost/disrupted status, with the breakpoint marked on each. */
export function DomainDiagram({ domains, fiveGene, threeGene, fiveLastAa, threeFirstAa }: Props) {
  const [hovered, setHovered] = useState<DomainCall | null>(null);
  const fiveDomains = domains.filter((d) => d.gene === fiveGene);
  const threeDomains = domains.filter((d) => d.gene === threeGene);
  const fiveLength = fiveDomains[0]?.partner_protein_length || fiveLastAa;
  const threeLength = threeDomains[0]?.partner_protein_length || threeFirstAa;

  const height = TRACK_GAP * 2 + 20;

  return (
    <div className="domain-diagram">
      <svg viewBox={`0 0 ${TRACK_WIDTH + MARGIN * 2} ${height}`} width="100%" role="img"
           aria-label="Domain retention diagram">
        <ProteinTrack
          label={fiveGene}
          proteinLength={fiveLength}
          breakpointAa={fiveLastAa}
          domains={fiveDomains}
          y={24}
          onHover={setHovered}
        />
        <ProteinTrack
          label={threeGene}
          proteinLength={threeLength}
          breakpointAa={threeFirstAa}
          domains={threeDomains}
          y={24 + TRACK_GAP}
          onHover={setHovered}
        />
      </svg>
      <div className="domain-legend">
        <span><i style={{ background: STATUS_COLOR.RETAINED }} /> Retained</span>
        <span><i style={{ background: STATUS_COLOR.DISRUPTED }} /> Disrupted</span>
        <span><i style={{ background: STATUS_COLOR.LOST, opacity: 0.5 }} /> Lost</span>
        <span><i className="junction-swatch" /> Breakpoint</span>
      </div>
      {hovered && (
        <div className="domain-tooltip">
          <strong>{hovered.name}</strong> ({hovered.gene} {hovered.start}-{hovered.end}) — {hovered.status}
        </div>
      )}
    </div>
  );
}
