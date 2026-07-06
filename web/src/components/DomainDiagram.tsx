import { useState } from "react";
import { canonicalizeDomains, colorFor, labelRows, type CanonDomain, type Status } from "../lib/domainDiagram";
import type { DomainCall } from "../lib/types";

interface Props {
  domains: DomainCall[];
  fiveGene: string;
  threeGene: string;
  fiveLastAa: number;
  threeFirstAa: number;
  fiveLength: number;
  threeLength: number;
  hybridCodon: boolean;
  fusionLength: number;
}

const TRACK_WIDTH = 640;
const TRACK_HEIGHT = 24;
const LABEL_ROW_HEIGHT = 14;
const MARGIN = 16;
const TRACK_TOP_PAD = 22;
const TRACK_BOTTOM_PAD = 14;
const TITLE_SPACE = 20;

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

function trackHeight(domains: CanonDomain[], proteinLength: number, hasBreakpointLabel: boolean) {
  const rows = labelRows(domains, proteinLength);
  const nRows = Math.max(0, ...rows.map((r) => r.row + 1));
  return TRACK_TOP_PAD + nRows * LABEL_ROW_HEIGHT + TRACK_HEIGHT + TRACK_BOTTOM_PAD + (hasBreakpointLabel ? 6 : 0);
}

export function DomainDiagram({
  domains,
  fiveGene,
  threeGene,
  fiveLastAa,
  threeFirstAa,
  fiveLength,
  threeLength,
  hybridCodon,
  fusionLength,
}: Props) {
  const [hovered, setHovered] = useState<CanonDomain | null>(null);

  const fiveDomains = canonicalizeDomains(domains, fiveGene);
  const threeDomains = canonicalizeDomains(domains, threeGene);

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

  const h1 = trackHeight(fiveDomains, fiveLength, true);
  const h2 = trackHeight(threeDomains, threeLength, true);
  const h3 = trackHeight(fusionDomains, fusionLength, false);
  const y1 = TITLE_SPACE;
  const y2 = y1 + h1;
  const y3 = y2 + h2;
  const totalHeight = y3 + h3;

  return (
    <div className="domain-diagram">
      <svg viewBox={`0 0 ${TRACK_WIDTH + MARGIN * 2} ${totalHeight}`} width="100%" role="img"
           aria-label="Domain retention diagram">
        <ProteinTrack
          label={fiveGene}
          proteinLength={fiveLength}
          breakpointAa={fiveLastAa}
          breakpointLabel={`breakpoint aa ${fiveLastAa}`}
          domains={fiveDomains}
          y={y1}
          onHover={setHovered}
        />
        <ProteinTrack
          label={threeGene}
          proteinLength={threeLength}
          breakpointAa={threeFirstAa}
          breakpointLabel={`breakpoint aa ${threeFirstAa}`}
          domains={threeDomains}
          y={y2}
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
        <span><i style={{ background: "#2f9e44" }} /> WD40 / β-propeller</span>
        <span><i style={{ background: "#e8590c" }} /> Kinase</span>
        <span><i style={{ background: "#0c8599" }} /> HELP</span>
        <span><i style={{ background: "#1971c2" }} /> MAM domain</span>
        <span><i style={{ background: "#495057" }} /> Other (stable per-name color)</span>
        <span className="legend-note">Retained = solid · Disrupted = dashed border, 85% opacity · Lost = 35% opacity</span>
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
