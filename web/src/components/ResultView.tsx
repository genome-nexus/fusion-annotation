import { useState } from "react";
import type { AnnotationResult } from "../lib/types";
import { computeDerivedInputs, type DerivedInputs, type PartnerDerived } from "../lib/derivedInputs";
import { DomainDiagram } from "./DomainDiagram";
import { DomainTable } from "./DomainTable";
import { TranscriptStructureDiagram } from "./TranscriptStructureDiagram";

interface Props {
  result: AnnotationResult;
  permalink: string;
}

type DiagramView = "domain" | "structure";

function CopyButton({ text }: { text: string }) {
  const [copied, setCopied] = useState(false);
  return (
    <button
      type="button"
      className="copy-btn"
      title="Copy"
      onClick={() => {
        navigator.clipboard?.writeText(text);
        setCopied(true);
        setTimeout(() => setCopied(false), 1200);
      }}
    >
      {copied ? "✓" : "⧉"}
    </button>
  );
}

function PartnerBreakpoints({ label, p, aa }: { label: string; p: PartnerDerived; aa: number }) {
  const rows: Array<{ name: string; value: string | null }> = [
    { name: "Genomic", value: p.genomic },
    { name: "Exon", value: p.exon != null ? `exon ${p.exon}` : null },
    { name: "CDS", value: `c.${p.cds}` },
    { name: "Protein", value: `aa ${aa}` },
  ];
  return (
    <div className="equiv-partner">
      <div className="equiv-partner-label">{label}</div>
      {rows.map(({ name, value }) =>
        value ? (
          <div key={name} className="equiv-row">
            <span className="equiv-kind">{name}</span>
            <code className="equiv-val">{value}</code>
            <CopyButton text={value} />
          </div>
        ) : null,
      )}
    </div>
  );
}

function EquivalentInputs({ derived, result }: { derived: DerivedInputs; result: AnnotationResult }) {
  const { interface: iface } = result;
  return (
    <details className="equiv-inputs">
      <summary>Equivalent breakpoint representations</summary>
      <div className="equiv-grid">
        <PartnerBreakpoints label={`5′ ${derived.five.gene}`} p={derived.five} aa={iface.five_last_aa} />
        <PartnerBreakpoints label={`3′ ${derived.three.gene}`} p={derived.three} aa={iface.three_first_aa} />
      </div>
    </details>
  );
}

export function ResultView({ result, permalink }: Props) {
  const { interface: iface, knowledge, resolved, warnings } = result;
  const [diagramView, setDiagramView] = useState<DiagramView>("domain");
  const hasTranscriptStructure = Boolean(resolved.five.structure || resolved.three.structure);
  const derived = computeDerivedInputs(result);

  return (
    <div className="result-view">
      <div className="permalink-row">
        <input readOnly value={permalink} onFocus={(e) => e.currentTarget.select()} />
        <button
          type="button"
          onClick={() => navigator.clipboard?.writeText(permalink)}
        >
          Copy permalink
        </button>
      </div>

      {warnings.length > 0 && (
        <ul className="warnings">
          {warnings.map((w) => (
            <li key={w}>⚠ {w}</li>
          ))}
        </ul>
      )}

      <h2 className="hgvsp">{iface.hgvsp_like}</h2>

      <dl className="summary-grid">
        <dt>Categorical fusion</dt>
        <dd>{iface.categorical_key}</dd>
        <dt>Frame</dt>
        <dd className={iface.in_frame ? "ok" : "bad"}>
          {iface.frame_status} (protein {iface.fusion_length} aa, internal stops {iface.internal_stops})
        </dd>
        <dt>Junction</dt>
        <dd>
          {iface.five_gene} {iface.five_last_aa_res}
          {iface.five_last_aa} :: {iface.three_gene} {iface.three_first_aa_res}
          {iface.three_first_aa}
          {iface.hybrid_codon && iface.junction_residue
            ? ` (hybrid codon → ${iface.junction_residue})`
            : ""}
        </dd>
        <dt>Genome build</dt>
        <dd>{resolved.genome_build}</dd>
        <dt>Transcripts used</dt>
        <dd>
          {resolved.five.transcript} ({resolved.five.transcript_source}) / {resolved.three.transcript} (
          {resolved.three.transcript_source})
        </dd>
      </dl>

      <h3>Visualization</h3>
      <div className="diagram-toggle" role="tablist" aria-label="Choose visualization">
        <button
          type="button"
          className={diagramView === "domain" ? "diagram-toggle-button active" : "diagram-toggle-button"}
          aria-pressed={diagramView === "domain"}
          onClick={() => setDiagramView("domain")}
        >
          Domain retention
        </button>
        {hasTranscriptStructure && (
          <button
            type="button"
            className={diagramView === "structure" ? "diagram-toggle-button active" : "diagram-toggle-button"}
            aria-pressed={diagramView === "structure"}
            onClick={() => setDiagramView("structure")}
          >
            Transcript structure
          </button>
        )}
      </div>
      {diagramView === "domain" ? (
        <DomainDiagram
          domains={iface.domains}
          fiveGene={iface.five_gene}
          threeGene={iface.three_gene}
          fiveLastAa={iface.five_last_aa}
          threeFirstAa={iface.three_first_aa}
          fiveLength={resolved.five.protein_length}
          threeLength={resolved.three.protein_length}
          hybridCodon={iface.hybrid_codon}
          fusionLength={iface.fusion_length}
        />
      ) : (
        <TranscriptStructureDiagram
          fiveGene={iface.five_gene}
          threeGene={iface.three_gene}
          fivePartner={resolved.five}
          threePartner={resolved.three}
        />
      )}

      <EquivalentInputs derived={derived} result={result} />

      <DomainTable domains={iface.domains} />

      <h3>Clinical knowledge</h3>
      <dl className="summary-grid">
        <dt>Oncogenic</dt>
        <dd>{knowledge.oncogenic ?? "—"}</dd>
        <dt>Therapies</dt>
        <dd>{knowledge.therapies.length ? knowledge.therapies.join(", ") : "—"}</dd>
        <dt>Diseases</dt>
        <dd>{knowledge.diseases.length ? knowledge.diseases.join(", ") : "—"}</dd>
        <dt>Sources</dt>
        <dd>{knowledge.sources.length ? knowledge.sources.join(", ") : "—"}</dd>
      </dl>

      <p className="disclaimer">
        This is a research/informatics tool, not a diagnostic device. Results should be reviewed by a
        qualified professional before they inform patient care.
      </p>
    </div>
  );
}
