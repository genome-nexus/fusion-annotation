import { useState } from "react";
import type {
  AnnotationResult,
  ApiError,
  FusionKnowledge,
  GeneCurationFusionResult,
  GeneCurationGeneResult,
  GeneFusionCurationContext,
} from "../lib/types";
import { computeDerivedInputs, type DerivedInputs, type PartnerDerived } from "../lib/derivedInputs";
import { civicEvidenceLink, civicMolecularProfileLink } from "../lib/externalLinks";
import { DomainDiagram } from "./DomainDiagram";
import { DomainTable } from "./DomainTable";
import { TranscriptStructureDiagram } from "./TranscriptStructureDiagram";

interface Props {
  result: AnnotationResult;
  permalink: string;
  fusionCurationResults?: GeneCurationFusionResult[] | null;
  geneCurationResults?: GeneCurationGeneResult[] | null;
  geneCurationLoading?: boolean;
  geneCurationEnabled?: boolean;
  geneCurationError?: ApiError | null;
  onCurateGenes?: () => void;
  onForceGeneCuration?: () => void;
  onExportGeneCurationCsv?: () => void;
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

function PartnerBreakpoints({ label, p, aa }: { label: string; p: PartnerDerived; aa: number | null }) {
  const rows: Array<{ name: string; value: string | null }> = [
    { name: "Genomic", value: p.genomic },
    { name: "Exon", value: p.exon != null ? `exon ${p.exon}` : null },
    { name: "CDS", value: p.cds != null ? `c.${p.cds}` : null },
    { name: "Protein", value: aa != null ? `aa ${aa}` : null },
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

function ClinicalEvidence({ knowledge }: { knowledge: FusionKnowledge }) {
  if (!knowledge.evidence || knowledge.evidence.length === 0) return null;
  return (
    <details className="clinical-evidence">
      <summary>
        CIViC evidence ({knowledge.evidence.length})
      </summary>
      <ul className="evidence-list">
        {knowledge.evidence.map((ev, idx) => {
          const id = ev.id as number | undefined;
          const type = ev.type as string | undefined;
          const level = ev.level as string | undefined;
          const disease = ev.disease as string | undefined;
          const therapies = (ev.therapies as string[] | undefined) || [];
          return (
            <li key={idx} className="evidence-item">
              {id ? (
                <a href={civicEvidenceLink(id).url} target="_blank" rel="noopener noreferrer">
                  <strong>Evidence #{id}</strong>
                </a>
              ) : (
                <strong>Evidence</strong>
              )}
              {type && <span className="evidence-type">{type}</span>}
              {level && <span className="evidence-level">Level {level}</span>}
              {disease && <span className="evidence-disease">{disease}</span>}
              {therapies.length > 0 && <span className="evidence-therapies">{therapies.join(", ")}</span>}
            </li>
          );
        })}
      </ul>
    </details>
  );
}

function curationPriority(item: GeneCurationGeneResult | GeneCurationFusionResult) {
  if (item.insufficient_evidence) {
    return {
      label: "Low priority",
      tone: "muted",
      title: "PubMed evidence was too sparse for a confident curation result.",
    };
  }
  if (item.cancer_associated === true) {
    return {
      label: "Review priority",
      tone: "attention",
      title: "Cancer-associated literature evidence was found; review before follow-up.",
    };
  }
  if (item.cancer_associated === false) {
    return {
      label: "Low priority",
      tone: "muted",
      title: "Current PubMed evidence does not support a cancer association.",
    };
  }
  return {
    label: "Needs review",
    tone: "neutral",
    title: "The curation result is uncertain and should be reviewed.",
  };
}

function curationEvidenceSignal(item: GeneCurationGeneResult | GeneCurationFusionResult) {
  if (item.insufficient_evidence) {
    return {
      label: "Sparse evidence",
      tone: "muted",
      title: "The curation model marked the retrieved literature as insufficient.",
    };
  }
  if (item.cancer_associated === true) {
    return {
      label: "Functional cancer evidence",
      tone: "evidence",
      title: "Cancer-associated evidence is shown without changing backend schema or tier logic.",
    };
  }
  if (item.cancer_associated === false) {
    return {
      label: "No cancer evidence",
      tone: "muted",
      title: "The curation result did not find supporting cancer evidence.",
    };
  }
  return {
    label: "Uncertain evidence",
    tone: "neutral",
    title: "The curation result did not make a binary cancer association call.",
  };
}

function curationBadges(item: GeneCurationGeneResult | GeneCurationFusionResult) {
  return [curationPriority(item), curationEvidenceSignal(item)];
}

function formatContextSide(context: GeneFusionCurationContext) {
  return context.side === "five_prime" ? "5' partner" : "3' partner";
}

function formatDomainList(domains?: string[]) {
  return domains && domains.length ? domains.join(", ") : "none";
}

function formatFusionSpecificity(value?: GeneFusionCurationContext["fusion_specificity"]) {
  if (value === "protein_domain_level") return "Protein/domain-level";
  if (value === "exon_level") return "Exon-level";
  return "Gene-pair only";
}

function pubmedUrl(pmid: string) {
  return `https://pubmed.ncbi.nlm.nih.gov/${encodeURIComponent(pmid)}/`;
}

function confidenceLabel(item?: GeneCurationGeneResult | GeneCurationFusionResult) {
  if (!item || item.insufficient_evidence) return "Low";
  if (item.cancer_associated === true) return "Higher";
  if (item.cancer_associated === false) return "Lower";
  return "Uncertain";
}

function renderFusionContext(context: GeneFusionCurationContext) {
  const breakpoint = context.side === "five_prime"
    ? {
        transcript: context.five_transcript,
        exon: context.five_exon,
        genomic: context.five_genomic,
        protein: context.five_protein_breakpoint,
      }
    : {
        transcript: context.three_transcript,
        exon: context.three_exon,
        genomic: context.three_genomic,
        protein: context.three_protein_breakpoint,
      };
  return (
    <div className="fusion-curation-context" key={`${context.gene}-${context.fusion}-${context.side}`}>
      <div className="fusion-curation-context-title">
        <strong>{context.fusion}</strong>
        <span>{formatContextSide(context)}</span>
      </div>
      <dl>
        <dt>Scope</dt>
        <dd>
          {formatFusionSpecificity(context.fusion_specificity)}
          {context.breakpoint_context_available ? "" : " · breakpoint context unavailable"}
        </dd>
        <dt>Partner</dt>
        <dd>{context.partner_gene}</dd>
        <dt>Breakpoint</dt>
        <dd>
          tx {breakpoint.transcript || "unknown"} · exon {breakpoint.exon || "unknown"} · genomic{" "}
          {breakpoint.genomic || "unknown"} · protein {breakpoint.protein || "unknown"}
        </dd>
        <dt>Domains</dt>
        <dd>
          retained: {formatDomainList(context.retained_domains)} · lost/disrupted:{" "}
          {formatDomainList([...(context.lost_domains || []), ...(context.disrupted_domains || [])])}
        </dd>
        <dt>Kinase</dt>
        <dd>
          {context.kinase_gene || "unknown"} · {context.kinase_domain_status || "unknown"}
        </dd>
      </dl>
      {context.limitations && context.limitations.length > 0 && (
        <ul className="fusion-curation-limitations">
          {context.limitations.map((limitation) => (
            <li key={limitation}>{limitation}</li>
          ))}
        </ul>
      )}
      {context.annotation_error && (
        <p className="fusion-curation-context-error">{context.annotation_error}</p>
      )}
    </div>
  );
}

function GeneInformationSection({
  result,
  fusionCurationResults,
  geneCurationResults,
  geneCurationLoading = false,
  geneCurationEnabled = false,
  geneCurationError = null,
  onCurateGenes,
  onForceGeneCuration,
  onExportGeneCurationCsv,
}: {
  result: AnnotationResult;
  fusionCurationResults?: GeneCurationFusionResult[] | null;
  geneCurationResults?: GeneCurationGeneResult[] | null;
  geneCurationLoading?: boolean;
  geneCurationEnabled?: boolean;
  geneCurationError?: ApiError | null;
  onCurateGenes?: () => void;
  onForceGeneCuration?: () => void;
  onExportGeneCurationCsv?: () => void;
}) {
  const genes = [
    result.resolved.five.gene,
    result.resolved.three.gene,
  ].filter((gene, index, values) => values.indexOf(gene) === index);
  const byGene = new Map(
    (geneCurationResults || []).map((item) => [item.gene.toUpperCase(), item]),
  );
  const fusionItem = (fusionCurationResults || []).find(
    (item) => item.fusion === result.interface.categorical_key,
  );
  const fusionSufficient = Boolean(
    fusionItem
      && !fusionItem.error
      && !fusionItem.insufficient_evidence
      && fusionItem.fusion_literature_identified !== false,
  );

  return (
    <section className="gene-info-section">
      <div className="gene-info-header">
        <div>
          <h3>Gene information</h3>
          <p>
            Literature-backed gene context from server-side curation, shown alongside the current fusion
            annotation.
          </p>
        </div>
        <div className="gene-info-actions">
          {geneCurationResults?.length || fusionCurationResults?.length ? (
            <button type="button" className="secondary-button" onClick={onExportGeneCurationCsv}>
              Export curation CSV
            </button>
          ) : null}
          {onCurateGenes && (
            <button
              type="button"
              className="secondary-button"
              disabled={!geneCurationEnabled || geneCurationLoading}
              onClick={onCurateGenes}
            >
              {geneCurationLoading ? "Loading literature..." : "Get literature info"}
            </button>
          )}
        </div>
      </div>

      {!geneCurationEnabled && (
        <div className="notice-box">
          Server-side gene curation is not configured for this deployment.
        </div>
      )}
      {geneCurationError && (
        <div className="error-box">
          <strong>Curation error {geneCurationError.status}:</strong> {geneCurationError.detail}
        </div>
      )}

      <div className="gene-info-list">
        <details className="gene-info-item" open={Boolean(fusionItem)}>
          <summary>
            <span className="gene-info-symbol">{result.interface.categorical_key}</span>
            <span className="gene-info-summary">
              {fusionItem
                ? `${fusionItem.fusion_literature_identified === false ? "Fusion literature not found" : "Fusion literature found"} · confidence ${confidenceLabel(fusionItem)}`
                : "No fusion-specific literature curation loaded"}
            </span>
          </summary>
          {fusionItem ? (
            fusionItem.error ? (
              <div className="error-box">{fusionItem.error}</div>
            ) : (
              <div className="gene-info-body">
                <div className="curation-badges" aria-label={`${fusionItem.fusion} review signals`}>
                  {curationBadges(fusionItem).map((badge) => (
                    <span className={`status-chip ${badge.tone}`} key={badge.label} title={badge.title}>
                      {badge.label}
                    </span>
                  ))}
                </div>
                <dl className="gene-curation-fields">
                  <dt>Fusion in literature</dt>
                  <dd>
                    {fusionItem.fusion_literature_identified == null
                      ? "Unknown"
                      : fusionItem.fusion_literature_identified ? "Yes" : "No"}
                  </dd>
                  <dt>Known driver signal</dt>
                  <dd>
                    {fusionItem.cancer_associated == null ? "Unknown" : fusionItem.cancer_associated ? "Yes" : "No"}
                  </dd>
                  <dt>Confidence</dt>
                  <dd>{confidenceLabel(fusionItem)}</dd>
                  <dt>Rationale</dt>
                  <dd>{fusionItem.rationale || "No rationale returned."}</dd>
                </dl>
                {fusionItem.fusion_contexts && fusionItem.fusion_contexts.length > 0 && (
                  <div className="fusion-curation-contexts">
                    {fusionItem.fusion_contexts.map(renderFusionContext)}
                  </div>
                )}
                <div className="pmid-row">
                  <strong>Supporting PMIDs</strong>
                  <span>
                    {fusionItem.supporting_pmids?.length
                      ? fusionItem.supporting_pmids.map((pmid, index) => (
                          <span key={pmid}>
                            {index > 0 && ", "}
                            <a href={pubmedUrl(pmid)} target="_blank" rel="noopener noreferrer">
                              {pmid}
                            </a>
                          </span>
                        ))
                      : "None selected"}
                  </span>
                </div>
                <div className="pmid-row">
                  <strong>Retrieved PMIDs</strong>
                  <span>{fusionItem.retrieved_pmids?.join(", ") || "None retrieved"}</span>
                </div>
              </div>
            )
          ) : (
            <p className="gene-info-empty">
              Run literature info to check whether this exact fusion is described in PubMed.
            </p>
          )}
        </details>

        {fusionSufficient && !geneCurationResults?.length && onForceGeneCuration && (
          <div className="notice-box gene-info-skip">
            Fusion-specific literature was sufficient, so per-gene retrieval was skipped to reduce PubMed and LLM
            usage.
            <button
              type="button"
              className="secondary-button"
              disabled={geneCurationLoading}
              onClick={onForceGeneCuration}
            >
              Get per-gene info
            </button>
          </div>
        )}

        {genes.map((gene) => {
          const item = byGene.get(gene.toUpperCase());
          const contexts = item?.fusion_contexts?.filter(
            (context) => context.fusion === result.interface.categorical_key,
          );
          return (
            <details className="gene-info-item" key={gene}>
              <summary>
                <span className="gene-info-symbol">{gene}</span>
                <span className="gene-info-summary">
                  {item
                    ? `${item.cancer_associated == null ? "Cancer association unknown" : item.cancer_associated ? "Cancer associated" : "No cancer association found"} · confidence ${confidenceLabel(item)}`
                    : "No literature curation loaded"}
                </span>
              </summary>
              {item ? (
                item.error ? (
                  <div className="error-box">{item.error}</div>
                ) : (
                  <div className="gene-info-body">
                    <div className="curation-badges" aria-label={`${gene} review signals`}>
                      {curationBadges(item).map((badge) => (
                        <span className={`status-chip ${badge.tone}`} key={badge.label} title={badge.title}>
                          {badge.label}
                        </span>
                      ))}
                    </div>
                    <dl className="gene-curation-fields">
                      <dt>Known driver signal</dt>
                      <dd>{item.cancer_associated == null ? "Unknown" : item.cancer_associated ? "Yes" : "No"}</dd>
                      <dt>Confidence</dt>
                      <dd>{confidenceLabel(item)}</dd>
                      <dt>Rationale</dt>
                      <dd>{item.rationale || "No rationale returned."}</dd>
                      <dt>Fusion knowledge</dt>
                      <dd>{result.knowledge.oncogenic ?? "No fusion-level knowledge-base signal returned."}</dd>
                    </dl>
                    {contexts && contexts.length > 0 && (
                      <div className="fusion-curation-contexts">
                        {contexts.map(renderFusionContext)}
                      </div>
                    )}
                    <div className="pmid-row">
                      <strong>Supporting PMIDs</strong>
                      <span>
                        {item.supporting_pmids?.length
                          ? item.supporting_pmids.map((pmid, index) => (
                              <span key={pmid}>
                                {index > 0 && ", "}
                                <a href={pubmedUrl(pmid)} target="_blank" rel="noopener noreferrer">
                                  {pmid}
                                </a>
                              </span>
                            ))
                          : "None selected"}
                      </span>
                    </div>
                    <div className="pmid-row">
                      <strong>Retrieved PMIDs</strong>
                      <span>{item.retrieved_pmids?.join(", ") || "None retrieved"}</span>
                    </div>
                  </div>
                )
              ) : (
                <p className="gene-info-empty">
                  {fusionSufficient
                    ? `Per-gene retrieval for ${gene} was skipped because fusion-specific literature was sufficient.`
                    : `Run literature info to retrieve gene-level AGCG context for ${gene}.`}
                </p>
              )}
            </details>
          );
        })}
      </div>
    </section>
  );
}

export function ResultView({
  result,
  permalink,
  fusionCurationResults,
  geneCurationResults,
  geneCurationLoading,
  geneCurationEnabled,
  geneCurationError,
  onCurateGenes,
  onForceGeneCuration,
  onExportGeneCurationCsv,
}: Props) {
  const { interface: iface, knowledge, resolved, warnings } = result;
  const [diagramView, setDiagramView] = useState<DiagramView>("domain");
  const hasTranscriptStructure = Boolean(resolved.five.structure || resolved.three.structure);
  const derived = computeDerivedInputs(result);
  const hasProteinBreakpoint = iface.five_last_aa != null && iface.three_first_aa != null;
  const frameClass = iface.in_frame == null ? "neutral" : iface.in_frame ? "ok" : "bad";

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
        <dd className={frameClass}>
          {iface.frame_status === "unknown"
            ? "Unknown; breakpoint required for protein reconstruction"
            : `${iface.frame_status} (protein ${iface.fusion_length} aa, internal stops ${iface.internal_stops})`}
        </dd>
        <dt>Junction</dt>
        <dd>
          {hasProteinBreakpoint
            ? (
                <>
                  {iface.five_gene} {iface.five_last_aa_res}
                  {iface.five_last_aa} :: {iface.three_gene} {iface.three_first_aa_res}
                  {iface.three_first_aa}
                  {iface.hybrid_codon && iface.junction_residue
                    ? ` (hybrid codon → ${iface.junction_residue})`
                    : ""}
                </>
              )
            : "Unknown; exon or genomic breakpoint not supplied"}
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

      {hasProteinBreakpoint ? (
        <EquivalentInputs derived={derived} result={result} />
      ) : (
        <div className="notice-box">
          Equivalent breakpoint representations are unavailable until exon or genomic breakpoints are supplied.
        </div>
      )}

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
        <dd>
          {knowledge.sources.length ? (
            <span>
              {knowledge.sources.map((source, idx) => {
                // Parse source strings like "CIViC MP 5" → Molecular Profile link
                const match = source.match(/^CIViC\s+MP\s+(\d+)$/i);
                if (match) {
                  const mpId = parseInt(match[1], 10);
                  return (
                    <span key={idx}>
                      {idx > 0 && ", "}
                      <a href={civicMolecularProfileLink(mpId).url} target="_blank" rel="noopener noreferrer">
                        {source}
                      </a>
                    </span>
                  );
                }
                return <span key={idx}>{idx > 0 ? ", " : ""}{source}</span>;
              })}
            </span>
          ) : (
            "—"
          )}
        </dd>
      </dl>

      <ClinicalEvidence knowledge={knowledge} />

      <GeneInformationSection
        result={result}
        fusionCurationResults={fusionCurationResults}
        geneCurationResults={geneCurationResults}
        geneCurationLoading={geneCurationLoading}
        geneCurationEnabled={geneCurationEnabled}
        geneCurationError={geneCurationError}
        onCurateGenes={onCurateGenes}
        onForceGeneCuration={onForceGeneCuration}
        onExportGeneCurationCsv={onExportGeneCurationCsv}
      />

      <p className="disclaimer">
        This is a research/informatics tool, not a diagnostic device. Results should be reviewed by a
        qualified professional before they inform patient care.
      </p>
    </div>
  );
}
