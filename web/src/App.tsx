import { useCallback, useEffect, useRef, useState } from "react";
import "./App.css";
import { BatchFusionForm } from "./components/BatchFusionForm";
import { ExampleFusions } from "./components/ExampleFusions";
import { FusionForm } from "./components/FusionForm";
import { ResultView } from "./components/ResultView";
import { VersionFootnote } from "./components/VersionFootnote";
import {
  annotateFusion,
  annotateFusionBatch,
  curateFusionGenes,
  getGeneCurationStatus,
  toSearchParams,
} from "./lib/api";
import { DEFAULT_PARAMS } from "./lib/defaultParams";
import { computeDerivedInputs, type DerivedInputs } from "./lib/derivedInputs";
import type {
  AnnotateParams,
  AnnotationResult,
  ApiError,
  BatchAnnotationItemResult,
  GeneFusionCurationContext,
  GeneCurationGeneResult,
  GeneCurationStatus,
} from "./lib/types";

/** Read the current URL's query string into AnnotateParams — the other half
 * of the stateless-permalink contract (writing happens in runAnnotation
 * below via history.pushState). No fields are persisted server-side: this
 * URL alone is enough to reproduce the annotation. */
function paramsFromLocation(): AnnotateParams {
  const search = new URLSearchParams(window.location.search);
  const params = { ...DEFAULT_PARAMS };
  for (const key of Object.keys(params) as (keyof AnnotateParams)[]) {
    const value = search.get(key);
    if (value !== null) (params as Record<string, unknown>)[key] = value;
  }
  return params;
}

function curationPriority(item: GeneCurationGeneResult) {
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

function curationEvidenceSignal(item: GeneCurationGeneResult) {
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

function curationBadges(item: GeneCurationGeneResult) {
  return [curationPriority(item), curationEvidenceSignal(item)];
}

function formatContextSide(context: GeneFusionCurationContext) {
  return context.side === "five_prime" ? "5' partner" : "3' partner";
}

function formatDomainList(domains?: string[]) {
  return domains && domains.length ? domains.join(", ") : "none";
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
      {context.annotation_error && (
        <p className="fusion-curation-context-error">{context.annotation_error}</p>
      )}
    </div>
  );
}

const CURATION_CSV_HEADERS = [
  "gene",
  "fusion",
  "gene_side",
  "partner_gene",
  "five_transcript",
  "three_transcript",
  "five_exon",
  "three_exon",
  "five_genomic",
  "three_genomic",
  "five_protein_breakpoint",
  "three_protein_breakpoint",
  "retained_domains",
  "lost_or_disrupted_domains",
  "kinase_gene",
  "kinase_gene_side",
  "kinase_domain_status",
  "cancer_associated",
  "rationale",
  "supporting_pmids",
  "retrieved_pmids",
  "insufficient_evidence",
  "error",
];

function csvCell(value: unknown) {
  if (value == null) return "";
  const text = Array.isArray(value) ? value.join("; ") : String(value);
  return `"${text.replaceAll('"', '""')}"`;
}

function curationCsvRows(items: GeneCurationGeneResult[]) {
  const rows = [CURATION_CSV_HEADERS.map(csvCell).join(",")];
  for (const item of items) {
    const contexts = item.fusion_contexts?.length
      ? item.fusion_contexts
      : [null];
    for (const context of contexts) {
      rows.push([
        item.gene,
        context?.fusion,
        context?.side,
        context?.partner_gene,
        context?.five_transcript,
        context?.three_transcript,
        context?.five_exon,
        context?.three_exon,
        context?.five_genomic,
        context?.three_genomic,
        context?.five_protein_breakpoint,
        context?.three_protein_breakpoint,
        context?.retained_domains,
        [...(context?.lost_domains || []), ...(context?.disrupted_domains || [])],
        context?.kinase_gene,
        context?.kinase_gene_side,
        context?.kinase_domain_status,
        item.cancer_associated == null ? "" : item.cancer_associated ? "TRUE" : "FALSE",
        item.rationale,
        item.supporting_pmids,
        item.retrieved_pmids,
        item.insufficient_evidence ? "TRUE" : "FALSE",
        item.error || context?.annotation_error,
      ].map(csvCell).join(","));
    }
  }
  return `${rows.join("\n")}\n`;
}

function downloadText(filename: string, text: string, type: string) {
  const blob = new Blob([text], { type });
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = filename;
  document.body.appendChild(link);
  link.click();
  link.remove();
  URL.revokeObjectURL(url);
}

function App() {
  const [formValues, setFormValues] = useState<AnnotateParams>(paramsFromLocation);
  const [result, setResult] = useState<AnnotationResult | null>(null);
  const [batchResults, setBatchResults] = useState<BatchAnnotationItemResult[] | null>(null);
  const [geneCurationStatus, setGeneCurationStatus] = useState<GeneCurationStatus | null>(null);
  const [geneCurationResults, setGeneCurationResults] = useState<GeneCurationGeneResult[] | null>(null);
  const [derived, setDerived] = useState<DerivedInputs | null>(null);
  const [error, setError] = useState<ApiError | null>(null);
  const [curationError, setCurationError] = useState<ApiError | null>(null);
  const [loading, setLoading] = useState(false);
  const [batchLoading, setBatchLoading] = useState(false);
  const [curationLoading, setCurationLoading] = useState(false);
  const [permalink, setPermalink] = useState(window.location.href);
  const requestSequence = useRef(0);
  const curationRequestSequence = useRef(0);
  const activeRequest = useRef<AbortController | null>(null);
  const activeCurationRequest = useRef<AbortController | null>(null);

  const runAnnotation = useCallback(async (params: AnnotateParams, shouldPushState = true) => {
    activeRequest.current?.abort();
    const controller = new AbortController();
    activeRequest.current = controller;
    const requestId = requestSequence.current + 1;
    requestSequence.current = requestId;

    setFormValues(params);
    setLoading(true);
    setError(null);
    setResult(null);
    setBatchResults(null);
    setBatchLoading(false);
    setGeneCurationResults(null);
    setCurationError(null);
    setDerived(null);

    if (shouldPushState) {
      // Write inputs to the URL query string *before* the request resolves,
      // so the address bar is always a valid, shareable permalink for this
      // lookup even if the caller navigates away before it finishes. Skipped
      // when we're replaying a URL that already reflects the intended state
      // (initial load, or a popstate from Back/Forward) — pushing again
      // there would create a duplicate history entry and trap the user.
      const search = toSearchParams(params);
      const newUrl = `${window.location.pathname}?${search.toString()}`;
      window.history.pushState(null, "", newUrl);
      setPermalink(window.location.href);
    }

    try {
      const annotation = await annotateFusion(params, controller.signal);
      if (requestSequence.current !== requestId) return;
      setResult(annotation);
      setDerived(computeDerivedInputs(annotation));
    } catch (err) {
      if (controller.signal.aborted || requestSequence.current !== requestId) return;
      setResult(null);
      setError(err as ApiError);
    } finally {
      if (activeRequest.current === controller) {
        activeRequest.current = null;
      }
      if (requestSequence.current === requestId) {
        setLoading(false);
      }
    }
  }, []);

  const runBatchAnnotation = useCallback(async (fusions: AnnotateParams[]) => {
    activeRequest.current?.abort();
    const controller = new AbortController();
    activeRequest.current = controller;
    const requestId = requestSequence.current + 1;
    requestSequence.current = requestId;

    setBatchLoading(true);
    setLoading(false);
    setError(null);
    setResult(null);
    setDerived(null);
    setBatchResults(null);
    setGeneCurationResults(null);
    setCurationError(null);

    try {
      const response = await annotateFusionBatch(fusions, controller.signal);
      if (requestSequence.current !== requestId) return;
      setBatchResults(response.results);
    } catch (err) {
      if (controller.signal.aborted || requestSequence.current !== requestId) return;
      setError(err as ApiError);
    } finally {
      if (activeRequest.current === controller) {
        activeRequest.current = null;
      }
      if (requestSequence.current === requestId) {
        setBatchLoading(false);
      }
    }
  }, []);

  const runGeneCuration = useCallback(async (fusions: AnnotateParams[]) => {
    activeCurationRequest.current?.abort();
    const controller = new AbortController();
    activeCurationRequest.current = controller;
    const requestId = curationRequestSequence.current + 1;
    curationRequestSequence.current = requestId;

    setCurationLoading(true);
    setCurationError(null);
    setGeneCurationResults(null);

    try {
      const response = await curateFusionGenes(fusions, controller.signal);
      if (curationRequestSequence.current !== requestId) return;
      setGeneCurationResults(response.genes);
    } catch (err) {
      if (controller.signal.aborted || curationRequestSequence.current !== requestId) return;
      setCurationError(err as ApiError);
    } finally {
      if (activeCurationRequest.current === controller) {
        activeCurationRequest.current = null;
      }
      if (curationRequestSequence.current === requestId) {
        setCurationLoading(false);
      }
    }
  }, []);

  const exportGeneCurationCsv = useCallback(() => {
    if (!geneCurationResults?.length) return;
    downloadText(
      "fusion_gene_curation.csv",
      curationCsvRows(geneCurationResults),
      "text/csv;charset=utf-8",
    );
  }, [geneCurationResults]);

  // On load (including a permalink being opened fresh), auto-run the
  // annotation if the URL already carries both gene names.
  useEffect(() => {
    const initial = paramsFromLocation();
    if (initial.five_gene && initial.three_gene) {
      runAnnotation(initial, false);
    }
    // Also react to browser back/forward between permalinks.
    const onPopState = () => {
      const params = paramsFromLocation();
      setFormValues(params);
      setPermalink(window.location.href);
      if (params.five_gene && params.three_gene) {
        runAnnotation(params, false);
      } else {
        setResult(null);
      }
    };
    window.addEventListener("popstate", onPopState);
    return () => {
      activeRequest.current?.abort();
      activeCurationRequest.current?.abort();
      window.removeEventListener("popstate", onPopState);
    };
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  useEffect(() => {
    const controller = new AbortController();
    getGeneCurationStatus(controller.signal)
      .then(setGeneCurationStatus)
      .catch(() => {
        setGeneCurationStatus({ enabled: false, model: "" });
      });
    return () => controller.abort();
  }, []);

  return (
    <div className="app">
      <header>
        <h1>fusion-annotation</h1>
        <p className="tagline">
          Protein-level gene-fusion annotation — reading frame, junction, domain retention, and curated
          clinical knowledge.
        </p>
      </header>

      <ExampleFusions onSelect={runAnnotation} disabled={loading} />

      <FusionForm initial={formValues} derived={derived} onSubmit={runAnnotation} loading={loading} />

      <BatchFusionForm
        genomeBuild={formValues.genome_build}
        onSubmit={runBatchAnnotation}
        loading={batchLoading}
        onCurate={runGeneCuration}
        curationLoading={curationLoading}
        curationEnabled={geneCurationStatus?.enabled ?? false}
      />

      {geneCurationStatus && !geneCurationStatus.enabled && (
        <div className="notice-box">
          Server-side gene curation is not configured for this deployment.
        </div>
      )}

      {error && (
        <div className="error-box">
          <strong>Error {error.status}:</strong> {error.detail}
        </div>
      )}

      {curationError && (
        <div className="error-box">
          <strong>Curation error {curationError.status}:</strong> {curationError.detail}
        </div>
      )}

      {result && <ResultView result={result} permalink={permalink} />}

      {batchResults && (
        <section className="batch-results">
          <h2>Batch results</h2>
          {batchResults.map((item, index) => {
            const label = `${item.input.five_gene}::${item.input.three_gene}`;
            return (
              <article className="batch-result-item" key={`${label}-${index}`}>
                <h3>{label}</h3>
                {item.error && <div className="error-box">{item.error}</div>}
                {item.result && (
                  <ResultView
                    result={item.result}
                    permalink={`${window.location.origin}${window.location.pathname}?${toSearchParams(item.input).toString()}`}
                  />
                )}
              </article>
            );
          })}
        </section>
      )}

      {geneCurationResults && (
        <section className="gene-curation-results">
          <div className="gene-curation-results-header">
            <div>
              <h2>Gene literature curation</h2>
              <p className="gene-curation-note">
                Server-side curation uses Genome Nexus fusion structure and PubMed abstracts for the batch.
              </p>
            </div>
            <button type="button" className="secondary-button" onClick={exportGeneCurationCsv}>
              Export curation CSV
            </button>
          </div>
          <div className="gene-curation-grid">
            {geneCurationResults.map((item) => (
              <article className="gene-curation-card" key={item.gene}>
                <div className="gene-curation-card-header">
                  <h3>{item.gene}</h3>
                  <div className="curation-badges" aria-label={`${item.gene} review signals`}>
                    {curationBadges(item).map((badge) => (
                      <span
                        className={`status-chip ${badge.tone}`}
                        key={badge.label}
                        title={badge.title}
                      >
                        {badge.label}
                      </span>
                    ))}
                  </div>
                </div>
                {item.error ? (
                  <div className="error-box">{item.error}</div>
                ) : (
                  <>
                    {item.fusion_contexts && item.fusion_contexts.length > 0 && (
                      <div className="fusion-curation-contexts">
                        {item.fusion_contexts.map(renderFusionContext)}
                      </div>
                    )}
                    <dl className="gene-curation-fields">
                      <dt>Cancer associated</dt>
                      <dd>{item.cancer_associated == null ? "Unknown" : item.cancer_associated ? "Yes" : "No"}</dd>
                      <dt>Rationale</dt>
                      <dd>{item.rationale || "No rationale returned."}</dd>
                    </dl>
                    <div className="pmid-row">
                      <strong>Supporting PMIDs</strong>
                      <span>{item.supporting_pmids?.join(", ") || "None selected"}</span>
                    </div>
                    <div className="pmid-row">
                      <strong>Retrieved PMIDs</strong>
                      <span>{item.retrieved_pmids?.join(", ") || "None retrieved"}</span>
                    </div>
                  </>
                )}
              </article>
            ))}
          </div>
        </section>
      )}

      <VersionFootnote />
    </div>
  );
}

export default App;
