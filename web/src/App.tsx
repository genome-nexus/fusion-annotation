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
  GeneCurationFusionResult,
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

const CURATION_CSV_HEADERS = [
  "row_type",
  "gene",
  "fusion",
  "fusion_literature_identified",
  "curation_source",
  "oncokb_gene_type",
  "oncokb_highest_sensitive_level",
  "oncokb_highest_resistance_level",
  "oncokb_url",
  "gene_side",
  "partner_gene",
  "breakpoint_context_available",
  "fusion_specificity",
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
  "limitations",
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

function curationCsvRows(
  fusionItems: GeneCurationFusionResult[],
  geneItems: GeneCurationGeneResult[],
) {
  const rows = [CURATION_CSV_HEADERS.map(csvCell).join(",")];
  for (const item of fusionItems) {
    const contexts = item.fusion_contexts?.length
      ? item.fusion_contexts
      : [null];
    for (const context of contexts) {
      rows.push([
        "fusion",
        context?.gene,
        item.fusion,
        item.fusion_literature_identified == null
          ? ""
          : item.fusion_literature_identified ? "TRUE" : "FALSE",
        "PubMed + LLM",
        "",
        "",
        "",
        "",
        context?.side,
        context?.partner_gene,
        context?.breakpoint_context_available == null
          ? ""
          : context.breakpoint_context_available ? "TRUE" : "FALSE",
        context?.fusion_specificity,
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
        context?.limitations,
        item.cancer_associated == null ? "" : item.cancer_associated ? "TRUE" : "FALSE",
        item.rationale,
        item.supporting_pmids,
        item.retrieved_pmids,
        item.insufficient_evidence ? "TRUE" : "FALSE",
        item.error || context?.annotation_error,
      ].map(csvCell).join(","));
    }
  }
  for (const item of geneItems) {
    const contexts = item.fusion_contexts?.length
      ? item.fusion_contexts
      : [null];
    for (const context of contexts) {
      rows.push([
        "gene",
        item.gene,
        context?.fusion,
        "",
        item.curation_source,
        item.oncokb_gene_type,
        item.oncokb_highest_sensitive_level,
        item.oncokb_highest_resistance_level,
        item.oncokb_url,
        context?.side,
        context?.partner_gene,
        context?.breakpoint_context_available == null
          ? ""
          : context.breakpoint_context_available ? "TRUE" : "FALSE",
        context?.fusion_specificity,
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
        context?.limitations,
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

function mergeFusionCurationResults(
  current: GeneCurationFusionResult[] | null,
  incoming: GeneCurationFusionResult[],
) {
  const byFusion = new Map((current || []).map((item) => [item.fusion, item]));
  for (const item of incoming) {
    byFusion.set(item.fusion, item);
  }
  return Array.from(byFusion.values());
}

function contextKey(context: NonNullable<GeneCurationGeneResult["fusion_contexts"]>[number]) {
  return `${context.fusion}|${context.side}|${context.partner_gene}`;
}

function mergeGeneCurationResults(
  current: GeneCurationGeneResult[] | null,
  incoming: GeneCurationGeneResult[],
) {
  const byGene = new Map((current || []).map((item) => [item.gene.toUpperCase(), item]));
  for (const item of incoming) {
    const key = item.gene.toUpperCase();
    const existing = byGene.get(key);
    if (!existing) {
      byGene.set(key, item);
      continue;
    }
    const contexts = new Map(
      (existing.fusion_contexts || []).map((context) => [contextKey(context), context]),
    );
    for (const context of item.fusion_contexts || []) {
      contexts.set(contextKey(context), context);
    }
    byGene.set(key, {
      ...existing,
      ...item,
      fusion_contexts: Array.from(contexts.values()),
    });
  }
  return Array.from(byGene.values());
}

type AppTab = "single" | "batch";

function App() {
  const [activeTab, setActiveTab] = useState<AppTab>("single");
  const [formValues, setFormValues] = useState<AnnotateParams>(paramsFromLocation);
  const [result, setResult] = useState<AnnotationResult | null>(null);
  const [batchResults, setBatchResults] = useState<BatchAnnotationItemResult[] | null>(null);
  const [selectedBatchIndex, setSelectedBatchIndex] = useState(0);
  const [geneCurationStatus, setGeneCurationStatus] = useState<GeneCurationStatus | null>(null);
  const [fusionCurationResults, setFusionCurationResults] = useState<GeneCurationFusionResult[] | null>(null);
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
    curationRequestSequence.current += 1;
    activeCurationRequest.current?.abort();
    activeCurationRequest.current = null;

    setFormValues(params);
    setLoading(true);
    setError(null);
    setResult(null);
    setBatchResults(null);
    setSelectedBatchIndex(0);
    setBatchLoading(false);
    setFusionCurationResults(null);
    setGeneCurationResults(null);
    setCurationError(null);
    setCurationLoading(false);
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
    curationRequestSequence.current += 1;
    activeCurationRequest.current?.abort();
    activeCurationRequest.current = null;

    setBatchLoading(true);
    setLoading(false);
    setError(null);
    setResult(null);
    setDerived(null);
    setBatchResults(null);
    setSelectedBatchIndex(0);
    setFusionCurationResults(null);
    setGeneCurationResults(null);
    setCurationError(null);
    setCurationLoading(false);

    try {
      const response = await annotateFusionBatch(fusions, controller.signal);
      if (requestSequence.current !== requestId) return;
      setBatchResults(response.results);
      setSelectedBatchIndex(Math.max(0, response.results.findIndex((item) => item.result)));
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

  const runGeneCuration = useCallback(async (
    fusions: AnnotateParams[],
    forceGeneCuration = false,
    genes?: string[],
  ) => {
    activeCurationRequest.current?.abort();
    const controller = new AbortController();
    activeCurationRequest.current = controller;
    const requestId = curationRequestSequence.current + 1;
    curationRequestSequence.current = requestId;
    const shouldMerge = Boolean(genes?.length);

    setCurationLoading(true);
    setCurationError(null);
    if (!shouldMerge) {
      setFusionCurationResults(null);
      setGeneCurationResults(null);
    }

    try {
      const response = await curateFusionGenes(fusions, forceGeneCuration, genes, controller.signal);
      if (curationRequestSequence.current !== requestId) return;
      if (shouldMerge) {
        setFusionCurationResults((current) => mergeFusionCurationResults(current, response.fusions || []));
        setGeneCurationResults((current) => mergeGeneCurationResults(current, response.genes));
      } else {
        setFusionCurationResults(response.fusions || []);
        setGeneCurationResults(response.genes);
      }
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
    if (!geneCurationResults?.length && !fusionCurationResults?.length) return;
    downloadText(
      "fusion_gene_curation.csv",
      curationCsvRows(fusionCurationResults || [], geneCurationResults || []),
      "text/csv;charset=utf-8",
    );
  }, [fusionCurationResults, geneCurationResults]);

  const selectedBatchItem = batchResults?.[selectedBatchIndex] ?? null;
  const batchInputs = batchResults?.map((item) => item.input) ?? [];

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

      <nav className="workflow-tabs" aria-label="Annotation workflow">
        <button
          type="button"
          className={activeTab === "single" ? "workflow-tab active" : "workflow-tab"}
          aria-pressed={activeTab === "single"}
          onClick={() => setActiveTab("single")}
        >
          Single fusion
        </button>
        <button
          type="button"
          className={activeTab === "batch" ? "workflow-tab active" : "workflow-tab"}
          aria-pressed={activeTab === "batch"}
          onClick={() => setActiveTab("batch")}
        >
          Batch annotation
        </button>
      </nav>

      {activeTab === "single" && (
        <>
          <ExampleFusions onSelect={runAnnotation} disabled={loading} />

          <FusionForm initial={formValues} derived={derived} onSubmit={runAnnotation} loading={loading} />
        </>
      )}

      {activeTab === "batch" && (
        <BatchFusionForm
          genomeBuild={formValues.genome_build}
          onSubmit={runBatchAnnotation}
          loading={batchLoading}
          onCurate={runGeneCuration}
          curationLoading={curationLoading}
          curationEnabled={geneCurationStatus?.enabled ?? false}
        />
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

      {activeTab === "single" && result && (
        <ResultView
          result={result}
          permalink={permalink}
          fusionCurationResults={fusionCurationResults}
          geneCurationResults={geneCurationResults}
          geneCurationLoading={curationLoading}
          geneCurationEnabled={geneCurationStatus?.enabled ?? true}
          geneCurationError={curationError}
          onCurateGenes={() => runGeneCuration([formValues])}
          onForceGeneCuration={() => runGeneCuration([formValues], true)}
          onCurateGene={(gene) => runGeneCuration([formValues], true, [gene])}
          onExportGeneCurationCsv={exportGeneCurationCsv}
        />
      )}

      {activeTab === "batch" && batchResults && (
        <section className="batch-results">
          <div className="batch-results-header">
            <div>
              <h2>Batch results</h2>
              <p>Select a completed run to inspect it with the same annotation view used for a single fusion.</p>
            </div>
            {geneCurationResults?.length || fusionCurationResults?.length ? (
              <button type="button" className="secondary-button" onClick={exportGeneCurationCsv}>
                Export curation CSV
              </button>
            ) : null}
          </div>
          <div className="batch-review-layout">
            <div className="batch-run-list" role="list" aria-label="Completed batch runs">
              {batchResults.map((item, index) => {
                const label = `${item.input.five_gene}::${item.input.three_gene}`;
                return (
                  <button
                    type="button"
                    role="listitem"
                    className={index === selectedBatchIndex ? "batch-run-item active" : "batch-run-item"}
                    key={`${label}-${index}`}
                    onClick={() => setSelectedBatchIndex(index)}
                  >
                    <span>{label}</span>
                    <small>{item.result ? "Annotated" : "Needs input review"}</small>
                  </button>
                );
              })}
            </div>
            <div className="batch-review-pane">
              {selectedBatchItem?.error && <div className="error-box">{selectedBatchItem.error}</div>}
              {selectedBatchItem?.result && (
                <ResultView
                  result={selectedBatchItem.result}
                  permalink={`${window.location.origin}${window.location.pathname}?${toSearchParams(selectedBatchItem.input).toString()}`}
                  fusionCurationResults={fusionCurationResults}
                  geneCurationResults={geneCurationResults}
                  geneCurationLoading={curationLoading}
                  geneCurationEnabled={geneCurationStatus?.enabled ?? true}
                  geneCurationError={curationError}
                  onCurateGenes={() => runGeneCuration(batchInputs)}
                  onForceGeneCuration={() => runGeneCuration(batchInputs, true)}
                  onCurateGene={(gene) => runGeneCuration([selectedBatchItem.input], true, [gene])}
                  onExportGeneCurationCsv={exportGeneCurationCsv}
                />
              )}
              {!selectedBatchItem && (
                <div className="notice-box">Run a batch to review individual fusion annotations here.</div>
              )}
            </div>
          </div>
        </section>
      )}

      <VersionFootnote />
    </div>
  );
}

export default App;
