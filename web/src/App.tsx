import { useCallback, useEffect, useRef, useState } from "react";
import "./App.css";
import { BatchFusionForm } from "./components/BatchFusionForm";
import { ExampleFusions } from "./components/ExampleFusions";
import { FusionForm } from "./components/FusionForm";
import { ResultView } from "./components/ResultView";
import { VersionFootnote } from "./components/VersionFootnote";
import { annotateFusion, annotateFusionBatch, toSearchParams } from "./lib/api";
import { DEFAULT_PARAMS } from "./lib/defaultParams";
import { computeDerivedInputs, type DerivedInputs } from "./lib/derivedInputs";
import type { AnnotateParams, AnnotationResult, ApiError, BatchAnnotationItemResult } from "./lib/types";

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

function App() {
  const [formValues, setFormValues] = useState<AnnotateParams>(paramsFromLocation);
  const [result, setResult] = useState<AnnotationResult | null>(null);
  const [batchResults, setBatchResults] = useState<BatchAnnotationItemResult[] | null>(null);
  const [derived, setDerived] = useState<DerivedInputs | null>(null);
  const [error, setError] = useState<ApiError | null>(null);
  const [loading, setLoading] = useState(false);
  const [batchLoading, setBatchLoading] = useState(false);
  const [permalink, setPermalink] = useState(window.location.href);
  const requestSequence = useRef(0);
  const activeRequest = useRef<AbortController | null>(null);

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

    try {
      const response = await annotateFusionBatch(fusions, controller.signal);
      if (requestSequence.current !== requestId) return;
      setBatchResults(response.results);
      const firstSuccessful = response.results.find((item) => item.result)?.result;
      setDerived(firstSuccessful ? computeDerivedInputs(firstSuccessful) : null);
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
      window.removeEventListener("popstate", onPopState);
    };
    // eslint-disable-next-line react-hooks/exhaustive-deps
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
      />

      {error && (
        <div className="error-box">
          <strong>Error {error.status}:</strong> {error.detail}
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

      <VersionFootnote />
    </div>
  );
}

export default App;
