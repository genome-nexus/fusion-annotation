import { useCallback, useEffect, useState } from "react";
import "./App.css";
import { FusionForm } from "./components/FusionForm";
import { ResultView } from "./components/ResultView";
import { annotateFusion, toSearchParams } from "./lib/api";
import type { AnnotateParams, AnnotationResult, ApiError } from "./lib/types";

const DEFAULT_PARAMS: AnnotateParams = {
  five_gene: "",
  three_gene: "",
  five_exon: "",
  three_exon: "",
  five_genomic: "",
  three_genomic: "",
  five_transcript: "",
  three_transcript: "",
  genome_build: "GRCh38",
  species: "homo_sapiens",
};

/** Read the current URL's query string into AnnotateParams — the other half
 * of the stateless-permalink contract (writing happens in runAnnotation
 * below via history.pushState). No fields are persisted server-side: this
 * URL alone is enough to reproduce the annotation. */
function paramsFromLocation(): AnnotateParams {
  const search = new URLSearchParams(window.location.search);
  const params = { ...DEFAULT_PARAMS };
  for (const key of Object.keys(params) as (keyof AnnotateParams)[]) {
    const value = search.get(key);
    if (value !== null) params[key] = value;
  }
  return params;
}

function App() {
  const [formValues, setFormValues] = useState<AnnotateParams>(paramsFromLocation);
  const [result, setResult] = useState<AnnotationResult | null>(null);
  const [error, setError] = useState<ApiError | null>(null);
  const [loading, setLoading] = useState(false);
  const [permalink, setPermalink] = useState(window.location.href);

  const runAnnotation = useCallback(async (params: AnnotateParams) => {
    setFormValues(params);
    setLoading(true);
    setError(null);

    // Write inputs to the URL query string *before* the request resolves, so
    // the address bar is always a valid, shareable permalink for this
    // lookup even if the caller navigates away before it finishes.
    const search = toSearchParams(params);
    const newUrl = `${window.location.pathname}?${search.toString()}`;
    window.history.pushState(null, "", newUrl);
    setPermalink(window.location.href);

    try {
      const annotation = await annotateFusion(params);
      setResult(annotation);
    } catch (err) {
      setResult(null);
      setError(err as ApiError);
    } finally {
      setLoading(false);
    }
  }, []);

  // On load (including a permalink being opened fresh), auto-run the
  // annotation if the URL already carries both gene names.
  useEffect(() => {
    const initial = paramsFromLocation();
    if (initial.five_gene && initial.three_gene) {
      runAnnotation(initial);
    }
    // Also react to browser back/forward between permalinks.
    const onPopState = () => {
      const params = paramsFromLocation();
      setFormValues(params);
      if (params.five_gene && params.three_gene) {
        runAnnotation(params);
      } else {
        setResult(null);
      }
    };
    window.addEventListener("popstate", onPopState);
    return () => window.removeEventListener("popstate", onPopState);
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

      <FusionForm initial={formValues} onSubmit={runAnnotation} loading={loading} />

      {error && (
        <div className="error-box">
          <strong>Error {error.status}:</strong> {error.detail}
        </div>
      )}

      {result && <ResultView result={result} permalink={permalink} />}
    </div>
  );
}

export default App;
