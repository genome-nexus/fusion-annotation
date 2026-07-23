import type {
  AnnotateParams,
  AnnotationResult,
  ApiError,
  BatchAnnotationResponse,
  GeneCurationResponse,
  GeneCurationStatus,
} from "./types";

// Base URL of the deployed REST API (api/app.py). In dev, Vite's proxy
// (vite.config.ts) forwards /api to a locally-running API, so this can stay
// empty; in production, set VITE_API_BASE_URL to the deployed Cloud Run URL.
const API_BASE_URL = import.meta.env.VITE_API_BASE_URL ?? "";

/** Strip empty-string/undefined fields so the query string (and thus the
 * permalink URL) only ever contains the fields the user actually filled in. */
function toSearchParams(params: AnnotateParams): URLSearchParams {
  const search = new URLSearchParams();
  for (const [key, value] of Object.entries(params)) {
    if (value !== undefined && value !== null && value !== "") {
      search.set(key, String(value));
    }
  }
  return search;
}

function cleanAnnotateParams(params: AnnotateParams) {
  const cleaned: Record<string, string | number | boolean> = {};
  for (const [key, value] of Object.entries(params)) {
    if (value === undefined || value === null || value === "") continue;
    if ((key === "five_exon" || key === "three_exon") && typeof value === "string") {
      cleaned[key] = Number(value);
    } else {
      cleaned[key] = value;
    }
  }
  return cleaned;
}

/** FastAPI/Pydantic validation errors (HTTP 422) return `detail` as a list of
 * {loc, msg, type} objects rather than a plain string; every other error
 * path in api/app.py (HTTPException) returns a string. Normalize both into a
 * single human-readable string so callers can render `error.detail` directly
 * without risking a "objects are not valid as a React child" crash. */
function formatDetail(detail: unknown, fallback: string): string {
  if (typeof detail === "string") return detail;
  if (Array.isArray(detail)) {
    return detail
      .map((err) => {
        const loc = Array.isArray(err?.loc) ? err.loc.join(".") : undefined;
        const msg = typeof err?.msg === "string" ? err.msg : JSON.stringify(err);
        return loc ? `${loc}: ${msg}` : msg;
      })
      .join(", ");
  }
  if (detail != null) return JSON.stringify(detail);
  return fallback;
}

export async function annotateFusion(
  params: AnnotateParams,
  signal?: AbortSignal,
): Promise<AnnotationResult> {
  const search = toSearchParams(params);
  const response = await fetch(`${API_BASE_URL}/api/annotate?${search.toString()}`, { signal });
  if (!response.ok) {
    let detail = response.statusText;
    try {
      const body = await response.json();
      detail = formatDetail(body?.detail, detail);
    } catch {
      // response body wasn't JSON — fall back to statusText
    }
    const error: ApiError = { status: response.status, detail };
    throw error;
  }
  return response.json();
}

export async function annotateFusionBatch(
  fusions: AnnotateParams[],
  signal?: AbortSignal,
): Promise<BatchAnnotationResponse> {
  const response = await fetch(`${API_BASE_URL}/api/annotate/batch`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ fusions: fusions.map(cleanAnnotateParams) }),
    signal,
  });
  if (!response.ok) {
    let detail = response.statusText;
    try {
      const body = await response.json();
      detail = formatDetail(body?.detail, detail);
    } catch {
      // response body wasn't JSON — fall back to statusText
    }
    const error: ApiError = { status: response.status, detail };
    throw error;
  }
  return response.json();
}

export async function getGeneCurationStatus(signal?: AbortSignal): Promise<GeneCurationStatus> {
  const response = await fetch(`${API_BASE_URL}/api/gene-curation/status`, { signal });
  if (!response.ok) {
    let detail = response.statusText;
    try {
      const body = await response.json();
      detail = formatDetail(body?.detail, detail);
    } catch {
      // response body wasn't JSON — fall back to statusText
    }
    const error: ApiError = { status: response.status, detail };
    throw error;
  }
  return response.json();
}

export async function curateFusionGenes(
  fusions: AnnotateParams[],
  forceGeneCuration = false,
  genes?: string[],
  tumorType?: string,
  signal?: AbortSignal,
): Promise<GeneCurationResponse> {
  const response = await fetch(`${API_BASE_URL}/api/gene-curation`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({
      fusions: fusions.map(cleanAnnotateParams),
      force_gene_curation: forceGeneCuration,
      genes: genes ?? [],
      tumor_type: tumorType || undefined,
    }),
    signal,
  });
  if (!response.ok) {
    let detail = response.statusText;
    try {
      const body = await response.json();
      detail = formatDetail(body?.detail, detail);
    } catch {
      // response body wasn't JSON — fall back to statusText
    }
    const error: ApiError = { status: response.status, detail };
    throw error;
  }
  return response.json();
}

export { toSearchParams };
