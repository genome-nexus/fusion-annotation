import type { AnnotateParams, AnnotationResult, ApiError } from "./types";

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

export async function annotateFusion(params: AnnotateParams): Promise<AnnotationResult> {
  const search = toSearchParams(params);
  const response = await fetch(`${API_BASE_URL}/api/annotate?${search.toString()}`);
  if (!response.ok) {
    let detail = response.statusText;
    try {
      const body = await response.json();
      detail = body.detail ?? detail;
    } catch {
      // response body wasn't JSON — fall back to statusText
    }
    const error: ApiError = { status: response.status, detail };
    throw error;
  }
  return response.json();
}

export { toSearchParams };
