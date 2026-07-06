# Copilot instructions for fusion-annotation

## Build, test, and lint commands

```bash
# Core package + tests
pip install -e .
pip install -e ".[test]"
pytest -v
pytest tests/test_eml4_alk.py::test_hgvsp_string -q

# REST API tests (needs api extras)
pip install -e ".[test,api]"
pytest tests/test_api_app.py::test_annotate_get_eml4_alk -q

# MCP server tests (needs server extras)
pip install -e ".[test,server]"
pytest tests/test_server_app.py::test_annotate_gene_fusion_include_diagram_false_skips_image -q

# Python lint
pip install -e ".[dev]"
ruff check .

# Web UI
cd web && npm ci
cd web && npm run lint
cd web && npm run build

# Regenerate the reference domain diagram used in docs/server parity work
pip install -e ".[docs]"
python docs/generate_domain_map.py
```

## High-level architecture

- `src/fusion_annotation/core.py` is the center of the system. `annotate_fusion()` orchestrates three layers: a VEP-like effect engine, an HGVS.p-like interface object, and an OncoKB-like knowledge layer. The canonical return shape is `{"interface", "knowledge", "resolved", "warnings"}`.
- Provider implementations are the main integration boundary. `GenomeNexusDataProvider` is the default fast human backend, `RestDataProvider` is the slower direct Ensembl/InterPro/CIViC fallback, and `MCPDataProvider` is only for runtimes that inject an MCP host callable.
- `src/fusion_annotation/provider_factory.py` is the shared provider-selection path for both HTTP-facing entry points. Keep backend selection rules in that module so `api/app.py` and `server/app.py` stay behaviorally aligned.
- `server/app.py` exposes the MCP tool over FastMCP and may attach a rendered PNG domain diagram in addition to the structured JSON result.
- `api/app.py` exposes the same annotation engine over `GET /api/annotate` and `POST /api/annotate`; the GET query string is intentionally a stateless permalink, so reopening the URL reruns annotation instead of loading stored state.
- `web/` is a Vite React SPA that mirrors the API/MCP schema closely. `web/src/lib/types.ts` matches the JSON returned by `annotate_fusion()`, and the current lookup is always encoded into the browser URL query string.
- `tests/fixtures/*.json` plus `StaticProvider` underpin the offline examples and most tests. The test suite prefers fixture-backed providers and monkeypatching over live network calls.

## Key conventions

- Keep the core package transport-agnostic and dependency-light. `src/fusion_annotation/` should not pick up FastAPI, Starlette, or MCP framework concerns unless the change is explicitly transport-specific.
- When changing annotation inputs or outputs, update all mirrored surfaces together: `core.py`, the MCP tool schema in `mcp_tool.py`, `api/app.py` request/response models, and `web/src/lib/types.ts` / form handling.
- Genomic breakpoints are preferred over exon-only breakpoints when both are available. The code intentionally echoes resolved transcripts and breakpoint provenance under `resolved` rather than hiding isoform ambiguity.
- Known oncogenic partner pairs that reconstruct out-of-frame should produce warnings, not silent fallbacks. Preserve the existing behavior that surfaces likely transcript/exon mismatch to the caller.
- Domain rendering has cross-surface parity requirements. The Python renderer used by docs/server and the TypeScript logic in `web/src/lib/domainDiagram.ts` use the same stable color mapping and domain de-duplication rules; keep them in sync.
- The default production path is Genome Nexus/UCSC/CIViC for human input. Non-human species and `FUSION_ANNOTATION_PROVIDER=rest` intentionally route through `RestDataProvider` instead.
- The public API and MCP server are intentionally stateless. Preserve permalink behavior, environment-driven provider selection, rate limiting / CORS / allowed-hosts controls, and fixture-backed offline tests.
- The web build expects `VITE_API_BASE_URL` at build time, and production assets are served from the GitHub Pages project-site base path `/fusion-annotation/`.
