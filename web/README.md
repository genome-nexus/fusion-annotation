# fusion-annotation web UI

A React + TypeScript SPA (Vite) for looking up a gene fusion in a browser,
backed by the REST API in [`../api`](../api). See the root
[README](../README.md#web-ui--rest-api) for how this fits together with the
MCP server.

## Local development

Start the API first (from the repo root):

```bash
pip install -e ".[api]"
python api/app.py          # serves on :8080
```

Then, in this directory:

```bash
npm install
npm run dev                # http://localhost:5173, proxies /api to :8080
```

`vite.config.ts` proxies `/api` and `/health` to `localhost:8080` in dev, so
`VITE_API_BASE_URL` can stay unset locally.

## Configuration

`VITE_API_BASE_URL` — base URL of the deployed API (e.g. a Cloud Run URL).
Vite inlines env vars into the built JS bundle at **build** time, so this
must be set before `npm run build`, not at container/hosting runtime (there
is no runtime container for the web UI — see below). See `.env.example`.

## Build & deploy

```bash
npm run build               # outputs to dist/ (base path /fusion-annotation/,
                            #   matching GitHub Pages project-site serving)
```

Deployed to **GitHub Pages** (static hosting, no container) at
`https://genome-nexus.github.io/fusion-annotation/`, via
[`deploy-on-release.yml`](../.github/workflows/deploy-on-release.yml)
(automatic on every published release) or
[`deploy-api-web.yml`](../.github/workflows/deploy-api-web.yml) (manual
trigger). Both build with `VITE_API_BASE_URL` set to the just-deployed API's
Cloud Run URL, then publish `dist/` to the `gh-pages` branch (Pages' classic
branch-based source — Settings → Pages → Source: "Deploy from a branch").
That means the live site always reflects exactly the commit that was last
released or manually redeployed, decoupled from wherever `main` has moved
on to since. `vite.config.ts` sets `base: '/fusion-annotation/'` for
production builds specifically because Pages serves a project site under
that subpath rather than at the domain root.

## Permalinks

The current lookup's inputs are always synced to the URL's query string
(`five_gene`, `three_gene`, `five_exon`, ... — see `src/lib/types.ts`), so
the address bar is a stateless, shareable permalink: opening it re-runs the
annotation against the live API rather than fetching a stored result.
