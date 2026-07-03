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
must be set before `npm run build`/`docker build`, not at container runtime.
See `.env.example`.

## Build & deploy

```bash
npm run build               # outputs to dist/
```

`Dockerfile` is a multi-stage build: compiles the SPA with `VITE_API_BASE_URL`
as a build arg, then serves the static output via nginx (`nginx.conf`) on
port 8080. Deploy to Cloud Run with:

```bash
export GCP_PROJECT=your-project-id
export API_URL=https://<deployed-api>.run.app
./deploy.sh
```

## Permalinks

The current lookup's inputs are always synced to the URL's query string
(`five_gene`, `three_gene`, `five_exon`, ... — see `src/lib/types.ts`), so
the address bar is a stateless, shareable permalink: opening it re-runs the
annotation against the live API rather than fetching a stored result.
