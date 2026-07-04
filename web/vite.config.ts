import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'

// https://vite.dev/config/
export default defineConfig(({ command }) => ({
  plugins: [react()],
  // GitHub Pages serves this as a project site at
  // https://genome-nexus.github.io/fusion-annotation/ (a subpath, not root),
  // so production asset URLs must be prefixed accordingly. The dev server
  // still runs at "/" so local `npm run dev` is unaffected.
  base: command === 'build' ? '/fusion-annotation/' : '/',
  server: {
    // In dev, forward /api and /health to a locally-running api/app.py
    // (uvicorn api/app.py, default port 8080) so VITE_API_BASE_URL can stay
    // unset locally; set it explicitly for a deployed API in production.
    // Override the target with DEV_API_PROXY_TARGET if 8080 is already
    // taken by something else on your machine (start the API with a
    // matching PORT env var).
    proxy: {
      '/api': process.env.DEV_API_PROXY_TARGET ?? 'http://localhost:8080',
      '/health': process.env.DEV_API_PROXY_TARGET ?? 'http://localhost:8080',
    },
  },
}))
