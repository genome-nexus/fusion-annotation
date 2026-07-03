import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'

// https://vite.dev/config/
export default defineConfig({
  plugins: [react()],
  server: {
    // In dev, forward /api and /health to a locally-running api/app.py
    // (uvicorn api/app.py, default port 8080) so VITE_API_BASE_URL can stay
    // unset locally; set it explicitly for a deployed API in production.
    proxy: {
      '/api': 'http://localhost:8080',
      '/health': 'http://localhost:8080',
    },
  },
})
