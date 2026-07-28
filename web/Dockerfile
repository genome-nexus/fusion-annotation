FROM node:20-alpine AS builder
WORKDIR /app
# Install deps first for layer caching
COPY web/package*.json ./
RUN npm ci
COPY web/ ./
# Build with base=/ so the app is served at the root path.
# (The vite.config.ts default uses /fusion-annotation/ for GitHub Pages;
# Docker runs at root so we override it here.)
RUN npm run build -- --base=/

FROM nginx:alpine
COPY --from=builder /app/dist /usr/share/nginx/html
COPY web/nginx.conf /etc/nginx/conf.d/default.conf
EXPOSE 80
