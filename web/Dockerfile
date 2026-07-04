# Multi-stage build: compile the Vite/React SPA, then serve the static
# output with a minimal nginx image. VITE_API_BASE_URL must be baked in at
# *build* time (Vite inlines env vars into the JS bundle), so it's a build
# arg here, not a runtime env var like the api/ and server/ containers.
FROM node:20-slim AS build
WORKDIR /app
COPY web/package.json web/package-lock.json ./
RUN npm ci
COPY web/ ./
ARG VITE_API_BASE_URL=""
ENV VITE_API_BASE_URL=${VITE_API_BASE_URL}
RUN npm run build

FROM nginx:1.27-alpine
COPY --from=build /app/dist /usr/share/nginx/html
COPY web/nginx.conf /etc/nginx/conf.d/default.conf
EXPOSE 8080
