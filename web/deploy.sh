#!/usr/bin/env bash
# Deploy the fusion-annotation web UI (React SPA) to Cloud Run.
#
# Usage:
#   export GCP_PROJECT=your-project-id
#   export API_URL=https://fusion-annotation-api-xxxxx.a.run.app   # from api/deploy.sh
#   ./web/deploy.sh
#
# Requires: gcloud CLI, authenticated, same permissions as api/deploy.sh and
# server/deploy.sh (Cloud Run Admin + Artifact Registry Writer + Cloud Build
# Editor on $GCP_PROJECT).
#
# Note: VITE_API_BASE_URL is baked into the JS bundle at *build* time (Vite
# inlines env vars), so changing the API URL later requires rebuilding and
# redeploying this image — it can't be changed via `gcloud run services
# update --update-env-vars` like the api/ and server/ containers.
set -euo pipefail

PROJECT="${GCP_PROJECT:?Set GCP_PROJECT to your GCP project id}"
REGION="${GCP_REGION:-us-central1}"
SERVICE="${SERVICE_NAME:-fusion-annotation-web}"
REPO="${AR_REPO:-fusion-annotation}"
API_URL="${API_URL:?Set API_URL to the deployed api/ service URL (see api/deploy.sh output)}"

cd "$(dirname "$0")/.."   # repo root

echo "== Project: $PROJECT   Region: $REGION   Service: $SERVICE =="
echo "== Building against API: $API_URL =="

gcloud config set project "$PROJECT" >/dev/null
gcloud services enable run.googleapis.com cloudbuild.googleapis.com \
    artifactregistry.googleapis.com --project "$PROJECT"

gcloud artifacts repositories describe "$REPO" --location="$REGION" >/dev/null 2>&1 || \
  gcloud artifacts repositories create "$REPO" --repository-format=docker \
    --location="$REGION" --description="fusion-annotation container images"

IMAGE="$REGION-docker.pkg.dev/$PROJECT/$REPO/$SERVICE:latest"

echo "== Building image via Cloud Build: $IMAGE =="
gcloud builds submit --config=web/cloudbuild.yaml \
  --substitutions="_IMAGE=$IMAGE,_VITE_API_BASE_URL=$API_URL" .

echo "== Deploying to Cloud Run =="
gcloud run deploy "$SERVICE" \
  --image "$IMAGE" \
  --region "$REGION" \
  --allow-unauthenticated \
  --max-instances=2 \
  --min-instances=0 \
  --memory=256Mi \
  --cpu=1 \
  --quiet

URL=$(gcloud run services describe "$SERVICE" --region "$REGION" --format='value(status.url)')

echo ""
echo "=================================================================="
echo "Deployed:  $URL"
echo ""
echo "Re-run api/deploy.sh with FUSION_ANNOTATION_CORS_ORIGINS=$URL to lock"
echo "down the API's CORS policy to this SPA's origin."
echo "=================================================================="
