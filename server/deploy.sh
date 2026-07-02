#!/usr/bin/env bash
# Deploy the fusion-annotation MCP server to Cloud Run.
#
# Usage:
#   export GCP_PROJECT=your-project-id
#   ./server/deploy.sh
#
# Requires: gcloud CLI, authenticated (`gcloud auth login` or a service
# account with Cloud Run Admin + Artifact Registry Writer + Cloud Build
# Editor on $GCP_PROJECT), and the Cloud Run / Cloud Build / Artifact
# Registry APIs enabled on the project (the script enables them if needed).
#
# Cost: Cloud Run's free tier covers 2M requests/month; --max-instances=2
# below caps a runaway bill. A handful of reviewers hitting this
# occasionally should cost $0.
set -euo pipefail

PROJECT="${GCP_PROJECT:?Set GCP_PROJECT to your GCP project id}"
REGION="${GCP_REGION:-us-central1}"
SERVICE="${SERVICE_NAME:-fusion-annotation-mcp}"
REPO="${AR_REPO:-fusion-annotation}"

cd "$(dirname "$0")/.."   # repo root

echo "== Project: $PROJECT   Region: $REGION   Service: $SERVICE =="

gcloud config set project "$PROJECT" >/dev/null
gcloud services enable run.googleapis.com cloudbuild.googleapis.com \
    artifactregistry.googleapis.com --project "$PROJECT"

# Artifact Registry repo for the container image (idempotent).
gcloud artifacts repositories describe "$REPO" --location="$REGION" >/dev/null 2>&1 || \
  gcloud artifacts repositories create "$REPO" --repository-format=docker \
    --location="$REGION" --description="fusion-annotation MCP server images"

IMAGE="$REGION-docker.pkg.dev/$PROJECT/$REPO/$SERVICE:latest"

echo "== Building image via Cloud Build: $IMAGE =="
gcloud builds submit --tag "$IMAGE" .

echo "== First deploy pass (to learn the assigned *.run.app hostname) =="
gcloud run deploy "$SERVICE" \
  --image "$IMAGE" \
  --region "$REGION" \
  --allow-unauthenticated \
  --max-instances=2 \
  --min-instances=0 \
  --memory=512Mi \
  --cpu=1 \
  --quiet

URL=$(gcloud run services describe "$SERVICE" --region "$REGION" --format='value(status.url)')
HOST=$(python3 -c "from urllib.parse import urlparse; print(urlparse('$URL').hostname)")

echo "== Locking transport security to host: $HOST =="
gcloud run services update "$SERVICE" \
  --region "$REGION" \
  --update-env-vars="FUSION_ANNOTATION_ALLOWED_HOSTS=$HOST" \
  --quiet

echo ""
echo "=================================================================="
echo "Deployed:  $URL"
echo "MCP endpoint:      $URL/mcp"
echo "Health check:      $URL/healthz"
echo ""
echo "In Claude.ai / Claude Desktop: Settings -> Connectors -> Add custom"
echo "connector -> URL: $URL/mcp"
echo "=================================================================="
