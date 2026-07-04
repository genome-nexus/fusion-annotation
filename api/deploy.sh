#!/usr/bin/env bash
# Deploy the fusion-annotation public REST API to Cloud Run.
#
# Usage:
#   export GCP_PROJECT=your-project-id
#   ./api/deploy.sh
#
# Requires: gcloud CLI, authenticated (`gcloud auth login` or a service
# account with Cloud Run Admin + Artifact Registry Writer + Cloud Build
# Editor on $GCP_PROJECT), and the Cloud Run / Cloud Build / Artifact
# Registry APIs enabled on the project (the script enables them if needed).
#
# Cost: Cloud Run's free tier covers 2M requests/month; --max-instances=2
# below caps a runaway bill, and FUSION_ANNOTATION_RATE_LIMIT (default
# 30/minute per IP) protects against a single abusive caller.
set -euo pipefail

PROJECT="${GCP_PROJECT:?Set GCP_PROJECT to your GCP project id}"
REGION="${GCP_REGION:-us-central1}"
SERVICE="${SERVICE_NAME:-fusion-annotation-api}"
REPO="${AR_REPO:-fusion-annotation}"
# Comma-separated list of allowed CORS origins (e.g. the deployed web/ SPA's
# origin). Defaults to "*" (open) for a public demo API.
CORS_ORIGINS="${FUSION_ANNOTATION_CORS_ORIGINS:-*}"

cd "$(dirname "$0")/.."   # repo root

echo "== Project: $PROJECT   Region: $REGION   Service: $SERVICE =="

gcloud config set project "$PROJECT" >/dev/null
gcloud services enable run.googleapis.com cloudbuild.googleapis.com \
    artifactregistry.googleapis.com --project "$PROJECT"

# Artifact Registry repo for the container image (idempotent; shared with the
# MCP server's images under the same repo, different image name).
gcloud artifacts repositories describe "$REPO" --location="$REGION" >/dev/null 2>&1 || \
  gcloud artifacts repositories create "$REPO" --repository-format=docker \
    --location="$REGION" --description="fusion-annotation container images"

IMAGE="$REGION-docker.pkg.dev/$PROJECT/$REPO/$SERVICE:latest"

echo "== Building image via Cloud Build: $IMAGE =="
gcloud builds submit --config=api/cloudbuild.yaml --substitutions="_IMAGE=$IMAGE" .

echo "== Deploying to Cloud Run =="
gcloud run deploy "$SERVICE" \
  --image "$IMAGE" \
  --region "$REGION" \
  --allow-unauthenticated \
  --max-instances=2 \
  --min-instances=0 \
  --memory=512Mi \
  --cpu=1 \
  --update-env-vars="^;^FUSION_ANNOTATION_CORS_ORIGINS=$CORS_ORIGINS" \
  --quiet

URL=$(gcloud run services describe "$SERVICE" --region "$REGION" --format='value(status.url)')

echo ""
echo "=================================================================="
echo "Deployed:  $URL"
echo "Annotate endpoint: $URL/api/annotate?five_gene=EML4&three_gene=ALK&five_exon=13&three_exon=20"
echo "API docs:          $URL/api/docs"
echo "Health check:      $URL/health"
echo ""
echo "Point the web/ SPA's VITE_API_BASE_URL at: $URL"
echo "Once the SPA is deployed, re-run with FUSION_ANNOTATION_CORS_ORIGINS"
echo "set to its origin to lock down CORS (defaults to * otherwise)."
echo "=================================================================="
