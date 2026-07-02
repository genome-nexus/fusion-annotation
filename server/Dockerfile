FROM python:3.12-slim

WORKDIR /app

# Install the package (core + server extras) first for layer caching.
COPY pyproject.toml README.md LICENSE ./
COPY src ./src
RUN pip install --no-cache-dir ".[server]"

COPY server ./server

ENV PORT=8080
EXPOSE 8080

# Cloud Run sets $PORT; uvicorn reads it via server/app.py's __main__ block.
CMD ["python", "server/app.py"]
