# UV + Python 3.14 slim image for Cloud Run / local container runs
FROM python:3.14-slim-bookworm

RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    && rm -rf /var/lib/apt/lists/*

COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /usr/local/bin/

WORKDIR /app

COPY pyproject.toml uv.lock README.md ./
COPY src ./src
COPY data ./data
COPY docs ./docs

# Runtime: project + app UI group (no dev tooling)
RUN uv sync --frozen --no-cache --no-dev --group app

ENV PYTHONPATH=/app/src:$PYTHONPATH

EXPOSE 8080

HEALTHCHECK --interval=30s --timeout=30s --start-period=5s --retries=3 \
    CMD curl -f http://localhost:8080/_stcore/health || exit 1

CMD ["uv", "run", "streamlit", "run", "src/Main.py", "--server.port=8080", "--server.address=0.0.0.0"]
