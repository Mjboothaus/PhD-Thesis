# UV + Python 3.14 slim image for Render / Cloud Run / local containers
FROM python:3.14-slim-bookworm

RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /usr/local/bin/

WORKDIR /app

COPY pyproject.toml uv.lock README.md ./
COPY src ./src
COPY data ./data
COPY docs ./docs

RUN mkdir -p .streamlit
COPY .streamlit/ ./.streamlit/

ENV UV_COMPILE_BYTECODE=1 \
    UV_LINK_MODE=copy \
    UV_NO_DEV=1 \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONPATH=/app/src \
    PATH=/app/.venv/bin:$PATH \
    PORT=8080 \
    STREAMLIT_SERVER_HEADLESS=true \
    STREAMLIT_BROWSER_GATHER_USAGE_STATS=false

# Install project + app group only (no default-groups / no dev tooling)
RUN uv sync --frozen --no-cache --no-dev --group app \
    && test -x /app/.venv/bin/streamlit

EXPOSE 8080

HEALTHCHECK --interval=30s --timeout=10s --start-period=40s --retries=3 \
    CMD curl -fsS "http://127.0.0.1:${PORT}/_stcore/health" || exit 1

# Use the image venv directly so runtime does not re-resolve default-groups (dev)
CMD ["sh", "-c", "exec streamlit run src/Main.py --server.port=${PORT:-8080} --server.address=0.0.0.0 --server.headless=true"]
