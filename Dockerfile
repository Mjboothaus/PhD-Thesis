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

# Streamlit config is optional; create empty dir so COPY always succeeds
RUN mkdir -p .streamlit
COPY .streamlit/ ./.streamlit/

ENV UV_COMPILE_BYTECODE=1 \
    UV_LINK_MODE=copy \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONPATH=/app/src \
    PORT=8080 \
    STREAMLIT_SERVER_HEADLESS=true \
    STREAMLIT_BROWSER_GATHER_USAGE_STATS=false

RUN uv sync --frozen --no-cache --no-dev --group app

EXPOSE 8080

HEALTHCHECK --interval=30s --timeout=10s --start-period=40s --retries=3 \
    CMD curl -fsS "http://127.0.0.1:${PORT}/_stcore/health" || exit 1

# Respect Render/Cloud Run $PORT; bind all interfaces; headless
CMD ["sh", "-c", "exec uv run streamlit run src/Main.py --server.port=${PORT:-8080} --server.address=0.0.0.0 --server.headless=true"]
