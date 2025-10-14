# Modern Dockerfile using UV and Python 3.13
FROM python:3.13-slim-bookworm

# Install system dependencies
RUN apt-get update && apt-get install -y \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Install UV
COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /usr/local/bin/

# Create app directory
WORKDIR /app

# Copy dependency files
COPY pyproject.toml uv.lock ./

# Install dependencies using UV
RUN uv sync --frozen --no-cache --no-dev

# Copy application code
COPY docs ./docs
COPY src ./src
COPY data ./data

# Expose port
EXPOSE 8080

# Run the application using UV
CMD ["uv", "run", "streamlit", "run", "src/Main.py", "--server.port=8080", "--server.address=0.0.0.0"]
