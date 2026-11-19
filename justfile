# Docs: https://just.systems/man/en/
# Recipes are intended to run unchanged on macOS, Windows, and Linux via uv.

set windows-shell := ["powershell.exe", "-NoLogo", "-Command"]

# Default recipe when just is called without arguments
default:
    @just --list

# Show available recipes
help:
    @just --list

# --- Development environment (uv only) ---

# Synchronise environment with pyproject.toml / uv.lock
sync:
    uv sync

# Initialise development environment (alias)
init: sync

# Refresh the lockfile from pyproject.toml
lock:
    uv lock

# Export a frozen requirements file (deploy escape hatch only)
export-reqs:
    uv export --frozen --no-dev --group app -o requirements.prod.txt

# --- Application ---

# Run the Streamlit app
app app_name="src/Main.py":
    uv run streamlit run {{app_name}} --server.address=localhost

# Run the experimental OO Streamlit entry (legacy; see docs/LEGACY.md)
app-opt app_name="src/oo_refactoring/Main.py":
    uv run streamlit run {{app_name}} --server.address=localhost

# Update Streamlit config snapshot
update-st-config:
    uv run streamlit config show > .streamlit/config.toml

# --- Bulk fluid (legacy pyOZ path) ---

# Run bulk-fluid (pyOz) - default 1 component LJ
bulk-fluid-pyoz input_file="lj/nrcg-lj-1comp.in":
    cd src/pyoz && uv run python pyoz.py -i tests/{{input_file}}

# --- Docker ---

docker-build project_name="phd_thesis":
    docker build . -t {{project_name}}

docker-run project_name="phd_thesis" server_port="8080":
    docker run -p {{server_port}}:{{server_port}} {{project_name}}

container: docker-build docker-run

# --- Testing and quality (ruff, ty, prek, pytest) ---

# Run unit tests
test:
    uv run pytest

# Run tests with coverage for charged_fluids
test-cov:
    uv run pytest --cov=charged_fluids --cov-report=term-missing --cov-report=html

# Type-check core package and tests
typecheck:
    uv run ty check

# Lint (check only)
lint:
    uv run ruff check src/charged_fluids tests

# Format check (no writes)
fmt-check:
    uv run ruff format --check src/charged_fluids tests

# Apply ruff formatter + autofixes
format:
    uv run ruff check --fix src/charged_fluids tests
    uv run ruff format src/charged_fluids tests

# Install git hooks via prek
hooks-install:
    uv run prek install

# Run all prek hooks on the tree
hooks:
    uv run prek run --all-files

# Format, typecheck, lint, test
quality: format typecheck lint test

# --- Blog notes ---

# List progressive blog notes
blog-list:
    @echo "Blog notes under docs/blog/"
    @ls docs/blog

# --- Clean ---

clean:
    rm -rf .pytest_cache
    rm -rf .mypy_cache
    rm -rf .ruff_cache
    rm -rf .coverage
    rm -rf htmlcov
    rm -rf .ty_cache
    rm -rf **/__pycache__
    rm -rf **/*.pyc
