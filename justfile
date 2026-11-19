# Docs: https://just.systems/man/en/
# Recipes call `uv` only so they run on macOS, Windows, and Linux without shell forks.
# Multipage Streamlit imports resolve via PYTHONPATH=src (also set in Dockerfile).

set windows-shell := ["powershell.exe", "-NoLogo", "-Command"]

# Ensure src/ is importable for components/sidebar/plotting in local runs
export PYTHONPATH := "src"

# Default recipe when just is called without arguments
default:
    @just --list

# Show available recipes
help:
    @just --list

# --- Development environment (uv only) ---

# Synchronise environment with pyproject.toml / uv.lock (default groups)
sync:
    uv sync

# Initialise development environment (alias)
init: sync

# App UI dependencies (Streamlit, Plotly, …)
sync-app:
    uv sync --group app

# Notebook group (Marimo + Plotly)
sync-notebook:
    uv sync --group notebook

# Refresh the lockfile from pyproject.toml
lock:
    uv lock

# Export a frozen requirements file (deploy escape hatch only)
export-reqs:
    uv export --frozen --no-dev --group app -o requirements.prod.txt

# --- Application ---

# Run the Streamlit app (localhost; uses PYTHONPATH=src)
app app_name="src/Main.py": sync-app
    uv run streamlit run {{app_name}} --server.address=localhost

# Run the experimental OO Streamlit entry (legacy; see docs/LEGACY.md)
app-opt app_name="src/oo_refactoring/Main.py": sync-app
    uv run streamlit run {{app_name}} --server.address=localhost

# Update Streamlit config snapshot
update-st-config:
    uv run streamlit config show > .streamlit/config.toml

# --- Bulk c(r) generation ---

# Generate bulk short-range c(r) into data/bulk/{symbol}.cr.dat
# Grid/closure defaults come from data/fluid_parameters.toml [symbol.bulk].
bulk-cr symbol="lj1":
    uv run python -m charged_fluids.bulk.cli {{symbol}}

# Dry-run bulk solve (no artefact write)
bulk-cr-dry symbol="lj1":
    uv run python -m charged_fluids.bulk.cli {{symbol}} --dry-run

# Optional explicit grid override
bulk-cr-grid symbol="lj1" n_points="1024" delta_r="0.05":
    uv run python -m charged_fluids.bulk.cli {{symbol}} --n-points {{n_points}} --delta-r {{delta_r}}

# List fluids that support bulk generation
bulk-cr-list:
    uv run python -m charged_fluids.bulk.cli --list

# Generate all supported bulk artefacts (lj, ionic, cf1 digitised)
bulk-cr-all:
    uv run python -m charged_fluids.bulk.cli --all

# Build CF1 water bulk c(r) from digitised thesis Fig 3.12 (alias of bulk-cr h2o)
bulk-cr-h2o:
    uv run python -m charged_fluids.bulk.cli h2o

# Legacy pyOZ path (vendored; see docs/LEGACY.md)
bulk-fluid-pyoz input_file="lj/nrcg-lj-1comp.in":
    uv run python src/pyoz/pyoz.py -i src/pyoz/tests/{{input_file}}

# --- Docker / deploy ---

docker-build project_name="phd_thesis":
    docker build . -t {{project_name}}

docker-run project_name="phd_thesis" server_port="8080":
    docker run --rm -e PORT={{server_port}} -p {{server_port}}:{{server_port}} {{project_name}}

container: docker-build docker-run

# Deploy via Render Blueprint (requires `render login` and linked repo)
# Uses render.yaml; Docker image is built on Render from Dockerfile.
deploy-render:
    render blueprints launch --confirm

# Trigger a new deploy of an existing Render service by name (if CLI is authed)
deploy-render-service service_name="phd-thesis":
    render deploys create {{service_name}} --confirm

# --- Testing and quality (ruff, ty, prek, pytest) ---

test:
    uv run pytest

test-cov:
    uv run pytest --cov=charged_fluids --cov-report=term-missing --cov-report=html

typecheck:
    uv run ty check

lint:
    uv run ruff check src/charged_fluids tests

fmt-check:
    uv run ruff format --check src/charged_fluids tests

format:
    uv run ruff check --fix src/charged_fluids tests
    uv run ruff format src/charged_fluids tests

hooks-install:
    uv run prek install

hooks:
    uv run prek run --all-files

quality: format typecheck lint test

# --- Marimo notebooks (Plotly figures) ---

marimo: sync-notebook
    uv run marimo edit notebooks/wall_fluid_explorer.py

marimo-bulk: sync-notebook
    uv run marimo edit notebooks/bulk_fluid_explorer.py

# --- Blog notes ---

blog-list:
    @uv run python -c "from pathlib import Path; print('\\n'.join(sorted(p.name for p in Path('docs/blog').glob('*'))))"

# --- Clean (portable via Python; avoids shell glob differences) ---

clean:
    uv run python -c "import pathlib, shutil; roots=['.pytest_cache','.mypy_cache','.ruff_cache','.ty_cache','htmlcov']; [shutil.rmtree(p, ignore_errors=True) for p in roots]; p=pathlib.Path('.coverage'); p.unlink(missing_ok=True); [shutil.rmtree(x, ignore_errors=True) for x in pathlib.Path('.').rglob('__pycache__')]; print('clean done')"
