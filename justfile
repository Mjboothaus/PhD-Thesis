# Docs: https://just.systems/man/en/

# Default recipe when just is called without arguments
default:
    @just --list

# Show available recipes
help:
    @just --list

# Development Environment Setup

# Initialize development environment with UV
init:
    uv pip install -r requirements.txt
    uv pip install -r requirements-dev.txt

# Update UV environment
update:
    uv pip freeze > requirements.txt
    uv pip install -r requirements.txt --upgrade

# Run the Streamlit app
app app_name="src/Main.py":
    uv run streamlit run {{app_name}} --server.address=localhost

# Run the app with optimized output visualizations
app-opt app_name="src/oo_refactoring/Main.py":
    uv run streamlit run {{app_name}} --server.address=localhost

# Update Streamlit config
update-st-config:
    uv run streamlit config show > .streamlit/config.toml

# Generate requirements files
reqs:
    uv pip freeze > requirements.txt
    uv pip freeze > requirements-dev.txt

# Install requirements for development
install-dev:
    uv pip install -r requirements-dev.txt


# Run bulk-fluid (pyOz) - default 1 component LJ
bulk-fluid-pyoz input_file="lj/nrcg-lj-1comp.in":
    cd src/pyoz && uv run python pyoz.py -i tests/{{input_file}}

# Docker and Deployment

# Build Docker image
docker-build project_name="phd_thesis":
    docker build . -t {{project_name}}

# Run Docker container
docker-run project_name="phd_thesis" server_port="8080":
    docker run -p {{server_port}}:{{server_port}} {{project_name}}

# Build and run Docker container
container: docker-build docker-run

# Deploy to Render.com
deploy-render:
    git push render main

# Testing and Quality

# Run all tests
test:
    uv run pytest

# Run tests with coverage
test-cov:
    uv run pytest --cov=src --cov-report=html

# Run type checking
type-check:
    uv run mypy src

# Run linting
lint:
    uv run black src
    uv run flake8 src

# Format code
format:
    uv run black src
    uv run isort src

# Run all quality checks
quality: format type-check lint test

# Clean up
clean:
    rm -rf .pytest_cache
    rm -rf .mypy_cache
    rm -rf .coverage
    rm -rf htmlcov
    rm -rf **/__pycache__
    rm -rf **/*.pyc
