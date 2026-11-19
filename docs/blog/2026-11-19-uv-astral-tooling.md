# From mixed requirements to uv-only + Astral quality gates

**Date:** 2026-11-19  
**Phase:** 0 — Foundations

## Context

The thesis reimplementation already had `pyproject.toml` and `uv.lock`, but the
repo still carried parallel `requirements*.txt` files and a black / isort /
flake8 / mypy quality stack. Python was declared as `>=3.13` while day-to-day
work ran on 3.14.

## Decisions

- Pin **Python 3.14** (`.python-version`, `requires-python = ">=3.14,<3.15"`).
- Use **uv dependency groups** (`dev`, `notebook`, `docs`) instead of optional
  extras for tooling.
- Replace black/isort/flake8/mypy with **ruff** (lint + format) and **ty**
  (types), driven by **prek** hooks.
- Keep **Just** as the cross-platform task runner; recipes call `uv run …` only.
- Introduce package **`charged_fluids`** under `src/` as the science core;
  Streamlit keeps thin import shims for now.
- Leave `oo_refactoring/`, dual bulk trees, and old nested tests documented in
  `docs/LEGACY.md` until Phase 2.

## Commands that matter

```bash
uv sync
just test
just format
just typecheck
just hooks-install
just quality
```

## Lock upgrade

As part of Phase 0 we ran `uv lock --upgrade` so the lockfile tracks current
latest wheels within the declared bounds (notably NumPy 2.5.x, SciPy 1.18.x,
pandas 3.x, Streamlit 1.67.x, Ruff 0.17.x, ty 0.0.107). Core science stays
independent of Streamlit via the `app` dependency group; `pyarrow>=26` avoids
source builds on Python 3.14.

## Gotchas

- `uv run pytest` previously failed while building `pyarrow` from source when
  Arrow CMake files were missing — prefer binary wheels and a clean `uv sync`
  on 3.14.
- prek understands `.pre-commit-config.yaml`; pin ruff/ty hook revs and refresh
  with `uv run prek update` when ready.
- Docker base image moved to `python:3.14-slim-bookworm` to match the pin.

## Open questions

- Whether Cloud Run still needs an exported `requirements.prod.txt` or can rely
  solely on `uv sync --frozen` in the image build.
- When to pull Marimo into the default dev group versus keeping it optional.
