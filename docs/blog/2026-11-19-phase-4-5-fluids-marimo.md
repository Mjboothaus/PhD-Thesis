# Phase 4–5 — LJ2 fluid and Marimo explorer

**Date:** 2026-11-19

## Fluids (Phase 4)

- **LJ2** (equal two-component) marked `ready` with bulk artefact
  `data/bulk/lj2.cr.dat` (pair c(r) from 1-comp OZ at total density).
- Shared Streamlit helper `components/fluid_page.py` drives LJ1, LJ2, KCl.
- Water / 2-2 remain stubs (`not_implemented`).

## Marimo (Phase 5)

- `notebooks/wall_fluid_explorer.py` — same `charged_fluids` core as Streamlit.
- `just marimo` → `uv sync --group notebook` + `marimo edit`.

## Fluid status matrix

| Symbol | Status | Bulk c(r) |
|--------|--------|-----------|
| lj1 | ready | `data/bulk/lj1.cr.dat` |
| lj2 | ready | `data/bulk/lj2.cr.dat` |
| kcl | experimental | legacy fallback |
| h2o / 2_2 | not_implemented | — |
