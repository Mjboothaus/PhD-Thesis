# Phase 3 — Bulk c(r) pipeline

**Date:** 2026-11-19

## Delivered

- New package module `charged_fluids.bulk`:
  - radial sine transforms
  - Picard HNC/PY bulk OZ (1-comp LJ solid; multi-comp charged experimental)
  - thesis long-range split helper: `c_short = c + beta u_coulomb` (wall-fluid units)
  - artefacts under `data/bulk/{symbol}.cr.dat` + `.meta.json`
- CLI: `uv run python -m charged_fluids.bulk.cli lj1`
- Just: `just bulk-cr lj1` / `just bulk-cr-all`
- `prepare_fluid_run` prefers `data/bulk/{symbol}.cr.dat` over legacy `data/*.data`

## LJ

Generated `data/bulk/lj1.cr.dat` converges quickly; short-range tail ~1e-5.

## KCl / issues #13 #17

Full melt-density multi-component OZ remains **experimental** (xfail). The
pipeline and long-range *definition* are in place; numerical stability at
molten-salt state points needs further Ng/Fourier work (legacy pyoz remains
reference). Legacy `nrcg-cr-approx-maybe.dat` tail |c|~0.01 — still the
production fallback until charged bulk generation is solid.

## Next

- Stabilise charged bulk OZ (Ng during iteration, better mixing)
- Regenerate KCl artefact and flip xfail
- Optional: wrap vendored `src/pyoz` as alternate backend
