# Phase 2 — Core package depth and thinner Streamlit pages

**Date:** 2026-11-19

## What changed

- **`data/fluid_parameters.toml`** is the single source of truth (issue #16),
  including potential blocks for LJ and ionic KCl.
- **`charged_fluids.setup.prepare_fluid_run`** builds fluid + grid + bulk c(r)
  integrals + pair potential + `Model` in one place.
- **`Model`** now always carries `hw`, `tw`, `c_short`, `f1`, `f2` (issue #21).
- LJ / KCl Streamlit pages are thin wrappers (no star-imports, shared setup).
- LJ bulk tab plots **`βu(r)`** so the potential is visible (issue #3).
- Solver UI reports **SciPy `nit`** explicitly (issue #14).

## Issue triage

See `docs/ISSUES_TRIAGE.md` — most open GitHub issues are four years old; several
are done or superseded; KCl bulk c(r) (#13/#17) remains the scientific blocker.

## Still deferred

- `src/oo_refactoring/` and dual bulk trees (documented legacy)
- On-the-fly PyOZ bulk generation (Phase 3)
