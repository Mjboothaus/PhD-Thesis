# Phase 1b — residual coverage without a Fortran oracle

**Date:** 2026-11-19  
**Phase:** 1b — `calc_tw` / residual characterisation

## Context

After Phase 1a the `charged_fluids` package already sat near 97% line coverage.
The remaining risk was the heart of the singlet equation — `calc_tw` and
`opt_func` — which the Streamlit UI only exercises via a full Newton–Krylov
solve.

## What we added

- Synthetic one- and two-component neutral setups (zero charge_pair, mild wall)
- `calc_tw` finiteness and sensitivity to non-zero f1/f2 kernels
- `opt_func` identity: residual ≡ tw_in − calc_tw(...)
- Coverage fillers: 1-D c(r) interpolate path, `get_memory_usage`, LJ2 params

## Why this matters

Without the original Fortran, the safest move is to lock *algebraic* properties
of the residual before any refactor of loops or pair indexing. Full physical
validation still waits on bulk c(r) quality (KCl xfail) and digitised figures.

## Commands

```bash
uv run pytest --cov=charged_fluids --cov-report=term-missing
```
