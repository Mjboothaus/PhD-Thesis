# Rebuilding confidence without the Fortran

**Date:** 2026-11-19  
**Phase:** 1a — Test harness and pure-function tests

## Context

The original PhD numerical code was Fortran 77. Those sources appear lost. The
Python port had only three unit tests nested under `src/tests/`, while the
Streamlit UI exercised Lennard-Jones somewhat successfully and left charged
fluids on shaky bulk inputs.

Without a bit-identical gold master, confidence has to come from:

1. Physical invariants (signs, limits, dimensions)
2. Analytic checks on simple grids
3. Characterisation of *current* behaviour (including known wrongness)
4. Digitised thesis figures later

## What we locked in

- Top-level `tests/` importing **`charged_fluids` only** (no Streamlit).
- Helpers: beta, epsilon, rho, charge, kappa (KCl ≈ 6.77 Å⁻¹), pair index.
- LJ potential minimum near \(r = 2^{1/6}\sigma\); wall potential decay.
- Load/interpolate `data/lj1-cr.data`; integral transforms on synthetic c(r).
- Optional `@pytest.mark.slow` LJ Newton–Krylov smoke solve on a small grid.
- `@pytest.mark.regression` for KCl bulk c(r) tail not vanishing — documents
  the charged-fluid gap called out in the README.

## Design notes

- Core modelling no longer imports Streamlit `@cache_data`.
- `set_fluid_parameters` deep-copies so tests cannot corrupt a global table.
- Integrals now assign `f1[k]` via `scipy.integrate.trapezoid` (clearer than
  broadcasting a scalar onto `f1[k:]`, behaviourally equivalent end state).

## Next

- Raise coverage on `calc_tw` / residual paths.
- Phase 2: deepen package boundaries and collapse legacy trees.
- Phase 3: bulk c(r) generation so the KCl regression can flip to a pass.
