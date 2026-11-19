# Legacy trees (do not expand until Phase 2)

These directories remain in the repository for reference and historical behaviour.
They are **out of scope** for the Phase 0/1 tooling and unit-test work.

| Path | Role | Plan |
|------|------|------|
| `src/oo_refactoring/` | Partial object-oriented Streamlit rewrite | Inventory in Phase 2; merge useful pieces or delete |
| `src/bulk/` | Bulk OZ / DFT experiment tree | Consolidate with `src/pyoz` in Phase 3 |
| `src/pyoz/` | Vendored PyOZ bulk solver | Single bulk path decision in Phase 3 |
| `src/tests/` | Old nested tests | Superseded by top-level `tests/` |
| `data/fluid_parameters.toml` | Parallel parameter table | Single source of truth in Phase 2 (`charged_fluids`) |
| `requirements*.txt` | Pre-uv freezes | Prefer `uv.lock`; export only if deploy needs a freeze |

Active science code lives in **`src/charged_fluids/`**. Streamlit entrypoints under
`src/Main.py` and `src/pages/` still import thin shims (`src/modelling.py`, etc.).
