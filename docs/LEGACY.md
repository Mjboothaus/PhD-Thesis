# Legacy trees

These directories remain for reference and historical behaviour. **Do not expand them.**
Active science code lives in **`src/charged_fluids/`**.

| Path | Role | Status |
|------|------|--------|
| `src/oo_refactoring/` | Partial object-oriented Streamlit rewrite | Legacy; use main Streamlit + `charged_fluids` |
| `src/bulk/` | Early bulk OZ / DFT experiment tree | Superseded by `charged_fluids.bulk` |
| `src/pyoz/` | Vendored PyOZ bulk solver | Reference only (`just bulk-fluid-pyoz`) |
| `src/tests/` | Old nested tests | Superseded by top-level `tests/` |
| `data/fluid_parameters.toml` | Fluid parameter table | **Current** single source of truth (loaded by `charged_fluids.parameters`) |
| `requirements*.txt` (removed) | Pre-uv freezes | Use `uv.lock`; `just export-reqs` writes a freeze on demand |

Streamlit entrypoints under `src/Main.py` and `src/pages/` import thin shims that delegate to `charged_fluids`.
