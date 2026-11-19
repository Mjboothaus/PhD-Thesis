# Legacy trees

These directories remain for reference and historical behaviour. **Do not expand them.**
Active science code lives in **`src/charged_fluids/`**.

| Path | Role | Status |
|------|------|--------|
| `src/oo_refactoring/` | Partial object-oriented Streamlit rewrite | Legacy; use main Streamlit + `charged_fluids` |
| `src/bulk/` | Early bulk OZ / DFT experiment tree | Superseded by `charged_fluids.bulk` |
| `src/pyoz/` | Vendored PyOZ bulk solver | Reference only (`just bulk-fluid-pyoz`) |
| `src/tests/` | Old nested tests | Superseded by top-level `tests/` |
| `notebooks/legacy/` | Pre-Marimo Jupyter notebooks | Archival; use `just marimo` / `just marimo-bulk` |
| `data/fluid_parameters.toml` | Fluid parameter table | **Current** single source of truth (loaded by `charged_fluids.parameters`) |
| `requirements*.txt` (removed) | Pre-uv freezes | Use `uv.lock`; `just export-reqs` writes a freeze on demand |

Streamlit entrypoints under `src/Main.py` and `src/pages/` import thin shims that delegate to `charged_fluids`.

## Deployment

- **Current production:** Render service `PhD-Thesis-MJB` — https://phd-thesis-mjb.onrender.com
- **Deprecated:** Google Cloud Run (and any former custom-domain path that depended on it). Do not target Cloud Run for new deploys.

## Notebooks

Active Marimo explorers (Plotly figures):

- `notebooks/wall_fluid_explorer.py` — wall singlet + bulk inputs (`just marimo`)
- `notebooks/bulk_fluid_explorer.py` — bulk artefacts / parameters / OZ (`just marimo-bulk`)

See `notebooks/legacy/README.md` for the mapping from old `.ipynb` files.
