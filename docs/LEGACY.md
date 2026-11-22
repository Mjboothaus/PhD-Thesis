# Legacy and redundant inventory

**Do not expand legacy trees.** Active science lives in **`src/charged_fluids/`**.  
Active UI: Streamlit under `src/Main.py` + `src/pages/` + `src/components/`.  
Active explorers: `notebooks/wall_fluid_explorer.py`, `notebooks/bulk_fluid_explorer.py`.  
**Production:** Render — https://phd-thesis-mjb.onrender.com (Cloud Run deprecated).

Reviewed: 2026-11-19 (post PRs #22–#36).  
Safe deletes applied: 2026-11-20 (`src/phd_thesis/`, `src/tests/`, `src/non_app_plotting.py`, `data/pyOZ_bulk_fluid/`; `logs/` gitignored).

---

## Already documented legacy (keep for now, do not grow)

| Path | Role | Recommendation |
|------|------|----------------|
| `src/oo_refactoring/` | Half-finished OO Streamlit rewrite | **Removed** — Krylov residual history lifted into `charged_fluids.ConvergenceHistory` / `solve_model` |
| `src/bulk/` | Early bulk OZ / DFT experiment | **Removed** — use `charged_fluids.bulk` |
| `src/pyoz/` (~15M) | Vendored PyOZ | **Reference only** (`just bulk-fluid-pyoz`) — keep until bulk science is trusted, then archive off-repo |
| `src/tests/` | Old nested tests | **Removed** 2026-11-20 — use top-level `tests/` |
| `notebooks/legacy/` | Pre-Marimo Jupyter | **Archival** — keep; see `notebooks/legacy/README.md` |

---

## Compatibility shims (still used by Streamlit edge)

| Path | Role | Recommendation |
|------|------|----------------|
| `src/modelling.py` | Re-exports `charged_fluids.modelling` | Keep until pages/components import `charged_fluids` only |
| `src/parameters.py` | Re-exports `charged_fluids.parameters` | Same |
| `src/numerics.py` | Re-exports `charged_fluids.numerics` | Same |
| `src/plotting.py` | Streamlit Plotly helpers | **Active** — not legacy |
| `src/sidebar.py` | Streamlit sidebar | **Active** |
| `src/helper_functions.py` | Markdown render helper | **Active** (Main + Theory pages) |
| `src/st_redirect.py` | Capture solver stdout in UI | **Active** |
| `src/components/*` | Shared fluid page shell | **Active** |
| `src/non_app_plotting.py` | Non-Streamlit plotting stub | **Removed** 2026-11-20 |
| `src/phd_thesis/` | Empty package | **Removed** 2026-11-20 |
| `src/utils/` | Logger + memory helpers for app | **Active** (calculation tab memory) |

---

## Data artefacts

| Path | Role | Recommendation |
|------|------|----------------|
| `data/bulk/*.cr.dat` + `.meta.json` | **Canonical** bulk c(r) | Keep |
| `data/fluid_parameters.toml` | **Canonical** fluid table | Keep |
| `data/fluid_parameters.json` | Old JSON parameters | **Redundant** if unused — prefer TOML only; verify then delete |
| `data/lj1-cr.data` | Legacy LJ1 c(r); still fixture/fallback | Keep until tests + TOML `cr_filename` point only at `data/bulk/lj1.cr.dat` |
| `data/nrcg-cr-approx-maybe.dat` | Old KCl approx dump | **Legacy** — no longer the KCl fixture; optional delete or move under `data/legacy/` |
| `data/Fig3.12_*` | Digitised CF1 figure inputs | Keep (h2o bulk build) |
| `data/pyOZ_bulk_fluid/` | pyOZ test mirror | **Removed** 2026-11-20 — use `src/pyoz/tests` if needed |
| `data/solver_out.txt` | Runtime solver log | Gitignored / local only — fine |
| `data/bulk/` (directory) | Current bulk dir | Keep |

---

## Top-level / tooling noise

| Path | Role | Recommendation |
|------|------|----------------|
| `.coverage` | Local coverage output | Should stay gitignored; delete from disk if present |
| `logs/` | App logs | **Gitignored** 2026-11-20 |
| `output/` | Empty placeholder | Keep or drop |
| `.devcontainer/` (2023) | Old Codespaces config | **Stale** (Python/tooling outdated) — refresh or remove |
| `.vscode/` | Editor settings | Personal; keep if useful |
| `requirements*.txt` | Already removed | Prefer `uv.lock` + `just export-reqs` |

---

## Deployment

| Item | Status |
|------|--------|
| **Render** `PhD-Thesis-MJB` | **Current production** |
| Google Cloud Run | **Deprecated** |
| `phd.databooth.com.au` | Not serving this app (404); optional custom domain on Render later |
| `render.yaml` + `Dockerfile` | Current deploy path |
| `just deploy-render-service` | CLI redeploy helper |

---

## Suggested cleanup PRs (ordered, optional)

1. **Safe deletes:** ~~done 2026-11-20~~ (`src/phd_thesis/`, `src/tests/`, `src/non_app_plotting.py`, `data/pyOZ_bulk_fluid/`, `logs/` ignored).
2. **Retire `oo_refactoring/` + `just app-opt`:** ~~done~~ — history callback ported to `charged_fluids`.
3. **Data hygiene:** move `nrcg-cr-approx-maybe.dat` → `data/legacy/`; point LJ1 tests/`cr_filename` at `data/bulk/lj1.cr.dat` only; drop or archive `fluid_parameters.json`.
4. **Remove experiment trees:** ~~`src/bulk/`~~ done.
5. **Defer `src/pyoz/`** until bulk OZ is fully trusted for all fluids.
6. **Collapse shims:** Streamlit imports `charged_fluids` directly; delete `src/modelling.py` / `parameters.py` / `numerics.py` shims.

---

## Active layout (canonical)

```
src/charged_fluids/     # science core
src/components/         # Streamlit shared UI
src/pages/              # Streamlit multipage
src/Main.py             # entry
data/fluid_parameters.toml
data/bulk/              # c(r) artefacts
tests/                  # pytest
notebooks/*_explorer.py # Marimo
```
