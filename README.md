# Charged fluids near interfaces: Integral equation theory (Ph.D. Thesis)

Modern Python reimplementation of the 1998 Ph.D. thesis *Charged fluids near interfaces: Integral equation theory*, with a UI-free science core, Streamlit app, and bulk Ornstein–Zernike pipeline.

## Web app

Production (Render): [https://phd-thesis-mjb.onrender.com](https://phd-thesis-mjb.onrender.com)  
Deploy: Docker via `render.yaml` / `Dockerfile` (`just deploy-render-service`).  
**Google Cloud Run is deprecated** for this app and is no longer the production target.

## Background

The original thesis included electronic resources at the University of Sydney (FORTRAN-era materials are no longer hosted). This repository rebuilds the wall-fluid integral-equation workflow and bulk \(c(r)\) generation in a tested, uv-managed stack.

![Online Resources cover page](https://github.com/Mjboothaus/PhD-Thesis/blob/deb94173c8ee3b2591c19219e2ea813f33181ae1/docs/Thesis_OnlineResources_Cover.jpg)

## Status

| Fluid | Wall solve | Bulk \(c_\mathrm{short}\) artefact | Notes |
|-------|------------|--------------------------------------|--------|
| **lj1** | Ready | `data/bulk/lj1.cr.dat` | One-component LJ |
| **lj2** | Ready | `data/bulk/lj2.cr.dat` | Two-component LJ |
| **kcl** | Experimental | `data/bulk/kcl.cr.dat` (Ng/PY + T anneal) | Melt bulk residual small with anneal; soft LJ+Coulomb core |
| **h2o** | Experimental | `data/bulk/h2o.cr.dat` (CF1 Ng/HNC anneal) | Bridge IE open; digitised Fig 3.12 opt-in |
| **2_2** | Experimental | `data/bulk/2_2.cr.dat` | Soft-core 2:2 (Ng/PY + anneal) |
| **1_1** | Experimental | `data/bulk/1_1.cr.dat` | Soft-core 1:1 (Ng/HNC + anneal) |

Core package: **`src/charged_fluids/`** (Python 3.14). Streamlit pages under `src/pages/` are thin shims. Shared preparation: `prepare_fluid_run` → `FluidRun`.

Architecture (one line): TOML parameters → bulk Picard OZ (HNC/PY, Ng, optional T anneal) → wall singlet residual (Newton–Krylov) → Streamlit / Marimo.

Physics status: **`docs/PHYSICS.md`**.

## Development (uv + Just)

```bash
just init            # uv sync (Python 3.14)
just test            # pytest
just quality         # ruff + ty + pytest
just app             # Streamlit (localhost)
just bulk-cr lj1     # write data/bulk/{symbol}.cr.dat (+ .meta.json; TOML defaults)
just bulk-cr-all     # all bulk-capable fluids
just bulk-cr-list    # list bulk symbols
just marimo          # wall–fluid Marimo analysis explorer (Plotly)
just marimo-bulk     # bulk c(r) / parameters explorer (Plotly)
just hooks-install   # prek git hooks
```

Bulk CLI (TOML supplies defaults for points/step/closure when flags are omitted):

```bash
uv run python -m charged_fluids.bulk.cli kcl
```

## Package layout

```
src/charged_fluids/
  modelling.py      # wall kernels, calc_tw, pair index
  setup.py          # FluidRun / prepare_fluid_run
  parameters.py     # TOML fluid table
  units.py          # unit conventions (Å, charge_pair, …)
  bulk/
    oz_solver.py    # Picard OZ (LJ, binary Coulomb+LJ + Ng)
    closures.py     # HNC / PY registry
    transforms.py   # FFT sine transforms (direct fallback)
    cli.py          # artefact generation CLI
data/fluid_parameters.toml   # single source of truth for fluids
data/bulk/                  # generated c(r) artefacts
tests/                      # pytest suite for charged_fluids
docs/blog/                  # progressive notes + capstone
docs/LEGACY.md              # historical trees (do not expand)
docs/FLUIDS.md               # how to add a fluid
```

## Bulk correlation heritage

Historical bulk solvers: [pyOZ](http://pyoz.vrbka.net) (discontinued), fork [ctk3b/pyoz](https://github.com/ctk3b/pyoz), and [elvissoares/PyOZ](https://github.com/elvissoares/PyOZ) (untested). Production path is **`charged_fluids.bulk`**; vendored pyOZ remains under `src/pyoz/` for reference (`just bulk-fluid-pyoz`).

## Documentation and blog

- Capstone before/after: `docs/blog/2026-11-19-capstone-before-after.md`
- Other notes: `docs/blog/`
- Legacy paths: `docs/LEGACY.md`
- Agent-oriented project notes: `WARP.md` (may lag; prefer this README + `charged_fluids` for current behaviour)

## Open work

- CF1 **bridge** bulk IE (coord_OH → 2); digitise thesis **Fig 3.2** wall curves for regression
- Full Pauling ionic core in bulk OZ (beyond soft LJ+Coulomb); melt g(r) vs MC
- Literature table benchmarks (Hansen–McDonald LJ; modern closures)
- Electrolyte thermodynamics (pressure / energy routes) beyond c_short artefacts
- Optional: collapse Streamlit import shims; custom domain on Render

## License

MIT
