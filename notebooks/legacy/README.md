# Legacy Jupyter notebooks

Archival notebooks from the pre-`charged_fluids` / pre-Marimo era (imports under
`src.*`, pyOZ paths, `numericalunits` experiments).

**Do not extend these.** Active explorers:

- `notebooks/wall_fluid_explorer.py` — wall singlet + bulk inputs (`just marimo`)
- `notebooks/bulk_fluid_explorer.py` — bulk artefacts / parameters / OZ (`just marimo-bulk`)

| Notebook | Historical role | Covered by |
|----------|-----------------|------------|
| `bulk_fluid_explore.ipynb` | Bulk c(r) / potential plots | `bulk_fluid_explorer.py` |
| `fluid_parameters.ipynb` | TOML / parameter dict experiments | `bulk_fluid_explorer.py` + `charged_fluids.parameters` |
| `lj.ipynb` | LJ potential CSV peek | bulk explorer / Streamlit LJ page |
| `check_units.ipynb` | `numericalunits` sandbox | `charged_fluids.units` |
| `numerics.ipynb` | Early `Discretisation` sketch | `charged_fluids.numerics` |
| `plot_convergence.ipynb` | Solver log parser | Streamlit / `src/plotting.py` |
| `singlet_kcl.ipynb` | Long-form KCl wall draft | `wall_fluid_explorer.py` + KCl Streamlit page |
