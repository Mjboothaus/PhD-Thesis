# Adding and extending fluids

Australian English. Parameter source of truth: `data/fluid_parameters.toml`.

## Status values

| Status | Meaning |
|--------|---------|
| `ready` | Wall + bulk path trusted for demos (LJ1, LJ2) |
| `experimental` | Runnable end-to-end; physics still under active check (KCl, CF1 water, electrolytes) |
| `not_implemented` | Parameters only; no bulk artefact / wall path |

## Minimum steps to add a fluid

1. **TOML section** `[symbol]` with `name`, `symbol`, `component`, `valence`,
   `temperature`, `concentration`, `epsilon_r`, `cr_filename`, `status`.
2. **`[symbol.potential]`** with `type`:
   - `lj` — neutral Lennard-Jones (`epsilon_lj`, `sigma_lj`)
   - `ionic` — soft-core + Coulomb (KCl / 1-1 / 2-2 style)
   - `cf1` — central-force water (O, H); bulk via live Ng/HNC OZ with temperature anneal (digitised Fig 3.12 still opt-in)
3. **`[symbol.bulk]`** grid/closure/Ng knobs as needed.
4. **Bulk artefact** via `just bulk-cr <symbol>` (or `just bulk-cr-h2o` for CF1).
5. **Streamlit page** calling `render_fluid_page("symbol")`.
6. **Tests** — parameters shape, bulk generate, `prepare_fluid_run`, optional wall smoke.
7. **Marimo** — add symbol to the fluid dropdown allow-list if desired.

## Helpers

```python
from charged_fluids import (
    clear_fluid_cache,
    fluids_with_status,
    wall_ready_symbols,
    prepare_fluid_run,
)

clear_fluid_cache()  # after editing TOML in the same process
print(wall_ready_symbols())
run = prepare_fluid_run("1_1", z_cutoff=12.0, n_point=101)
```

## Wall smoke

Experimental fluids should pass a short Newton–Krylov smoke (finite residual step)
even when full convergence is not required. See `tests/test_fluid_wall_smoke.py`.

## Bulk generation CLI

```bash
just bulk-cr-list                 # symbols with bulk support
just bulk-cr lj1                  # TOML defaults for grid/closure
just bulk-cr-grid kcl 512 0.05    # explicit grid
just bulk-cr-all                  # every supported fluid
uv run python -m charged_fluids.bulk.cli --dry-run 1_1
uv run python -m charged_fluids.bulk.cli h2o   # live CF1 Ng/HNC anneal (default)
# Opt-in digitised thesis bulk: set [h2o.bulk] source = "digitised_fig3_12"
```

Success for ionic melts and CF1 is ``converged`` **or** a small ``|c_short|`` tail (usable wall input).
Bare HNC CF1 under-associates vs Duh–Haymet bridge IE (``coord_OH`` ≲ 2); wall physics remains experimental.

