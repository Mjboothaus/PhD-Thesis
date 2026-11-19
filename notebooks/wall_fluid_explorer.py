"""Marimo notebook: charged fluids near interfaces (shared charged_fluids core).

Run:
  uv sync --group notebook
  uv run marimo edit notebooks/wall_fluid_explorer.py
"""

import marimo

__generated_with = "0.13.0"
app = marimo.App(width="medium")


@app.cell
def _intro():
    import marimo as mo

    mo.md(
        r"""
        # Fluids near interfaces - research notebook

        Uses the same **`charged_fluids`** package as the Streamlit app:
        TOML parameters, bulk c(r) artefacts, Newton-Krylov wall solver.
        """
    )
    return (mo,)


@app.cell
def _controls(mo):
    from charged_fluids import list_fluid_symbols

    symbols = [s for s in list_fluid_symbols() if s in {"lj1", "lj2", "kcl"}]
    symbol = mo.ui.dropdown(options=symbols, value="lj1", label="Fluid")
    n_point = mo.ui.slider(101, 801, value=201, step=50, label="n_point")
    z_cutoff = mo.ui.slider(10.0, 40.0, value=20.0, step=1.0, label="z_cutoff (A)")
    run_solve = mo.ui.run_button(label="Run wall-fluid solve")
    mo.hstack([symbol, n_point, z_cutoff, run_solve], justify="start")
    return n_point, run_solve, symbol, z_cutoff


@app.cell
def _prepare(symbol, n_point, z_cutoff, mo):
    from charged_fluids import prepare_fluid_run

    try:
        run = prepare_fluid_run(
            symbol.value,
            z_cutoff=float(z_cutoff.value),
            n_point=int(n_point.value),
            tolerance=1e-6,
            max_iteration=30,
        )
        msg = (
            f"**{run.fluid.name}** · status `{run.fluid.status}` · "
            f"c(r) `{run.cr_path.name}` · kappa={run.kappa:.4g}"
        )
    except Exception as exc:  # noqa: BLE001
        run = None
        msg = f"Setup failed: `{exc}`"
    mo.md(msg)
    return (run,)


@app.cell
def _bulk_plot(run, mo):
    import matplotlib.pyplot as plt

    if run is None:
        mo.md("_No run to plot._")
        return
    z = run.model.z
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.5))
    for j in range(run.model.c_short.shape[1]):
        ax[0].plot(z, run.model.c_short[:, j], label=f"pair {j}")
    ax[0].set_xlim(0, min(10, z[-1]))
    ax[0].set_xlabel("r (A)")
    ax[0].set_ylabel("c_short")
    ax[0].legend()
    ax[0].set_title("Bulk short-range c(r)")

    # beta u may be multi-pair
    bu = run.beta_u
    if bu.ndim == 1:
        ax[1].plot(z, bu)
    else:
        for j in range(bu.shape[1]):
            ax[1].plot(z, bu[:, j], label=f"pair {j}")
        ax[1].legend()
    ax[1].set_xlim(0, min(10, z[-1]))
    ax[1].set_ylim(-2, 20)
    ax[1].set_xlabel("r (A)")
    ax[1].set_ylabel(r"$\beta u$")
    ax[1].set_title("Pair potential")
    fig.tight_layout()
    mo.mpl.interactive(fig)
    return


@app.cell
def _solve(run, run_solve, mo):
    import numpy as np
    from charged_fluids import calc_hw, opt_func, solve_model

    if run is None or not run_solve.value:
        mo.md("_Press **Run wall-fluid solve** when ready._")
        return (None,)

    if run.fluid.status not in {"ready", "experimental"}:
        mo.md(f"Fluid status `{run.fluid.status}` - solve skipped.")
        return (None,)

    sol = solve_model(
        opt_func,
        run.model.tw,
        run.fluid,
        run.model,
        run.discrete,
        run.beta_phiw,
        run.beta_psi_charge,
    )
    hw = calc_hw(sol.x, run.fluid.n_component, run.beta_phiw)
    g = hw + 1.0
    status = "converged" if sol.success else "failed"
    mo.md(
        f"Solve **{status}** after {getattr(sol, 'nit', '?')} iterations. "
        f"g(z->∞) ≈ {g[-1, 0]:.4f}"
    )
    return (g,)


@app.cell
def _wall_plot(run, g, mo):
    import matplotlib.pyplot as plt

    if run is None or g is None:
        return
    z = run.model.z
    fig, ax = plt.subplots(figsize=(6, 3.5))
    for i, name in enumerate(run.fluid.component):
        ax.plot(z, g[:, i], label=rf"$g_{{w,{name}}}$")
    ax.set_xlabel("z (A)")
    ax.set_ylabel("g_wi(z)")
    ax.set_xlim(0, min(15, z[-1]))
    ax.legend()
    ax.set_title("Wall-fluid radial distribution")
    fig.tight_layout()
    mo.mpl.interactive(fig)
    return


if __name__ == "__main__":
    app.run()
