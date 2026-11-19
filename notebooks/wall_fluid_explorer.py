"""Marimo notebook: charged fluids near interfaces (shared charged_fluids core).

Run:
  just marimo
  # or
  uv sync --group notebook
  uv run marimo edit notebooks/wall_fluid_explorer.py
"""

import marimo

__generated_with = "0.25.2"
app = marimo.App(width="medium")

with app.setup:
    import matplotlib.pyplot as plt


@app.cell
def _intro():
    import marimo as mo

    mo.md(
        r"""
        # Fluids near interfaces — research notebook

        Same **`charged_fluids`** core as the Streamlit app and bulk CLI:

        - TOML parameters (`data/fluid_parameters.toml`)
        - Bulk \(c_\mathrm{short}(r)\) artefacts under `data/bulk/`
        - Singlet wall residual via Newton–Krylov

        Unit conventions: `charged_fluids.units` (Å grids, joule potentials, `charge_pair` scale).
        """
    )
    return (mo,)


@app.cell
def _controls(mo):
    from charged_fluids import list_fluid_symbols

    symbols = [s for s in list_fluid_symbols() if s in {"lj1", "lj2", "kcl", "h2o", "2_2"}]
    symbol = mo.ui.dropdown(options=symbols, value="lj1", label="Fluid")
    n_point = mo.ui.slider(101, 801, value=201, step=50, label="n_point")
    z_cutoff = mo.ui.slider(10.0, 40.0, value=20.0, step=1.0, label="z_cutoff (A)")
    psi_0 = mo.ui.slider(-200.0, 200.0, value=0.0, step=10.0, label="psi_0 (mV)")
    tolerance = mo.ui.number(start=1e-9, stop=1e-3, step=1e-7, value=1e-6, label="tolerance")
    max_iter = mo.ui.slider(5, 100, value=30, step=5, label="max_iteration")
    run_solve = mo.ui.run_button(label="Run wall-fluid solve")
    _ = mo.vstack(
        [
            mo.hstack([symbol, n_point, z_cutoff], justify="start"),
            mo.hstack([psi_0, tolerance, max_iter, run_solve], justify="start"),
        ]
    )
    return max_iter, n_point, psi_0, run_solve, symbol, tolerance, z_cutoff


@app.cell
def _prepare(symbol, n_point, z_cutoff, psi_0, tolerance, max_iter, mo):
    from charged_fluids import prepare_fluid_run

    run = None
    try:
        run = prepare_fluid_run(
            symbol.value,
            z_cutoff=float(z_cutoff.value),
            n_point=int(n_point.value),
            tolerance=float(tolerance.value),
            max_iteration=int(max_iter.value),
            psi_0_mv=float(psi_0.value),
        )
        fluid = run.fluid
        charged = fluid.symbol not in {"lj1", "lj2"}
        lines = [
            f"**{fluid.name}** (`{fluid.symbol}`) · status `{fluid.status}`",
            f"c(r) artefact: `{run.cr_path.name}`",
            f"T = {fluid.temperature} K · conc = {fluid.concentration} M · "
            f"components = {', '.join(fluid.component)}",
            f"n_point = {run.discrete.n_point} · z_cutoff = {run.model.z[-1]:.3g} A · "
            f"psi_0 = {run.psi_0_mv:.3g} mV",
        ]
        if charged:
            lines.append(f"Inverse Debye length kappa = **{run.kappa:.4g}** 1/A")
        msg = "  \n".join(lines)
    except Exception as exc:  # noqa: BLE001 — surface setup errors in the UI
        msg = f"**Setup failed:** `{type(exc).__name__}: {exc}`"
    _ = mo.md(msg)
    return (run,)


@app.cell
def _bulk_metrics(run, mo):
    if run is not None:
        f = run.fluid
        _ = mo.md(
            f"""
            ### Bulk / wall setup metrics

            | Quantity | Value |
            |---|---|
            | Temperature (K) | {f.temperature} |
            | Concentration (M) | {f.concentration} |
            | n_component / n_pair | {f.n_component} / {f.n_pair} |
            | Potential type | `{run.potential.get("type", "?")}` |
            | kappa (1/A) | {run.kappa:.6g} |
            | c_short shape | `{run.model.c_short.shape}` |
            """
        )
    else:
        _ = mo.md("_No fluid run._")
    return


@app.cell
def _bulk_plots(run, mo):
    if run is not None:
        z_grid = run.model.z
        is_charged = run.fluid.symbol not in {"lj1", "lj2"}
        bu_ylim = [0, 10000] if is_charged else [-1, 10]
        r_max = min(10.0, float(z_grid[-1]))

        fig_bulk, axes = plt.subplots(2, 2, figsize=(10, 7))
        ax_c, ax_u, ax_f1, ax_f2 = axes.ravel()

        for pair_i in range(run.model.c_short.shape[1]):
            ax_c.plot(z_grid, run.model.c_short[:, pair_i], label=f"pair {pair_i}")
        ax_c.set_xlim(0, r_max)
        ax_c.set_xlabel("r (A)")
        ax_c.set_ylabel("c_short")
        ax_c.set_title("Bulk short-range c(r)")
        ax_c.legend(fontsize=8)

        bu = run.beta_u
        if bu.ndim == 1:
            ax_u.plot(z_grid, bu)
        else:
            for pair_i in range(bu.shape[1]):
                ax_u.plot(z_grid, bu[:, pair_i], label=f"pair {pair_i}")
            ax_u.legend(fontsize=8)
        ax_u.set_xlim(0, r_max)
        ax_u.set_ylim(*bu_ylim)
        ax_u.set_xlabel("r (A)")
        ax_u.set_ylabel(r"$\beta u$")
        ax_u.set_title("Pair potential")

        for pair_i in range(run.model.f1.shape[1]):
            ax_f1.plot(z_grid, run.model.f1[:, pair_i], label=f"pair {pair_i}")
            ax_f2.plot(z_grid, run.model.f2[:, pair_i], label=f"pair {pair_i}")
        ax_f1.set_xlim(0, r_max)
        ax_f1.set_xlabel("z (A)")
        ax_f1.set_ylabel("f1")
        ax_f1.set_title(r"$f_1(z)=\int_z^\infty z'\,c_\mathrm{short}\,dz'$")
        ax_f1.legend(fontsize=8)
        ax_f2.set_xlim(0, r_max)
        ax_f2.set_xlabel("z (A)")
        ax_f2.set_ylabel("f2")
        ax_f2.set_title(r"$f_2(z)=\int_z^\infty (z')^2 c_\mathrm{short}\,dz'$")
        ax_f2.legend(fontsize=8)

        fig_bulk.tight_layout()
        _ = mo.mpl.interactive(fig_bulk)
    else:
        _ = mo.md("_No run to plot._")
    return


@app.cell
def _solve(run, run_solve, mo):
    from charged_fluids import calc_hw, opt_func, solve_model

    g = None
    hw = None
    sol = None

    if run is None or not run_solve.value:
        _ = mo.md("_Press **Run wall-fluid solve** when ready._")
    elif run.fluid.status not in {"ready", "experimental"}:
        _ = mo.md(f"Fluid status `{run.fluid.status}` — solve skipped.")
    else:
        with mo.status.spinner("Solving wall-fluid residual (Newton–Krylov)…"):
            sol = solve_model(
                opt_func,
                run.model.tw.copy(),
                run.fluid,
                run.model,
                run.discrete,
                run.beta_phiw,
                run.beta_psi_charge,
            )
        hw = calc_hw(sol.x, run.fluid.n_component, run.beta_phiw)
        g = hw + 1.0
        solve_status = "converged" if sol.success else "failed"
        nit = getattr(sol, "nit", "?")
        g_inf = ", ".join(f"{g[-1, i]:.4f}" for i in range(g.shape[1]))
        _ = mo.md(
            f"Solve **{solve_status}** after {nit} iterations.  \n"
            f"Message: `{getattr(sol, 'message', '')}`  \n"
            f"g(z→∞) ≈ [{g_inf}]"
        )
    return g, hw, sol


@app.cell
def _wall_plot(run, g, mo):
    if run is not None and g is not None:
        z_wall = run.model.z
        fig_wall, ax_wall = plt.subplots(figsize=(7, 3.8))
        for comp_i, name in enumerate(run.fluid.component):
            ax_wall.plot(z_wall, g[:, comp_i], label=rf"$g_{{w,{name}}}$")
        ax_wall.axhline(1.0, color="grey", ls="--", lw=0.8)
        ax_wall.set_xlabel("z (A)")
        ax_wall.set_ylabel(r"$g_{wi}(z)$")
        ax_wall.set_xlim(0, min(15.0, float(z_wall[-1])))
        ax_wall.legend()
        ax_wall.set_title("Wall–fluid distribution")
        fig_wall.tight_layout()
        _ = mo.mpl.interactive(fig_wall)
    return


@app.cell
def _bulk_oz_section(mo):
    _ = mo.md(
        r"""
        ---
        ## Optional: regenerate bulk \(c(r)\) (dry-run)

        Uses `charged_fluids.bulk.generate.generate_bulk_cr` with the same path as
        `just bulk-cr` / `python -m charged_fluids.bulk.cli`. Does **not** overwrite
        `data/bulk/` unless you uncheck dry-run (careful).
        """
    )
    bulk_symbol = mo.ui.dropdown(options=["lj1", "lj2", "kcl", "2_2"], value="lj1", label="Bulk fluid")
    bulk_n = mo.ui.slider(256, 2048, value=512, step=256, label="bulk n_points")
    dry_run = mo.ui.checkbox(value=True, label="Dry-run (do not write artefacts)")
    run_bulk = mo.ui.run_button(label="Run bulk OZ")
    _ = mo.hstack([bulk_symbol, bulk_n, dry_run, run_bulk], justify="start")
    return bulk_n, bulk_symbol, dry_run, run_bulk


@app.cell
def _bulk_oz_run(bulk_symbol, bulk_n, dry_run, run_bulk, mo):
    from charged_fluids.bulk.generate import generate_bulk_cr

    result = None
    if not run_bulk.value:
        _ = mo.md("_Press **Run bulk OZ** to solve Picard HNC/PY (may take a few seconds)._")
    else:
        with mo.status.spinner(f"Bulk OZ for {bulk_symbol.value}…"):
            try:
                result, cr_path, meta_path = generate_bulk_cr(
                    bulk_symbol.value,
                    n_points=int(bulk_n.value),
                    write=not bool(dry_run.value),
                )
            except Exception as exc:  # noqa: BLE001
                _ = mo.md(f"**Bulk OZ failed:** `{type(exc).__name__}: {exc}`")
                result = None
            else:
                oz_status = "converged" if result.converged else "NOT converged"
                tail = float(abs(result.c_short[:, -max(5, result.r.size // 20) :]).max())
                paths = ""
                if cr_path is not None:
                    paths = f"  \nWrote `{cr_path}` and `{meta_path}`"
                _ = mo.md(
                    f"**{bulk_symbol.value}**: {oz_status} after {result.n_iter} iterations  \n"
                    f"residual = `{result.residual:.3e}` · "
                    f"closure = `{result.meta.get('closure')}`  \n"
                    f"|c_short| tail max = `{tail:.3e}`{paths}"
                )
    return (result,)


@app.cell
def _bulk_oz_plot(result, mo):
    if result is not None:
        fig_oz, ax_oz = plt.subplots(1, 2, figsize=(10, 3.5))
        for pair_j in range(result.c_short.shape[0]):
            ax_oz[0].plot(result.r, result.c_short[pair_j], label=f"pair {pair_j}")
            ax_oz[1].plot(result.r, result.g[pair_j], label=f"pair {pair_j}")
        ax_oz[0].set_xlabel("r (A)")
        ax_oz[0].set_ylabel("c_short")
        ax_oz[0].set_title("Generated c_short(r)")
        ax_oz[0].set_xlim(0, min(15.0, float(result.r[-1])))
        ax_oz[0].legend(fontsize=8)
        ax_oz[1].set_xlabel("r (A)")
        ax_oz[1].set_ylabel("g(r)")
        ax_oz[1].set_title("g(r)")
        ax_oz[1].set_xlim(0, min(15.0, float(result.r[-1])))
        ax_oz[1].axhline(1.0, color="grey", ls="--", lw=0.8)
        ax_oz[1].legend(fontsize=8)
        fig_oz.tight_layout()
        _ = mo.mpl.interactive(fig_oz)
    return


if __name__ == "__main__":
    app.run()
