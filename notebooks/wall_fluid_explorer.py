"""Marimo notebook: charged fluids near interfaces (shared charged_fluids core).

Run:
  just marimo
  # or
  uv sync --group notebook
  uv run marimo edit notebooks/wall_fluid_explorer.py

Display rule (marimo): the last *expression* in a cell is what appears in the UI.
Do not assign the only visual to a throwaway name without referring to it again.

Figures use Plotly (interactive zoom/pan). One fluid selection drives wall and bulk.
"""

from __future__ import annotations

import marimo

__generated_with = "0.25.2"
app = marimo.App(width="full", app_title="Charged fluids near interfaces")


with app.setup:
    import numpy as np
    import plotly.graph_objects as go
    import plotly.io as pio
    from plotly.subplots import make_subplots

    pio.templates.default = "plotly_white"

    # Fluids that can appear in the shared selector
    FLUID_OPTIONS = ("lj1", "lj2", "kcl", "h2o", "2_2", "1_1")
    # Bulk OZ regenerate is not defined for digitised CF1 water
    BULK_OZ_SYMBOLS = frozenset({"lj1", "lj2", "kcl", "2_2", "1_1"})


def _pair_labels(components: list[str]) -> list[str]:
    labels: list[str] = []
    for i, a in enumerate(components):
        for b in components[i:]:
            labels.append(f"{a}–{b}" if a != b else f"{a}–{a}")
    return labels


def _series_from_array(x, y_2d):
    """Yield (name_suffix_index, y) for each pair series."""
    y = np.asarray(y_2d)
    if y.ndim == 1:
        yield 0, y
        return
    if y.shape[0] == len(x):
        for i in range(y.shape[1]):
            yield i, y[:, i]
    else:
        for i in range(y.shape[0]):
            yield i, y[i]


def _add_series(fig, x, y_2d, *, row=None, col=None, labels=None, legendgroup=None):
    for i, yi in _series_from_array(x, y_2d):
        name = labels[i] if labels and i < len(labels) else f"pair {i}"
        tr = go.Scatter(
            x=x,
            y=yi,
            mode="lines",
            name=name,
            legendgroup=legendgroup or name,
            hovertemplate="%{x:.3g} Å<br>%{y:.4g}<extra>" + name + "</extra>",
        )
        if row is None:
            fig.add_trace(tr)
        else:
            fig.add_trace(tr, row=row, col=col)


def _style_fig(fig, *, height=420, title=None):
    fig.update_layout(
        title=dict(text=title, x=0.02, xanchor="left") if title else None,
        height=height,
        margin=dict(l=56, r=24, t=56 if title else 36, b=48),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1.0,
            font=dict(size=11),
        ),
        hovermode="x unified",
    )
    fig.update_xaxes(showgrid=True, zeroline=False)
    fig.update_yaxes(showgrid=True, zeroline=False)
    return fig


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell
def _(mo):
    from charged_fluids import list_fluid_symbols

    available = [s for s in list_fluid_symbols() if s in FLUID_OPTIONS]
    if not available:
        available = list(FLUID_OPTIONS)

    symbol = mo.ui.dropdown(
        options=available,
        value="lj1" if "lj1" in available else available[0],
        label="Fluid",
    )
    n_point = mo.ui.slider(101, 801, value=201, step=50, label="Wall grid points")
    z_cutoff = mo.ui.slider(10.0, 40.0, value=20.0, step=1.0, label="z cutoff (Å)")
    psi_0 = mo.ui.slider(-200.0, 200.0, value=0.0, step=10.0, label="ψ₀ (mV)")
    tolerance = mo.ui.number(
        start=1e-9, stop=1e-3, step=1e-7, value=1e-6, label="Solver tolerance"
    )
    max_iter = mo.ui.slider(5, 100, value=30, step=5, label="Max iterations")
    bulk_n = mo.ui.slider(256, 2048, value=512, step=256, label="Bulk OZ points")
    dry_run = mo.ui.checkbox(value=True, label="Dry-run bulk write")
    run_solve = mo.ui.run_button(label="Solve wall–fluid")
    run_bulk = mo.ui.run_button(label="Run bulk OZ")

    header = mo.vstack(
        [
            mo.md(
                r"""
                # Fluids near interfaces
                Shared **`charged_fluids`** core · Plotly figures · one fluid for wall **and** bulk
                """
            ),
            mo.callout(
                mo.md(
                    "Choose a **fluid** once below. Wall solve and optional bulk OZ both "
                    "use that selection. Prefer dry-run when regenerating \(c(r)\)."
                ),
                kind="info",
            ),
        ],
        gap=0.5,
    )
    header
    return (
        bulk_n,
        dry_run,
        max_iter,
        n_point,
        psi_0,
        run_bulk,
        run_solve,
        symbol,
        tolerance,
        z_cutoff,
    )


@app.cell
def _(
    bulk_n,
    dry_run,
    max_iter,
    mo,
    n_point,
    psi_0,
    run_bulk,
    run_solve,
    symbol,
    tolerance,
    z_cutoff,
):
    controls = mo.vstack(
        [
            mo.md("### Fluid & numerics"),
            mo.hstack(
                [symbol, n_point, z_cutoff, psi_0],
                justify="start",
                gap=1.0,
                wrap=True,
            ),
            mo.hstack(
                [tolerance, max_iter, bulk_n, dry_run],
                justify="start",
                gap=1.0,
                wrap=True,
            ),
            mo.hstack([run_solve, run_bulk], justify="start", gap=1.0),
        ],
        gap=0.75,
    )
    controls
    return


@app.cell
def _(max_iter, mo, n_point, psi_0, symbol, tolerance, z_cutoff):
    from charged_fluids import prepare_fluid_run

    run = None
    err = None
    try:
        run = prepare_fluid_run(
            symbol.value,
            z_cutoff=float(z_cutoff.value),
            n_point=int(n_point.value),
            tolerance=float(tolerance.value),
            max_iteration=int(max_iter.value),
            psi_0_mv=float(psi_0.value),
        )
    except Exception as exc:  # noqa: BLE001 — surface in UI
        err = f"{type(exc).__name__}: {exc}"

    if err is not None:
        summary = mo.callout(mo.md(f"**Setup failed** — `{err}`"), kind="danger")
    else:
        f = run.fluid
        charged = f.symbol not in {"lj1", "lj2"}
        stats = [
            mo.stat(label="Fluid", value=f"{f.name}", caption=f"`{f.symbol}` · {f.status}"),
            mo.stat(
                label="T / conc",
                value=f"{f.temperature:g} K",
                caption=f"{f.concentration} M · {', '.join(f.component)}",
            ),
            mo.stat(
                label="Grid",
                value=f"{run.discrete.n_point}",
                caption=f"z ≤ {run.model.z[-1]:.3g} Å · ψ₀ {run.psi_0_mv:g} mV",
            ),
            mo.stat(
                label="κ" if charged else "Potential",
                value=f"{run.kappa:.4g} Å⁻¹" if charged else str(run.potential.get("type", "—")),
                caption=f"c_short {tuple(run.model.c_short.shape)} · `{run.cr_path.name}`",
            ),
        ]
        summary = mo.hstack(stats, justify="space-between", gap=1.0, wrap=True)
    summary
    return (run,)


@app.cell
def _(mo, run):
    if run is None:
        bulk_fig = mo.md("_Fix setup errors to see bulk inputs._")
    else:
        z = run.model.z
        labels = _pair_labels(list(run.fluid.component))
        is_charged = run.fluid.symbol not in {"lj1", "lj2"}
        bu_ylim = [0, 10000] if is_charged else [-1, 10]
        r_max = min(10.0, float(z[-1]))

        fig = make_subplots(
            rows=2,
            cols=2,
            subplot_titles=(
                r"Bulk \(c_\mathrm{short}(r)\)",
                r"Pair potential \(\beta u\)",
                r"\(f_1(z)\)",
                r"\(f_2(z)\)",
            ),
            vertical_spacing=0.14,
            horizontal_spacing=0.10,
        )
        _add_series(fig, z, run.model.c_short, row=1, col=1, labels=labels)
        _add_series(fig, z, run.beta_u, row=1, col=2, labels=labels)
        _add_series(fig, z, run.model.f1, row=2, col=1, labels=labels)
        _add_series(fig, z, run.model.f2, row=2, col=2, labels=labels)

        fig.update_xaxes(range=[0, r_max], title_text="r, z (Å)")
        fig.update_yaxes(title_text="c_short", row=1, col=1)
        fig.update_yaxes(title_text="βu", range=bu_ylim, row=1, col=2)
        fig.update_yaxes(title_text="f₁", row=2, col=1)
        fig.update_yaxes(title_text="f₂", row=2, col=2)
        # avoid duplicate legend entries across subplots
        names_seen: set[str] = set()
        for tr in fig.data:
            if tr.name in names_seen:
                tr.showlegend = False
            else:
                names_seen.add(tr.name)
        _style_fig(fig, height=560, title=f"Bulk inputs · {run.fluid.symbol}")
        bulk_fig = mo.ui.plotly(fig)
    bulk_fig
    return


@app.cell
def _(mo, run, run_solve):
    from charged_fluids import calc_hw, opt_func, solve_model

    g = None
    solve_banner = mo.md("")

    if run is None:
        solve_banner = mo.md("_Wall solve unavailable until setup succeeds._")
    elif not run_solve.value:
        solve_banner = mo.callout(
            mo.md("Press **Solve wall–fluid** to run Newton–Krylov on the singlet residual."),
            kind="neutral",
        )
    elif run.fluid.status not in {"ready", "experimental"}:
        solve_banner = mo.callout(
            mo.md(f"Fluid status `{run.fluid.status}` — solve skipped."),
            kind="warn",
        )
    else:
        with mo.status.spinner("Solving wall–fluid residual…"):
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
        ok = bool(sol.success)
        nit = getattr(sol, "nit", "?")
        g_inf = ", ".join(f"{g[-1, i]:.4f}" for i in range(g.shape[1]))
        kind = "success" if ok else "danger"
        solve_banner = mo.callout(
            mo.md(
                f"**{'Converged' if ok else 'Failed'}** after {nit} iterations.  \n"
                f"`{getattr(sol, 'message', '')}`  \n"
                f"\(g(z\\to\\infty) \\approx [{g_inf}]\)"
            ),
            kind=kind,
        )
    solve_banner
    return (g,)


@app.cell
def _(g, mo, run):
    if run is None or g is None:
        wall_fig = mo.md(r"_Wall \(g_{wi}(z)\) will appear here after a solve._")
    else:
        z = run.model.z
        fig = go.Figure()
        for i, name in enumerate(run.fluid.component):
            fig.add_trace(
                go.Scatter(
                    x=z,
                    y=g[:, i],
                    mode="lines",
                    name=rf"g<sub>w,{name}</sub>",
                    hovertemplate="%{x:.3g} Å<br>%{y:.4g}<extra>" + name + "</extra>",
                )
            )
        fig.add_hline(y=1.0, line_dash="dot", line_color="#888", line_width=1)
        fig.update_xaxes(title_text="z (Å)", range=[0, min(15.0, float(z[-1]))])
        fig.update_yaxes(title_text=r"g<sub>wi</sub>(z)")
        _style_fig(fig, height=400, title=f"Wall–fluid distribution · {run.fluid.symbol}")
        wall_fig = mo.ui.plotly(fig)
    wall_fig
    return


@app.cell
def _(bulk_n, dry_run, mo, run_bulk, symbol):
    from charged_fluids.bulk.generate import generate_bulk_cr

    result = None
    sym = symbol.value

    if sym not in BULK_OZ_SYMBOLS:
        oz_banner = mo.callout(
            mo.md(
                f"**`{sym}`** has no live bulk OZ path here "
                f"(e.g. digitised CF1). Use the stored artefact above, "
                f"or switch fluid to one of `{', '.join(sorted(BULK_OZ_SYMBOLS))}`."
            ),
            kind="warn",
        )
    elif not run_bulk.value:
        mode = "dry-run" if dry_run.value else "**will write** `data/bulk/`"
        oz_banner = mo.callout(
            mo.md(
                f"Optional: regenerate bulk \(c(r)\) for **`{sym}`** "
                f"({int(bulk_n.value)} points, {mode}). "
                "Same entry point as `just bulk-cr`."
            ),
            kind="neutral",
        )
    else:
        with mo.status.spinner(f"Bulk OZ · {sym}…"):
            try:
                result, cr_path, meta_path = generate_bulk_cr(
                    sym,
                    n_points=int(bulk_n.value),
                    write=not bool(dry_run.value),
                )
            except Exception as exc:  # noqa: BLE001
                oz_banner = mo.callout(
                    mo.md(f"**Bulk OZ failed:** `{type(exc).__name__}: {exc}`"),
                    kind="danger",
                )
                result = None
            else:
                ok = bool(result.converged)
                tail = float(abs(result.c_short[:, -max(5, result.r.size // 20) :]).max())
                paths = ""
                if cr_path is not None:
                    paths = f"  \nWrote `{cr_path}` · `{meta_path}`"
                oz_banner = mo.callout(
                    mo.md(
                        f"**{sym}**: {'converged' if ok else 'not converged'} "
                        f"· {result.n_iter} iters · residual `{result.residual:.3e}`  \n"
                        f"closure `{result.meta.get('closure')}` · "
                        f"|c_short| tail max `{tail:.3e}`{paths}"
                    ),
                    kind="success" if ok else "warn",
                )
    oz_banner
    return (result,)


@app.cell
def _(mo, result, run, symbol):
    if result is None:
        oz_fig = mo.md("")
    else:
        labels = None
        if run is not None and run.fluid.symbol == symbol.value:
            labels = _pair_labels(list(run.fluid.component))
        r_max = min(15.0, float(result.r[-1]))
        fig = make_subplots(
            rows=1,
            cols=2,
            subplot_titles=(r"Generated \(c_\mathrm{short}(r)\)", r"\(g(r)\)"),
            horizontal_spacing=0.10,
        )
        _add_series(fig, result.r, result.c_short, row=1, col=1, labels=labels)
        _add_series(fig, result.r, result.g, row=1, col=2, labels=labels)
        fig.add_hline(y=1.0, line_dash="dot", line_color="#888", line_width=1, row=1, col=2)
        names_seen: set[str] = set()
        for tr in fig.data:
            if tr.name in names_seen:
                tr.showlegend = False
            else:
                names_seen.add(tr.name)
        fig.update_xaxes(range=[0, r_max], title_text="r (Å)")
        fig.update_yaxes(title_text="c_short", row=1, col=1)
        fig.update_yaxes(title_text="g(r)", row=1, col=2)
        _style_fig(fig, height=380, title=f"Bulk OZ · {symbol.value}")
        oz_fig = mo.ui.plotly(fig)
    oz_fig
    return


@app.cell
def _(mo):
    footer = mo.accordion(
        {
            "Units & companions": mo.md(
                r"""
                - Length: Å · concentration: mol L⁻¹ · potentials → \(\beta u\)
                - Coulomb scale: `charge_pair / (2 r)` (see `charged_fluids.units`)
                - Bulk-only explorer: `just marimo-bulk`
                - Legacy Jupyter: `notebooks/legacy/`
                """
            )
        }
    )
    footer
    return


if __name__ == "__main__":
    app.run()
