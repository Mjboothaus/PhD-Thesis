r"""Marimo notebook: wall–fluid singlet analysis (shared charged_fluids core).

Interactive explorer for integral-equation wall profiles \(g_{wi}(z)\), bulk
inputs \(c_\mathrm{short}\), pair potentials, and optional bulk OZ regenerate.

Run:
  just marimo
  # or
  uv sync --group notebook
  uv run marimo edit notebooks/wall_fluid_explorer.py

Display rule (marimo): the last *expression* in a cell is what appears in the UI.
Figures use Plotly (interactive zoom/pan). One fluid selection drives wall and bulk.
"""

from __future__ import annotations

import marimo

__generated_with = "0.25.2"
app = marimo.App(width="full", app_title="Wall–fluid analysis · charged fluids")


with app.setup:
    import numpy as np
    import plotly.graph_objects as go
    import plotly.io as pio
    from plotly.subplots import make_subplots

    pio.templates.default = "plotly_white"

    # Fluids offered in the shared selector (must also exist in TOML)
    FLUID_OPTIONS = ("lj1", "lj2", "kcl", "h2o", "2_2", "1_1")


def _pair_labels(components: list[str]) -> list[str]:
    labels: list[str] = []
    for i, a in enumerate(components):
        for b in components[i:]:
            labels.append(f"{a}–{b}" if a != b else f"{a}–{a}")
    return labels


def _series_from_array(x, y_2d):
    """Yield (series_index, y) for each pair column."""
    y = np.asarray(y_2d)
    if y.ndim == 1:
        yield 0, y
        return
    # Wall loaders store c_short as (n_z, n_pair); bulk OZ as (n_pair, n_r)
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
            x=np.asarray(x, dtype=float),
            y=np.asarray(yi, dtype=float),
            mode="lines",
            name=name,
            legendgroup=legendgroup or name,
            hovertemplate="%{x:.3g} Å<br>%{y:.4g}<extra>" + name + "</extra>",
        )
        if row is None:
            fig.add_trace(tr)
        else:
            fig.add_trace(tr, row=row, col=col)


def _dedupe_legend(fig) -> None:
    seen: set[str] = set()
    for tr in fig.data:
        name = getattr(tr, "name", None) or ""
        if name in seen:
            tr.showlegend = False
        else:
            seen.add(name)


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


def _contact_and_bulk_g(z: np.ndarray, g: np.ndarray) -> list[dict]:
    """Per-component contact (first finite g>0) and far-field values."""
    rows: list[dict] = []
    g = np.asarray(g, dtype=float)
    z = np.asarray(z, dtype=float)
    if g.ndim == 1:
        g = g.reshape(-1, 1)
    n_comp = g.shape[1]
    for i in range(n_comp):
        gi = g[:, i]
        finite = np.isfinite(gi) & (gi > 1.0e-12)
        if finite.any():
            j = int(np.argmax(finite))
            g_c = float(gi[j])
            z_c = float(z[j])
        else:
            g_c, z_c = float("nan"), float("nan")
        near = z <= min(5.0, float(z[-1]))
        if near.any() and np.isfinite(gi[near]).any():
            jpk = int(np.nanargmax(np.where(near, gi, np.nan)))
            g_pk, z_pk = float(gi[jpk]), float(z[jpk])
        else:
            g_pk, z_pk = float("nan"), float("nan")
        rows.append(
            {
                "g_contact": g_c,
                "z_contact": z_c,
                "g_peak": g_pk,
                "z_peak": z_pk,
                "g_inf": float(gi[-1]) if gi.size else float("nan"),
            }
        )
    return rows


def _reshape_tw(tw: np.ndarray, n_z: int, n_component: int) -> np.ndarray:
    """Return tw as (n_z, n_component)."""
    tw = np.asarray(tw, dtype=float)
    if tw.ndim == 1:
        if tw.size == n_z * n_component:
            return tw.reshape(n_component, n_z).T
        return tw.reshape(n_z, n_component)
    if tw.shape == (n_z, n_component):
        return tw
    if tw.shape == (n_component, n_z):
        return tw.T
    return tw.reshape(n_z, n_component)


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell
def _(mo):
    from charged_fluids import list_fluid_symbols
    from charged_fluids.bulk.generate import list_bulk_symbols

    available = [s for s in list_fluid_symbols() if s in FLUID_OPTIONS]
    if not available:
        available = list(FLUID_OPTIONS)
    bulk_symbols = frozenset(list_bulk_symbols())

    symbol = mo.ui.dropdown(
        options=available,
        value="lj1" if "lj1" in available else available[0],
        label="Fluid",
    )
    n_point = mo.ui.slider(81, 801, value=201, step=40, label="Wall grid points")
    z_cutoff = mo.ui.slider(8.0, 40.0, value=20.0, step=1.0, label="z cutoff (Å)")
    z_plot_max = mo.ui.slider(4.0, 25.0, value=12.0, step=0.5, label="Plot z_max (Å)")
    psi_0 = mo.ui.slider(-200.0, 200.0, value=0.0, step=10.0, label="ψ₀ (mV)")
    tolerance = mo.ui.number(
        start=1e-9, stop=1e-3, step=1e-7, value=1e-6, label="Solver tolerance"
    )
    max_iter = mo.ui.slider(5, 120, value=40, step=5, label="Max iterations")
    bulk_n = mo.ui.slider(256, 2048, value=512, step=256, label="Bulk OZ points")
    dry_run = mo.ui.checkbox(value=True, label="Dry-run bulk write")
    run_solve = mo.ui.run_button(label="Solve wall–fluid")
    run_bulk = mo.ui.run_button(label="Run bulk OZ")

    header = mo.vstack(
        [
            mo.md(
                r"""
                # Wall–fluid analysis
                Shared **`charged_fluids`** core · Plotly · one fluid for wall **and** bulk

                Newton–Krylov singlet residual \(\delta t_w(z)=0\) → wall–species
                distribution functions \(g_{wi}(z)\). Bulk \(c_\mathrm{short}(r)\) is
                loaded from `data/bulk/{symbol}.cr.dat` (same path as Streamlit).
                """
            ),
            mo.callout(
                mo.md(
                    r"Choose a **fluid** once. Press **Solve wall–fluid** for profiles; "
                    r"optional **Run bulk OZ** regenerates \(c(r)\) (prefer dry-run). "
                    r"CF1 water (`h2o`) has a live Ng/HNC bulk path with temperature anneal."
                ),
                kind="info",
            ),
        ],
        gap=0.5,
    )
    header
    return (
        bulk_n,
        bulk_symbols,
        dry_run,
        max_iter,
        n_point,
        psi_0,
        run_bulk,
        run_solve,
        symbol,
        tolerance,
        z_cutoff,
        z_plot_max,
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
    z_plot_max,
):
    controls = mo.vstack(
        [
            mo.md("### Fluid & numerics"),
            mo.hstack(
                [symbol, n_point, z_cutoff, z_plot_max, psi_0],
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
            mo.stat(
                label="Fluid",
                value=f"{f.name}",
                caption=f"`{f.symbol}` · {f.status}",
            ),
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
                value=(
                    f"{run.kappa:.4g} Å⁻¹"
                    if charged
                    else str(run.potential.get("type", "—"))
                ),
                caption=f"c_short {tuple(run.model.c_short.shape)} · `{run.cr_path.name}`",
            ),
        ]
        blocks = [mo.hstack(stats, justify="space-between", gap=1.0, wrap=True)]
        if f.status == "experimental":
            blocks.append(
                mo.callout(
                    mo.md(
                        f"**Experimental** fluid `{f.symbol}`: wall solve will run, but "
                        "bulk thermodynamics / residual quality may still be imperfect "
                        "(KCl melt, CF1 bare HNC, soft-core electrolytes)."
                    ),
                    kind="warn",
                )
            )
        summary = mo.vstack(blocks, gap=0.5)
    summary
    return (run,)


@app.cell
def _(mo, run, z_plot_max):
    if run is None:
        bulk_fig = mo.md("_Fix setup errors to see bulk inputs._")
    else:
        z = run.model.z
        labels = _pair_labels(list(run.fluid.component))
        is_charged = run.fluid.symbol not in {"lj1", "lj2"}
        bu_ylim = [0, 10000] if is_charged else [-1, 10]
        r_max = min(float(z_plot_max.value), float(z[-1]))

        fig = make_subplots(
            rows=2,
            cols=2,
            subplot_titles=(
                r"Bulk \(c_\mathrm{short}\) (wall grid)",
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
        _dedupe_legend(fig)
        _style_fig(fig, height=560, title=f"Bulk inputs · {run.fluid.symbol}")
        bulk_fig = mo.ui.plotly(fig)
    bulk_fig
    return


@app.cell
def _(mo, run, run_solve):
    from charged_fluids import calc_hw, opt_func, solve_model

    g = None
    tw = None
    sol = None
    solve_banner = mo.md("")

    if run is None:
        solve_banner = mo.md("_Wall solve unavailable until setup succeeds._")
    elif not run_solve.value:
        solve_banner = mo.callout(
            mo.md(
                r"Press **Solve wall–fluid** to run Newton–Krylov on the singlet residual "
                r"\(F[t_w]=0\)."
            ),
            kind="neutral",
        )
    elif run.fluid.status not in {"ready", "experimental"}:
        solve_banner = mo.callout(
            mo.md(f"Fluid status `{run.fluid.status}` — solve skipped."),
            kind="warn",
        )
    else:
        with mo.status.spinner("Solving wall–fluid residual…"):
            try:
                sol = solve_model(
                    opt_func,
                    run.model.tw.copy(),
                    run.fluid,
                    run.model,
                    run.discrete,
                    run.beta_phiw,
                    run.beta_psi_charge,
                )
            except Exception as exc:  # noqa: BLE001
                solve_banner = mo.callout(
                    mo.md(f"**Solver raised:** `{type(exc).__name__}: {exc}`"),
                    kind="danger",
                )
                sol = None
            else:
                tw = np.asarray(sol.x, dtype=float)
                hw = calc_hw(tw, run.fluid.n_component, run.beta_phiw)
                g = hw + 1.0
                ok = bool(sol.success)
                nit = getattr(sol, "nit", "?")
                nfev = getattr(sol, "nfev", "?")
                fun = np.asarray(getattr(sol, "fun", []), dtype=float)
                f_norm = float(np.linalg.norm(fun)) if fun.size else float("nan")
                g_inf = ", ".join(f"{g[-1, i]:.4f}" for i in range(g.shape[1]))
                kind = "success" if ok else "warn"
                solve_banner = mo.callout(
                    mo.md(
                        f"**{'Converged' if ok else 'Finished (not fully converged)'}** "
                        f"after **{nit}** iterations · **{nfev}** residual evals.  \n"
                        f"`{getattr(sol, 'message', '')}`  \n"
                        rf"\(\|F\|_2 = {f_norm:.4e}\) · "
                        rf"\(g(z\to\infty) \approx [{g_inf}]\)"
                    ),
                    kind=kind,
                )
    solve_banner
    return g, sol, tw


@app.cell
def _(g, mo, run, tw, z_plot_max):
    if run is None or g is None:
        wall_block = mo.md(r"_Wall \(g_{wi}(z)\) will appear here after a solve._")
    else:
        z = run.model.z
        z_max = min(float(z_plot_max.value), float(z[-1]))
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
        fig.update_xaxes(title_text="z (Å)", range=[0, z_max])
        fig.update_yaxes(title_text=r"g<sub>wi</sub>(z)")
        _style_fig(fig, height=420, title=f"Wall–fluid distribution · {run.fluid.symbol}")

        rows = _contact_and_bulk_g(z, g)
        table_lines = [
            "| Species | z_contact (Å) | g_contact | z_peak (Å) | g_peak | g(∞) |",
            "|---|---:|---:|---:|---:|---:|",
        ]
        for i, row in enumerate(rows):
            sp = run.fluid.component[i] if i < len(run.fluid.component) else str(i)
            table_lines.append(
                f"| {sp} | {row['z_contact']:.3g} | {row['g_contact']:.4g} | "
                f"{row['z_peak']:.3g} | {row['g_peak']:.4g} | {row['g_inf']:.4g} |"
            )
        diagnostics = mo.md(
            "### Wall profile diagnostics\n\n"
            + "\n".join(table_lines)
            + "\n\n"
            + r"Contact = first grid point with \(g>0\); peak scanned for \(z\le 5\) Å. "
            + r"Ideal bulk limit is \(g\to 1\)."
        )

        tw_fig = mo.md("")
        if tw is not None:
            tw_arr = _reshape_tw(tw, len(z), run.fluid.n_component)
            tfig = go.Figure()
            for i, name in enumerate(run.fluid.component):
                tfig.add_trace(
                    go.Scatter(
                        x=z,
                        y=tw_arr[:, i],
                        mode="lines",
                        name=rf"t<sub>w,{name}</sub>",
                    )
                )
            tfig.update_xaxes(title_text="z (Å)", range=[0, z_max])
            tfig.update_yaxes(title_text=r"t<sub>w</sub>(z)")
            _style_fig(tfig, height=320, title=r"Indirect wall correlation \(t_w\)")
            tw_fig = mo.ui.plotly(tfig)

        wall_block = mo.vstack(
            [mo.ui.plotly(fig), diagnostics, tw_fig],
            gap=0.75,
        )
    wall_block
    return


@app.cell
def _(bulk_n, bulk_symbols, dry_run, mo, run_bulk, symbol):
    from charged_fluids.bulk.generate import generate_bulk_cr

    result = None
    sym = symbol.value

    if sym not in bulk_symbols:
        oz_banner = mo.callout(
            mo.md(
                f"**`{sym}`** is not in the bulk-generation allow-list "
                f"(`list_bulk_symbols()` → `{', '.join(sorted(bulk_symbols))}`). "
                "Wall analysis still uses any stored artefact under `data/bulk/`."
            ),
            kind="warn",
        )
    elif not run_bulk.value:
        mode = "dry-run" if dry_run.value else "**will write** `data/bulk/`"
        note_h2o = ""
        if sym == "h2o":
            note_h2o = (
                " CF1 uses Ng/HNC **temperature annealing** (can take tens of seconds "
                "at the full TOML grid; fewer bulk points are faster for exploration)."
            )
        oz_banner = mo.callout(
            mo.md(
                rf"Optional: regenerate bulk \(c(r)\) for **`{sym}`** "
                f"({int(bulk_n.value)} points, {mode}). "
                f"Same entry point as `just bulk-cr`.{note_h2o}"
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
                extra = ""
                if result.meta.get("coord_OH_1.3A") is not None:
                    extra = (
                        f"  \nCF1 coords: OH≈{result.meta['coord_OH_1.3A']:.3g}, "
                        f"HH≈{result.meta.get('coord_HH_1.9A', float('nan')):.3g} "
                        f"(bare HNC under-associates vs full molecules)."
                    )
                oz_banner = mo.callout(
                    mo.md(
                        f"**{sym}**: {'converged' if ok else 'not converged'} "
                        f"· {result.n_iter} iters · residual `{result.residual:.3e}`  \n"
                        f"closure `{result.meta.get('closure')}` · "
                        f"|c_short| tail max `{tail:.3e}`{paths}{extra}"
                    ),
                    kind="success" if ok or tail < 1e-2 else "warn",
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
        _dedupe_legend(fig)
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
            "Workflow (wall analysis)": mo.md(
                r"""
                1. **TOML** (`data/fluid_parameters.toml`) → densities, charges, potential type.
                2. **Bulk** \(c_\mathrm{short}(r)\) artefact under `data/bulk/` (Picard OZ or CF1 anneal).
                3. **Wall setup** `prepare_fluid_run` interpolates bulk \(c\) onto the \(z\) grid and
                   builds \(f_1,f_2\) kernels and wall potential \(\phi_w\).
                4. **Solve** Newton–Krylov on \(F[t_w]=0\) (`solve_model` / `opt_func`).
                5. **Output** \(g_{wi}(z)=h_{wi}(z)+1\) with
                   \(h_{wi}=\mathrm{e}^{-\beta\phi_{wi}+t_{wi}}-1\).

                Streamlit uses the same core (`components/fluid_page.py`). Production demo:
                [phd-thesis-mjb.onrender.com](https://phd-thesis-mjb.onrender.com).
                """
            ),
            "Units & companions": mo.md(
                r"""
                - Length: Å · concentration: mol L⁻¹ · potentials → \(\beta u\)
                - Coulomb scale: `charge_pair / (2 r)` (see `charged_fluids.units`)
                - Bulk-only explorer: `just marimo-bulk`
                - CLI bulk regenerate: `just bulk-cr <symbol>`
                - Legacy Jupyter: `notebooks/legacy/`
                """
            ),
            "CF1 water notes": mo.md(
                r"""
                - Components ordered **O, H** (pairs OO, OH, HH).
                - Live bulk: Ng Coulomb split + HNC with temperature anneal (~6000 K → 298 K).
                - Digitised thesis Fig 3.12 bulk remains opt-in:
                  `[h2o.bulk] source = "digitised_fig3_12"`.
                - Bare HNC under-associates (coord_OH ≲ 2); Duh–Haymet bridge IE is future work.
                - Thesis Fig 3.2 wall benchmark (ψ₀ = 0) is still experimental in this stack.
                """
            ),
        }
    )
    footer
    return


if __name__ == "__main__":
    app.run()
