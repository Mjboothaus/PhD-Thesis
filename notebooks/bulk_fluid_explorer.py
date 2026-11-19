"""Marimo notebook: bulk c(r) artefacts, fluid parameters, and optional OZ.

Consolidates the useful parts of the legacy Jupyter notebooks
(`bulk_fluid_explore`, `fluid_parameters`, `lj`, `check_units` ideas) against
the modern `charged_fluids` package.

Run:
  just marimo-bulk
  # or
  uv run marimo edit notebooks/bulk_fluid_explorer.py
"""

import marimo

__generated_with = "0.25.2"
app = marimo.App(width="medium")


with app.setup:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots


def _add_pair_series(fig, x, y, *, row=None, col=None, labels=None):
    import numpy as np

    y = np.asarray(y)
    if y.ndim == 1:
        series = [y]
        names = labels or ["series 0"]
    elif y.shape[0] == len(x):
        series = [y[:, i] for i in range(y.shape[1])]
        names = labels or [f"pair {i}" for i in range(y.shape[1])]
    else:
        series = [y[i] for i in range(y.shape[0])]
        names = labels or [f"pair {i}" for i in range(y.shape[0])]
    for name, yi in zip(names, series, strict=False):
        tr = go.Scatter(x=x, y=yi, mode="lines", name=name)
        if row is None:
            fig.add_trace(tr)
        else:
            fig.add_trace(tr, row=row, col=col)


def _pair_labels(components: list[str]) -> list[str]:
    labels: list[str] = []
    for i, a in enumerate(components):
        for b in components[i:]:
            labels.append(f"{a}-{b}")
    return labels


@app.cell
def _intro():
    import marimo as mo

    intro = mo.md(
        r"""
        # Bulk fluid explorer

        Inspect **stored** bulk \(c_\mathrm{short}(r)\) artefacts and TOML fluid
        parameters. Optionally re-run Picard OZ (dry-run by default).

        For wall–fluid singlet solves see `notebooks/wall_fluid_explorer.py`
        (`just marimo`).

        Legacy Jupyter notebooks under `notebooks/legacy/` are archival only.
        """
    )
    intro
    return (mo,)


@app.cell
def _controls(mo):
    from charged_fluids import list_fluid_symbols

    all_syms = list(list_fluid_symbols())
    symbol = mo.ui.dropdown(
        options=all_syms,
        value="lj1" if "lj1" in all_syms else all_syms[0],
        label="Fluid",
    )
    r_max = mo.ui.slider(5.0, 40.0, value=12.0, step=1.0, label="plot r_max (A)")
    controls = mo.vstack(
        [
            mo.md("### Artefact & parameters"),
            mo.hstack([symbol, r_max], justify="start"),
        ]
    )
    controls
    return r_max, symbol


@app.cell
def _params_table(mo, symbol):
    from charged_fluids import load_fluid_table, set_fluid_parameters

    try:
        fluid = set_fluid_parameters(symbol.value)
        if fluid is None:
            raise KeyError(symbol.value)
        table = load_fluid_table()
        raw = table.get(symbol.value, {})
        bulk = raw.get("bulk", {}) if isinstance(raw, dict) else {}
        rows = [
            ("name", fluid.name),
            ("symbol", fluid.symbol),
            ("status", fluid.status),
            ("components", ", ".join(fluid.component)),
            ("n_component / n_pair", f"{fluid.n_component} / {fluid.n_pair}"),
            ("temperature (K)", fluid.temperature),
            ("concentration (M)", fluid.concentration.tolist()),
            ("epsilon_r", fluid.epsilon_r),
            ("cr_filename", fluid.cr_filename or "—"),
            ("potential.type", fluid.potential.get("type", "—")),
        ]
        if isinstance(bulk, dict):
            for k, v in bulk.items():
                rows.append((f"bulk.{k}", v))
        body = "\n".join(f"| `{k}` | {v} |" for k, v in rows)
        params_md = mo.md(
            f"""
            ### Fluid parameters (`{symbol.value}`)

            | Key | Value |
            |---|---|
            {body}
            """
        )
    except Exception as exc:  # noqa: BLE001
        params_md = mo.md(f"**Parameter load failed:** `{type(exc).__name__}: {exc}`")
        fluid = None
    params_md
    return (fluid,)


@app.cell
def _load_artefact(mo, symbol):
    import json
    from pathlib import Path

    import numpy as np

    from charged_fluids.bulk.artefacts import default_bulk_dir, read_meta

    cr = None
    meta = None
    bulk_dir = default_bulk_dir()
    path = bulk_dir / f"{symbol.value}.cr.dat"
    meta_path = bulk_dir / f"{symbol.value}.meta.json"
    try:
        if path.is_file():
            data = np.loadtxt(path)
            if data.ndim == 1:
                cr = {"r": data, "c_short": np.zeros_like(data)}
            else:
                cr = {"r": data[:, 0], "c_short": data[:, 1:]}
        if meta_path.is_file():
            meta = read_meta(meta_path)
        msg = mo.md(
            f"Artefact `{path}` · "
            f"{'found' if path.is_file() else '**missing**'} · "
            f"meta {'found' if meta_path.is_file() else 'missing'}"
        )
    except Exception as exc:  # noqa: BLE001
        msg = mo.md(f"**Artefact load failed:** `{type(exc).__name__}: {exc}`")
        cr = None
        meta = None
    msg
    return cr, meta


@app.cell
def _meta_view(meta, mo):
    if not meta:
        out = mo.md("_No `.meta.json` for this fluid._")
    else:
        import json

        pretty = json.dumps(meta, indent=2, default=str)
        out = mo.md(f"### Artefact metadata\n\n```json\n{pretty}\n```")
    out
    return


@app.cell
def _cr_plot(cr, fluid, mo, r_max):
    if cr is None:
        out = mo.md("_No c(r) data to plot._")
    else:
        r = cr["r"]
        c = cr["c_short"]
        labels = None
        if fluid is not None:
            labels = _pair_labels(list(fluid.component))
            n_series = c.shape[1] if c.ndim == 2 and c.shape[0] == len(r) else (
                c.shape[0] if c.ndim == 2 else 1
            )
            if len(labels) != n_series:
                labels = None
        fig = go.Figure()
        _add_pair_series(fig, r, c, labels=labels)
        fig.update_layout(
            title=f"Stored c_short(r) — {getattr(fluid, 'symbol', '?')}",
            xaxis_title="r (A)",
            yaxis_title="c_short",
            xaxis_range=[0, float(r_max.value)],
            height=420,
            margin=dict(t=50, b=40),
        )
        out = fig
    out
    return


@app.cell
def _units_note(mo):
    note = mo.md(
        r"""
        ---
        ### Units (modern stack)

        See `charged_fluids.units` — not the old `numericalunits` playground.

        | Quantity | Convention |
        |---|---|
        | Length | Ångström on wall/bulk grids |
        | Temperature | kelvin |
        | Concentration | mol L⁻¹ |
        | Coulomb scale | `charge_pair / (2 r)` in the thesis reduced form |
        | Potentials | joule (then \(\beta u\)) |

        Legacy notebooks that experimented with `numericalunits` live under
        `notebooks/legacy/check_units.ipynb`.
        """
    )
    note
    return


@app.cell
def _oz_controls(mo, symbol):
    n_points = mo.ui.slider(256, 2048, value=512, step=256, label="n_points")
    dry_run = mo.ui.checkbox(value=True, label="Dry-run (do not write)")
    run_oz = mo.ui.run_button(label="Run bulk OZ")
    box = mo.vstack(
        [
            mo.md(
                r"""
                ---
                ## Optional bulk OZ regenerate

                Same entry point as `just bulk-cr` /
                `python -m charged_fluids.bulk.cli`.
                """
            ),
            mo.hstack([n_points, dry_run, run_oz], justify="start"),
            mo.md(f"Target fluid: **`{symbol.value}`**"),
        ]
    )
    box
    return dry_run, n_points, run_oz


@app.cell
def _oz_run(dry_run, mo, n_points, run_oz, symbol):
    from charged_fluids.bulk.generate import generate_bulk_cr

    result = None
    if not run_oz.value:
        msg = mo.md("_Press **Run bulk OZ** when ready._")
    else:
        with mo.status.spinner(f"Bulk OZ {symbol.value}…"):
            try:
                result, cr_path, meta_path = generate_bulk_cr(
                    symbol.value,
                    n_points=int(n_points.value),
                    write=not bool(dry_run.value),
                )
            except Exception as exc:  # noqa: BLE001
                msg = mo.md(f"**Failed:** `{type(exc).__name__}: {exc}`")
                result = None
            else:
                status = "converged" if result.converged else "NOT converged"
                paths = ""
                if cr_path is not None:
                    paths = f"  \nWrote `{cr_path}`, `{meta_path}`"
                msg = mo.md(
                    f"**{symbol.value}**: {status} · iters={result.n_iter} · "
                    f"residual=`{result.residual:.3e}` · "
                    f"closure=`{result.meta.get('closure')}`{paths}"
                )
    msg
    return (result,)


@app.cell
def _oz_plot(mo, r_max, result):
    if result is None:
        out = mo.md("_OZ plots appear after a run._")
    else:
        fig = make_subplots(
            rows=1,
            cols=2,
            subplot_titles=("c_short(r)", "g(r)"),
            horizontal_spacing=0.08,
        )
        _add_pair_series(fig, result.r, result.c_short, row=1, col=1)
        _add_pair_series(fig, result.r, result.g, row=1, col=2)
        fig.add_hline(y=1.0, line_dash="dash", line_color="grey", line_width=1, row=1, col=2)
        xmax = min(float(r_max.value), float(result.r[-1]))
        fig.update_xaxes(range=[0, xmax], title_text="r (A)")
        fig.update_yaxes(title_text="c_short", row=1, col=1)
        fig.update_yaxes(title_text="g(r)", row=1, col=2)
        fig.update_layout(height=400, margin=dict(t=50, b=40))
        out = fig
    out
    return


if __name__ == "__main__":
    app.run()
