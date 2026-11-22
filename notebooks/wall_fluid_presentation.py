# /// script
# requires-python = ">=3.10,<3.15"
#
# [tool.marimo-studio]
# default = "explainer"
# ///

import marimo

__generated_with = "0.25.2"
app = marimo.App(width="medium", app_title="Wall–fluid presentation")


with app.setup:
    from pathlib import Path

    import numpy as np
    import plotly.graph_objects as go
    import plotly.io as pio
    from plotly.subplots import make_subplots

    pio.templates.default = "plotly_white"

    REPO = Path(__file__).resolve().parents[1]
    DATA = REPO / "data"


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell
def title(mo):
    title_md = mo.md(
        r"""
        # Charged fluids near interfaces
        ### Integral equation theory — wall–fluid presentation

        Modern Python reimplementation of a 1998 PhD.  
        Production app: [phd-thesis-mjb.onrender.com](https://phd-thesis-mjb.onrender.com)
        """
    )
    title_md
    return (title_md,)


@app.cell
def agenda(mo):
    agenda_md = mo.md(
        r"""
        ## Agenda

        1. **Problem** — fluids at a charged wall  
        2. **Method** — bulk OZ → wall singlet residual  
        3. **Stack** — `charged_fluids`, Streamlit, Marimo  
        4. **Fluids** — LJ, KCl, CF1 water, electrolytes  
        5. **Thesis Fig 3.2** — digitised wall–O / wall–H benchmark  
        6. **Open work** — CF1 wall residual vs Fig 3.2  
        """
    )
    agenda_md
    return (agenda_md,)


@app.cell
def problem(mo):
    problem_md = mo.md(
        r"""
        ## The physical problem

        How do ions and solvent structure next to a **charged surface**?

        - Singlet wall–fluid distribution functions \(g_{wi}(z)\)
        - Input: bulk short-range direct correlations \(c_\mathrm{short}(r)\)
        - Closure + Ornstein–Zernike theory (thesis: integral equation route)

        Thesis: *Charged fluids near interfaces* (Booth, 1998).
        """
    )
    problem_md
    return (problem_md,)


@app.cell
def method(mo):
    method_md = mo.md(
        r"""
        ## Method (one line)

        ```text
        TOML parameters
          → bulk Picard OZ (HNC/PY, Ng, optional T anneal)
          → data/bulk/{symbol}.cr.dat
          → prepare_fluid_run (c_short, f₁, f₂, wall φ)
          → Newton–Krylov  F[t_w] = 0
          → g_wi(z) = h_wi(z) + 1
        ```

        Coulomb scale is shared: \(\beta u_c = \mathrm{charge\_pair}/(2r)\).
        """
    )
    method_md
    return (method_md,)


@app.cell
def stack(mo):
    stack_md = mo.md(
        r"""
        ## Software stack

        | Layer | Role |
        |-------|------|
        | **`charged_fluids`** | UI-free science core |
        | **Streamlit** | Production multipage app (Render) |
        | **Marimo** | Research explorers + this presentation |
        | **uv + Just** | Deps and tasks |

        Interactive solves: `just marimo` → `wall_fluid_explorer.py`.  
        This notebook is **presentation-only** (Studio pilot).
        """
    )
    stack_md
    return (stack_md,)


@app.cell
def fluids_table(mo):
    fluids_md = mo.md(
        r"""
        ## Fluids in the modern stack

        | Symbol | Status | Bulk |
        |--------|--------|------|
        | **lj1** / **lj2** | ready | HNC LJ |
        | **kcl** | experimental | Ng/PY + T anneal |
        | **h2o** (CF1) | experimental | Ng/HNC + T anneal |
        | **1_1** / **2_2** | experimental | Ng + anneal, soft RPM |

        Detail: `docs/PHYSICS.md` · plans: `docs/ROADMAP.md`.
        """
    )
    fluids_md
    return (fluids_md,)


@app.cell
def load_fig32():
    from charged_fluids import load_fig32_profiles

    fig32 = load_fig32_profiles(DATA)
    z32 = fig32["z"]
    g_wO = fig32["g_wO"]
    g_wH = fig32["g_wH"]
    return fig32, g_wH, g_wO, z32


@app.cell
def fig32_intro(mo):
    fig32_intro_md = mo.md(
        r"""
        ## Thesis benchmark — Fig 3.2

        Case 1, \(\psi_0 = 0\) mV: normalised wall density profiles for **oxygen**
        and **hydrogen** sites (CF1 water).

        Digitised from `data/Fig3.2.jpg` into CSVs; peaks ~**7** (O) and ~**8.5** (H)
        near \(z \approx 1\) Å. Bulk \(c(r)\) is held fixed for all CF1 wall cases.
        """
    )
    fig32_intro_md
    return (fig32_intro_md,)


@app.cell
def fig32_plot(g_wH, g_wO, z32):
    _fig = make_subplots(
        rows=2,
        cols=1,
        shared_xaxes=True,
        subplot_titles=(r"Wall–oxygen \(g_{wO}(z)\)", r"Wall–hydrogen \(g_{wH}(z)\)"),
        vertical_spacing=0.12,
    )
    _fig.add_trace(
        go.Scatter(
            x=z32,
            y=g_wO,
            mode="lines",
            name="thesis g_wO",
            line=dict(color="#111", width=2),
        ),
        row=1,
        col=1,
    )
    _fig.add_trace(
        go.Scatter(
            x=z32,
            y=g_wH,
            mode="lines",
            name="thesis g_wH",
            line=dict(color="#111", width=2),
        ),
        row=2,
        col=1,
    )
    _fig.add_hline(y=1.0, line_dash="dot", line_color="#999", row=1, col=1)
    _fig.add_hline(y=1.0, line_dash="dot", line_color="#999", row=2, col=1)
    _fig.update_xaxes(title_text="z (Å)", range=[0, 15], row=2, col=1)
    _fig.update_yaxes(title_text="g_wO", range=[-0.5, 9], row=1, col=1)
    _fig.update_yaxes(title_text="g_wH", range=[-0.5, 9], row=2, col=1)
    _fig.update_layout(
        height=520,
        showlegend=False,
        title=dict(text="Digitised thesis Fig 3.2 · ψ₀ = 0 mV", x=0.02, xanchor="left"),
        margin=dict(l=56, r=24, t=56, b=48),
    )
    fig32_figure = _fig
    fig32_figure
    return (fig32_figure,)


@app.cell
def fig32_stats(g_wH, g_wO, mo, z32):
    i_o = int(np.argmax(g_wO))
    i_h = int(np.argmax(g_wH))
    peak_o = float(g_wO[i_o])
    peak_h = float(g_wH[i_h])
    z_o = float(z32[i_o])
    z_h = float(z32[i_h])
    ginf_o = float(g_wO[-1])
    ginf_h = float(g_wH[-1])
    stats_md = mo.md(
        f"""
        ### Digitised peak summary

        | Species | Peak g | z (Å) | g at 15 Å |
        |---------|-------:|------:|----------:|
        | O | {peak_o:.2f} | {z_o:.2f} | {ginf_o:.2f} |
        | H | {peak_h:.2f} | {z_h:.2f} | {ginf_h:.2f} |

        API: `charged_fluids.load_fig32_profiles` / `compare_to_fig32`.
        """
    )
    stats_md
    return (stats_md,)


@app.cell
def open_work(mo):
    open_md = mo.md(
        r"""
        ## Open work (parked roadmap)

        **P0** — Live CF1 wall Newton–Krylov does **not** yet recover Fig 3.2
        contact peaks (residual blow-up / flat \(g\)). Digitised curves are the
        regression target.

        **P1** — CF1 bulk bridge (Duh–Haymet) so coord_OH → ~2.  
        **P2** — Pauling ionic bulk core; electrolyte thermodynamics.  
        **P3** — LJ table benchmarks; shim/data hygiene.

        Full plan: `docs/ROADMAP.md`.
        """
    )
    open_md
    return (open_md,)


@app.cell
def try_explorer(mo):
    try_md = mo.md(
        r"""
        ## Try it yourself

        ```bash
        just marimo              # interactive wall explorer
        just marimo-bulk         # bulk c(r) explorer
        just marimo-presentation # this notebook
        just app                 # Streamlit production UI
        ```

        Marimo Studio views (pilot): explainer + lecture under
        `studio/wall_fluid_presentation/` (or `notebooks/studio/...`).
        """
    )
    try_md
    return (try_md,)


@app.cell
def closing(mo):
    closing_md = mo.md(
        r"""
        ## Thanks

        Science core: **`charged_fluids`**.  
        UI: Streamlit (deploy) · Marimo (explore + present).  
        Production: https://phd-thesis-mjb.onrender.com
        """
    )
    closing_md
    return (closing_md,)


if __name__ == "__main__":
    app.run()
