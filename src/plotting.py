import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st


def plotly_line(x, y, column_names, y_label="", legend_label="", title=None, xliml=None, yliml=None):
    df = pd.DataFrame(data=np.column_stack((x, y)), columns=column_names)
    fig = px.line(df, x=column_names[0], y=column_names[1:], title=title,
                  labels={column_names[0]: column_names[0],
                          "value": y_label,
                          "variable": legend_label})
    fig.update_xaxes(range=xliml)
    fig.update_yaxes(range=yliml)
    return fig


def plot_wall_curves(n_component, z, z_plots, component):
    for plot_title, details in z_plots.items():
        y_col_labels = ["z"]
        fn_label = details["fn_label"]
        y_col_labels.extend(f"{fn_label}_w{component[i]}" for i in range(n_component))
        try:
            xliml = details["xlim"]
        except Exception:
            xliml = None
        try:
            yliml = details["ylim"]
        except Exception:
            yliml = None

        fig = plotly_line(z, details["plot_fn"], y_col_labels, y_label=details["plot_name"], legend_label="",
                          xliml=xliml, yliml=yliml, title=plot_title)
        #show(fig) - Jupyter
        st.plotly_chart(fig)


def plot_bulk_curves(n_component, r, r_plots, component):
    for plot_title, details in r_plots.items():
        y_col_labels = ["r"]
        fn_label = details["fn_label"]
        for i in range(n_component):
            y_col_labels.extend(
                f"{fn_label}_{component[i]}{component[j]}" for j in range(i, n_component))
        xliml = details["xlim"] if "xlim" in details.keys() else None
        yliml = details["ylim"] if "ylim" in details.keys() else None
 
        fig = plotly_line(r, details["plot_fn"], y_col_labels, y_label=details["plot_name"], legend_label="",
                          xliml=xliml, yliml=yliml, title=plot_title)
        #fig.show() - Jupyter
        st.plotly_chart(fig)


    
def plot_convergence_history(history, log_y=True):
    """Plot :class:`charged_fluids.ConvergenceHistory` from ``solve_model``."""
    if history is None or not getattr(history, "iterations", None):
        st.info("No convergence history recorded.")
        return None
    iters = np.asarray(history.iterations, dtype=float)
    y = np.column_stack(
        [
            np.asarray(history.residuals, dtype=float),
            np.asarray(history.max_abs, dtype=float),
            np.asarray(history.mean_abs, dtype=float),
        ]
    )
    fig = plotly_line(
        iters,
        y,
        ["Iteration #", "|F|_2", "max |F|", "mean |F|"],
        y_label="residual",
        title="Newton–Krylov residual history",
    )
    if log_y:
        fig.update_yaxes(type="log")
    st.plotly_chart(fig, use_container_width=True)
    return fig


def plot_convergence(filename=None, log_y=True, history=None):
    """Plot solver residual vs iteration.

    Prefer structured ``history`` from ``solve_model``; fall back to parsing a
    SciPy ``disp=True`` log file (legacy path).
    """
    if history is not None and getattr(history, "iterations", None):
        return plot_convergence_history(history, log_y=log_y)
    if not filename:
        st.info("No convergence log or history to plot.")
        return None
    with open(filename) as solver_out:
        content = solver_out.readlines()
    conv_table = list(
        map(
            np.float32,
            [
                line.replace(":", ",")
                .replace("|F(x)| = ", "")
                .split(";")[0]
                .split(",")
                for line in content
                if line != "\n"
            ],
        )
    )
    if not conv_table:
        st.info("Convergence log was empty.")
        return None
    convergence_df = pd.DataFrame(
        conv_table, columns=["Iteration #", "|F(x)|"], dtype=np.float32
    )
    fig = plotly_line(
        x=convergence_df["Iteration #"],
        y=convergence_df["|F(x)|"],
        column_names=["Iteration #", "|F(x)|"],
    )
    if log_y:
        fig.update_yaxes(type="log")
    st.plotly_chart(fig, use_container_width=True)
    return fig
