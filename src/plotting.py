import altair as alt
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st


def make_simple_plot(x, y, xlabel=None, ylabel=None, title=None, xliml=None, yliml=None, output="jupyter"):
    if xliml is None:
        xliml = [0, 10]
    if yliml is None:
        yliml = [-20, 40]
    fig, ax = plt.subplots()
    ax.plot(x, y)
    ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
    plt.xlim(xliml)
    plt.ylim(yliml)
    ax.grid()
    if output == "jupyter":
        plt.show()
    else:
        return fig


def fast_plot(x, y):
    return make_simple_plot(x=x, y=y)


def plotly_line(x, y, column_names, y_label="", legend_label="", title=None, xliml=None, yliml=None):
    df = pd.DataFrame(data=np.column_stack((x, y)), columns=column_names)
    fig = px.line(df, x=column_names[0], y=column_names[1:], title=title, 
        labels={column_names[0]: column_names[0],
        "value": y_label,
        "variable": legend_label})
    fig.update_xaxes(range=xliml)
    fig.update_yaxes(range=yliml)
    return fig



def plot_wall_curves(n_component, z, z_plots):
    for plot_title, details in z_plots.items():
        y_col_labels = ["z"]
        fn_label = details["fn_label"]
        y_col_labels.extend(f"{fn_label}_{i}" for i in range(n_component))
        try:
            xliml=details["xlim"]
        except Exception:
            xliml=None
        try:
            yliml=details["ylim"]
        except Exception:
            yliml=None

        fig = plotly_line(z, details["plot_fn"], y_col_labels, y_label=details["plot_name"], legend_label="",
                             xliml=xliml, yliml=yliml, title=plot_title)
        st.plotly_chart(fig)


def plot_bulk_curves(n_component, r, r_plots):
    for plot_title, details in r_plots.items():
        y_col_labels = ["r"]
        fn_label = details["fn_label"]
        for i in range(n_component):
            y_col_labels.extend(f"{fn_label}_{i}{j}" for j in range(i, n_component))
        try:
            xliml=details["xlim"]
        except Exception:
            xliml=None
        try:
            yliml=details["ylim"]
        except Exception:
            yliml=None

        fig = plotly_line(r, details["plot_fn"], y_col_labels, y_label=details["plot_name"], legend_label="",
                            xliml=xliml, yliml=yliml, title=plot_title)
        st.plotly_chart(fig)
