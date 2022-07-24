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
    #df = pd.DataFrame(data=np.c_[x, y], columns=column_names)
    fig = px.line(df, x=column_names[0], y=column_names[1:], title=title, 
        labels={column_names[0]: column_names[0],
        "value": y_label,
        "variable": legend_label})
    fig.update_xaxes(range=xliml)
    fig.update_yaxes(range=yliml)
    return fig
