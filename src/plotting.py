import matplotlib.pyplot as plt
import altair as alt
import pandas as pd
import numpy as np
import streamlit as st


def make_simple_plot(x, y, xlabel=None, ylabel=None, title=None, xliml=None, yliml=None):
    if xliml is None:
        xliml = [0, 10]
    if yliml is None:
        yliml = [-20, 40]
    _, ax = plt.subplots()
    ax.plot(x, y)
    ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
    plt.xlim(xliml)
    plt.ylim(yliml)
    ax.grid()
    plt.show()
    return title


def fast_plot(x, y):
    return make_simple_plot(x=x, y=y)


import plotly.express as px

def plot_plotly_line(x, y, column_names):
    df = pd.DataFrame(data=np.column_stack((x, y)), columns=column_names)
    #df = pd.DataFrame(data=np.c_[x, y], columns=column_names)
    fig = px.line(df, x=column_names[0], y=column_names[1:])
    return fig