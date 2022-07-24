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

def plotly_line(x, y, x_header, y_header, title=None):
    df = pd.DataFrame(data=np.column_stack((x, y)), columns=["r", "u1", "u2", "u3"])
    df = pd.melt(df, id_vars='r', value_vars=df.columns[:-1])
    st.write(df)
    #fig = px.line(x=df.r, y=df.u1)
    return fig

# df = pd.DataFrame(data=np.column_stack((r, beta*u)), columns=["r", "u1", "u2", "u3"])