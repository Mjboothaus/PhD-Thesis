import matplotlib.pyplot as plt


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