import matplotlib.pyplot as plt


def make_simple_plot(x, y, xlabel, ylabel, title, xliml = None, yliml = None):
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
    return make_simple_plot(x=x, y=y, xlabel="", ylabel="", title="Fast plot")