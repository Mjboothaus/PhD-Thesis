import matplotlib.pyplot as plt


def make_simple_plot(x, y, xlabel, ylabel, title, xliml=[0, 10], yliml=[-20, 40]):
    _, ax = plt.subplots()
    ax.plot(x, y)
    ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
    plt.xlim(xliml)
    plt.ylim(yliml)
    ax.grid()
    plt.show()
    return title