from pathlib import Path

import matplotlib.pyplot as plt
from IPython.display import Markdown, display
from streamlit import markdown

def make_simple_plot(x, y, xlabel, ylabel, title):
    _, ax = plt.subplots()
    ax.plot(x, y)
    ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
    plt.xlim([0, 10])
    plt.ylim([-20, 40])
    ax.grid()
    plt.show()
    return title

def read_render_markdown_file(markdown_file, output="jupyter"):
    try:
        md_text = Path(markdown_file).read_text()
        if output == "jupyter":
            display(Markdown(md_text))
        else:
            markdown(md_text, unsafe_allow_html=False)
    except Exception:
        print(f"Error with markdown file: {markdown_file}")
        return None