from pathlib import Path

from IPython.display import Markdown, display
from streamlit import cache_data, markdown


@cache_data
def st_read_markdown_file(markdown_file):
    return Path(markdown_file).read_text()


def read_markdown_file(markdown_file):
    return Path(markdown_file).read_text()


def read_render_markdown_file(markdown_file, output="jupyter"):
    if output == "jupyter":
        try:
            md_text = read_markdown_file(markdown_file)
            display(Markdown(md_text))
        except Exception:
            print(f"Error with markdown file: {markdown_file}")
            return None
    else:
        try:
            md_text = st_read_markdown_file(markdown_file)
            markdown(md_text, unsafe_allow_html=True)
        except Exception:
            print(f"Error with markdown file: {markdown_file}")
            return None
