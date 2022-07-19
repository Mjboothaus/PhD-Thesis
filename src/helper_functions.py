from IPython.display import Markdown, display
from streamlit import markdown
from pathlib import Path


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
