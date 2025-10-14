from pathlib import Path
from streamlit import markdown, cache_data
from typing import Optional

@cache_data
def read_markdown_file(markdown_file: str) -> Optional[str]:
    """Read a markdown file and return its content."""
    try:
        return Path(markdown_file).read_text()
    except Exception as e:
        print(f"Error reading markdown file {markdown_file}: {e}")
        return None

def read_render_markdown_file(markdown_file: str, output: str = "streamlit") -> None:
    """Read and render a markdown file.
    
    Args:
        markdown_file: Path to the markdown file
        output: Output format ('streamlit' or 'jupyter')
    """
    md_text = read_markdown_file(markdown_file)
    if md_text is None:
        return
        
    if output == "jupyter":
        try:
            # Only import IPython if needed for Jupyter output
            from IPython.display import Markdown, display
            display(Markdown(md_text))
        except ImportError:
            print("IPython not available. Install with 'uv pip install ipython' for Jupyter support.")
            print(md_text)  # Fallback to plain text
        except Exception as e:
            print(f"Error displaying markdown: {e}")
    else:
        try:
            markdown(md_text, unsafe_allow_html=True)
        except Exception as e:
            print(f"Error displaying markdown in Streamlit: {e}")

