from pathlib import Path

import streamlit as st

from helper_functions import read_render_markdown_file

APP_TITLE = "Charged fluids near interfaces"
SUB_TITLE = "Integral Equation Theory"

st.set_page_config(
    page_title=APP_TITLE,
    layout="wide",
    menu_items={
        "About":
        "Created with love & care at DataBooth - www.databooth.com.au"
    })

def create_app_header(app_title, subtitle=None):
    st.header(app_title)
    if subtitle is not None:
        st.subheader(subtitle)
    return None

create_app_header(APP_TITLE, SUB_TITLE)

read_render_markdown_file("docs/app_main_inc.md", output="streamlit")

# Following was an attempt to get LaTeX working in Plotly labels - unsuccessful

# mathjax = """
# <script type="text/javascript" async
#   src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.7/MathJax.js?config=TeX-MML-AM_CHTML">
# </script>
# """

# st.markdown(mathjax, unsafe_allow_html=True)