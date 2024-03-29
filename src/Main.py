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

read_render_markdown_file("docs/app_main.md", output="streamlit")
