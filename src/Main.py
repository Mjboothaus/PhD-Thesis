"""Streamlit multipage entry — charged fluids near interfaces."""

from __future__ import annotations

import sys
from pathlib import Path

import streamlit as st

# Ensure `src/` is on sys.path so multipage imports resolve locally and in Docker.
_SRC = Path(__file__).resolve().parent
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from helper_functions import read_render_markdown_file

APP_TITLE = "Charged fluids near interfaces"
SUB_TITLE = "Integral equation theory"

st.set_page_config(
    page_title=APP_TITLE,
    layout="wide",
    menu_items={
        "About": "Created with love & care at DataBooth — www.databooth.com.au"
    },
)


def create_app_header(app_title: str, subtitle: str | None = None) -> None:
    st.header(app_title)
    if subtitle is not None:
        st.subheader(subtitle)


create_app_header(APP_TITLE, SUB_TITLE)
read_render_markdown_file("docs/app_main.md", output="streamlit")
