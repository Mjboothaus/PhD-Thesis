from pathlib import Path

import streamlit as st
from scipy.constants import elementary_charge, epsilon_0

from helper_functions import read_render_markdown_file

APP_TITLE = "Singlet integral equation for liquids"

st.set_page_config(
    page_title=APP_TITLE,
    layout="wide",
    menu_items={
        "About":
        "Created with love & care at DataBooth - www.databooth.com.au"
    })

# Other constants - Molten KCl (page 65)

# Temperature
# Concentration
# Bulk direct correlation functions
# Short range wall potential
# Surface potential or charge

# Solving for t(z) = ln g + b \phi(z)
# Newton-GMRES algorithm


def create_app_header(app_title, subtitle=None):
    st.header(app_title)
    if subtitle is not None:
        st.subheader(subtitle)
    return None


create_app_header(APP_TITLE, "Michael J. Booth - Ph.D (Science)")

# "Charged fluids near interfaces: Integral equation theory"

read_render_markdown_file("docs/app_main.md", output="streamlit")
