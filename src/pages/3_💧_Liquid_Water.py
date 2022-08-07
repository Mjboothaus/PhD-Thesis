import streamlit as st
from parameters import set_fluid_parameters
from sidebar import create_sidebar

# from helper_functions import read_render_markdown_file

fluid_symbol = "h2o"

fluid = set_fluid_parameters(fluid_symbol)

if fluid is not None:
    create_sidebar(fluid)
    st.write(fluid)
else:
    st.error("Invalid choice of fluid")

st.markdown("Liquid water at a charged interface: Coming soon...")
