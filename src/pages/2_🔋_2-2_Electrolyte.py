import streamlit as st
from parameters import Fluid, set_fluid_parameters
from sidebar import create_sidebar

# from helper_functions import read_render_markdown_file

fluid_symbol = "2_2"

fluid = set_fluid_parameters(fluid_symbol)

if fluid is not None:
    create_sidebar(fluid)
    st.write(fluid)
else:
    st.error("Invalid choice of fluid")

st.markdown("2-2 Aqueous electrolyte at a charged interface: Coming soon...")