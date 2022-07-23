import streamlit as st
from model_parameters import Fluid, set_fluid_parameters
from numerical_parameters import create_sidebar

# from helper_functions import read_render_markdown_file

fluid_symbol = "kcl"

fluid = set_fluid_parameters(fluid_symbol)

if fluid is not None:
    create_sidebar(fluid)
    st.write(fluid)
else:
    st.error("Invalid choice of fluid")

st.markdown("Molten potassium chloride at a charged interface: Now in progress...")