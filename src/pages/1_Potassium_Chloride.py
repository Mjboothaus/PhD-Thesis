
import streamlit as st
from helper_functions import read_render_markdown_file

def create_sidebar():
    st.sidebar.text("Inputs:")
    st.sidebar.number_input("Temperature (C):", value=25, min_value=-298, max_value=2000)
    st.sidebar.number_input(label=r"$\psi_0$ (mV):", value=0, min_value=-1000, max_value=1000)
    nPoint = st.sidebar.number_input("nPoint:", value=4096)
    grid_size = st.sidebar.number_input("Grid size (Angstrom):", value=0.006)
    z_cutoff = round((nPoint - 1) * grid_size, 4)
    st.sidebar.text(f"Cutoff (A): {z_cutoff}")

    return None

create_sidebar()