
import streamlit as st
from helper_functions import read_render_markdown_file

def create_sidebar():
    st.sidebar.text("Inputs:")
    st.sidebar.number_input(label=r"$\psi_0$ (mV):", value=0, min_value=-1000, max_value=1000)
    z_cutoff = st.sidebar.number_input("z_cutoff (A)", value=25, min_value=20, max_value=50)
    nPoint = st.sidebar.number_input("nPoint:", value=2001)
    grid_size = z_cutoff / (nPoint - 1)
    st.sidebar.text(f"Grid size (A): {grid_size}")
    return None

create_sidebar()