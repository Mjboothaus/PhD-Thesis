import streamlit as st

def create_sidebar(fluid):
    st.sidebar.text(fluid.name)
    st.sidebar.text("Inputs:")
    st.sidebar.number_input(label=r"$\psi_0$ (mV):",
                            value=0, min_value=-1000, max_value=1000, key=fluid.index)
    z_cutoff = st.sidebar.number_input(
        "z_cutoff (A)", value=25, min_value=20, max_value=50, key=fluid.index)
    nPoint = st.sidebar.number_input("nPoint:", value=2001, key=fluid.index)
    grid_size = z_cutoff / (nPoint - 1)
    st.sidebar.text(f"Grid size (A): {grid_size}")
    return None