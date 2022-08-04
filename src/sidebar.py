import streamlit as st


def create_sidebar(fluid):
    # st.sidebar.text(fluid.name) - clear from menu bar
    st.sidebar.markdown("__Inputs:__")
    st.sidebar.markdown(f"Temperature (K): {fluid.temperature}")
    st.sidebar.markdown(f"Concentration $M / dm^3$: {fluid.concentration}")

    # st.sidebar.latex("\psi_0 (mV):")

    if fluid.symbol not in ["lj1", "lj2"]:
        psi_0 = st.sidebar.number_input(
            label=("\psi_0 (mV):"), value=0, min_value=-1000, max_value=1000)
    else:
        psi_0 = 0.0

    st.sidebar.markdown("__Numerical parameters:__")

    z_cutoff = st.sidebar.number_input(
        label="z_cutoff (A)", value=25, min_value=20, max_value=50)
    n_point = st.sidebar.number_input(label="nPoint:", value=2001)
    grid_size = z_cutoff / (n_point - 1)
    st.sidebar.markdown(f"Grid size (A): {grid_size}")

    st.sidebar.markdown("__Optimisation parameters:__")
    max_iteration = st.sidebar.number_input(
        "Maximum iterations", min_value=1, max_value=1000, value=1)
    tolerance = st.sidebar.number_input(
        "Convergence tolerance", min_value=1e-12, max_value=1e-6, value=1e-9)

    return z_cutoff, int(n_point), psi_0, tolerance, int(max_iteration)
