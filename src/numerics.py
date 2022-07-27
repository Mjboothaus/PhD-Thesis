from dataclasses import dataclass

import numpy as np
import streamlit as st


@dataclass
class Discretisation:
    n_point: int
    n_component: int
    n_pair: int
    z_cutoff: float
    grid_size: float
    tolerance: float
    max_iteration: int

# sidebar doesn't strictly belong with numerics... might be convenient though


def create_sidebar(fluid):
    # st.sidebar.text(fluid.name) - clear from menu bar
    st.sidebar.text(f"Temperature (K): {fluid.temperature}")
    st.sidebar.text(f"Concentration (m/dm3): {fluid.concentration[0]}")

    st.sidebar.text("Inputs:")
    psi_0 = st.sidebar.number_input(label=r"$\psi_0$ (mV):", value=0, min_value=-1000, max_value=1000)

    st.sidebar.markdown("Numerical parameters:")

    z_cutoff = st.sidebar.number_input("z_cutoff (A)", value=25, min_value=20, max_value=50)
    n_point = st.sidebar.number_input("nPoint:", value=2001)
    grid_size = z_cutoff / (n_point - 1)
    st.sidebar.text(f"Grid size (A): {grid_size}")

    st.sidebar.markdown("Optimisation parameters:")
    max_iteration = st.sidebar.number_input("Maximum iterations", min_value=10, max_value=1000, value=20)
    tolerance = st.sidebar.number_input("Convergence tolerance", min_value=1e-12, max_value=1e-6, value=1e-9)

    return z_cutoff, int(n_point), psi_0, tolerance, int(max_iteration)


def set_num_parameters(n_point, z_cutoff, n_component, n_pair, tolerance, max_iteration):
    grid_size = z_cutoff / (n_point - 1)
    return Discretisation(n_point, n_component, n_pair, z_cutoff, grid_size, tolerance, max_iteration)
