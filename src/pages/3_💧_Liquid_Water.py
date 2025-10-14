from pathlib import Path

import numpy as np
import pandas as pd
import streamlit as st
from modelling import *
from numerics import set_num_parameters
from parameters import fluid_specific_parameters, set_fluid_parameters
from plotting import plot_bulk_curves, plot_wall_curves
from sidebar import create_sidebar
from components.fluid_calculation import create_calculation_tab

# Initialize fluid parameters
fluid_symbol = "h2o"
fluid = set_fluid_parameters(fluid_symbol)

if fluid is not None:
    z_cutoff, n_point, psi_0, tolerance, max_iteration = create_sidebar(fluid)
    
    # Set up numerical parameters
    d = set_num_parameters(n_point, z_cutoff, fluid.n_component,
                          fluid.n_pair, tolerance, max_iteration)
    
    n_component = d.n_component
    n_pair = d.n_pair
    
    z = np.linspace(0.0, z_cutoff, n_point)
    z_index = np.arange(0, n_point, dtype=int)
    wall_zeros = np.zeros((n_point, n_component))
    fluid_zeros = np.zeros((n_point, n_pair))
    
    # Model setup
    model = Model(z=z, z_index=z_index, hw=wall_zeros,
                 c_short=fluid_zeros, f1=fluid_zeros, f2=fluid_zeros)
    
    # Basic parameters
    beta_phiw = np.zeros_like(wall_zeros)  # Placeholder
    beta_psi_charge = 0.0  # Placeholder
    
    st.subheader(f"Charged fluids near an interface: {fluid.name}")
    
    tab1, tab2, tab0 = st.tabs(["Calculation", "Output graphs", "Bulk properties"])
    
    with tab1:
    solution, hw_solution, calculation_run = create_calculation_tab(
        fluid=fluid,
        model=model,
        d=d,
        beta_phiw=beta_phiw,
        beta_psi_charge=beta_psi_charge,
        n_component=n_component,
        n_pair=n_pair,
        z=z,
        status="in_progress"
    )
    
    with tab2:
        st.info("Output graphs will be available once the model implementation is complete.")
    
    with tab0:
        st.markdown("#")
        col1, col2, col3 = st.columns(3)
        col1.metric("Temperature (K)", f"{fluid.temperature}")
        col2.metric("Concentration (M/dm3)", f"{fluid.concentration[0]}")
        col3.metric("Components", f"{', '.join(fluid.component)}")
        
        st.info("Bulk properties will be available once the model implementation is complete.")
else:
    st.error("Invalid choice of fluid")

