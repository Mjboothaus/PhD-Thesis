"""2-2 Electrolyte fluid near a charged interface."""
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
fluid_symbol = "2_2"
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
    
    # Calculate fluid parameters
    fluid.beta = calc_beta(fluid.temperature)
    fluid.epsilon = calc_epsilon(fluid.epsilon_r)
    fluid.charge = calc_charge(fluid.valence)
    fluid.charge_pair = calc_charge_pair(
        fluid.beta, fluid.charge, fluid.epsilon, n_component, n_pair)
    fluid.rho = calc_rho(fluid.concentration)
    
    kappa = calc_kappa(fluid.beta, fluid.charge, fluid.rho, fluid.epsilon)
    st.sidebar.markdown(f"__Inverse Debye length (1/A): {kappa:.3f}__")
    
    # Model setup
    model = Model(z=z, z_index=z_index, hw=wall_zeros,
                  c_short=fluid_zeros, f1=fluid_zeros, f2=fluid_zeros)
    
    # Calculate potentials
    beta_phiw = fluid.beta * calc_phiw(z, n_point, n_component)
    beta_psi = fluid.beta * psi_0 * 1.0e-3  # Convert to mV (in Volts)
    beta_psi_charge = -beta_psi * fluid.charge
    
    # Bulk fluid correlation function
    CR_PATH = f"{Path.cwd().as_posix()}/data/{fluid.cr_filename}"
    
    try:
        c_short, _ = load_and_interpolate_cr(Path(CR_PATH), n_point, n_pair, z)
        
        # Calculate integrals
        f1 = integral_z_infty_dr_r_c_short(c_short, n_pair, n_point, z)
        f2 = integral_z_infty_dr_r2_c_short(c_short, n_pair, n_point, z)
        model.f1 = f1
        model.f2 = f2
        
        f1_integrand = calc_f1_integrand(c_short, n_pair, z, n_point)
        f2_integrand = calc_f2_integrand(c_short, n_pair, z, n_point)
        
        bulk_data_available = True
        
    except Exception as e:
        st.error(
            f"Bulk fluid correlation data not available at: {Path(CR_PATH).parent.as_posix()}"
        )
        bulk_data_available = False
    
    st.subheader(f"Charged fluids near an interface: {fluid.name}")
    
    tab1, tab2, tab0 = st.tabs(["Calculation", "Output graphs", "Bulk properties"])
    
    with tab1:
        # Use shared calculation component
        solution, hw_solution, calculation_run = create_calculation_tab(
            fluid=fluid,
            model=model,
            d=d,
            beta_phiw=beta_phiw,
            beta_psi_charge=beta_psi_charge,
            n_component=n_component,
            n_pair=n_pair,
            z=z,
            status="ready" if bulk_data_available else "in_progress"
        )
    
    with tab2:
        st.markdown("#")
        if calculation_run:
            if solution is not None:
                z_plots = {
                    "Solution: g_{wi}(z)": {
                        "fn_label": "g",
                        "plot_fn": hw_solution + 1,
                        "plot_name": "Solution: g(z)"
                    }
                }
                plot_wall_curves(n_component, z, z_plots, fluid.component)
            else:
                st.info("_Solution not found: no output available._")
        else:
            st.markdown("Select the Calculation tab and press the __[Run calculation]__ button.")
    
    with tab0:
        st.markdown("#")
        col1, col2, col3 = st.columns(3)
        col1.metric("Temperature (K)", f"{fluid.temperature}")
        col2.metric("Concentration (M/dm3)", f"{fluid.concentration[0]}")
        col3.metric("Components", f"{', '.join(fluid.component)}")
        
        if bulk_data_available:
            r_plots = {
                "c_short": {
                    "fn_label": "c_short",
                    "plot_fn": c_short,
                    "plot_name": "c_short(r)"
                },
                "f1": {
                    "fn_label": "f1",
                    "plot_fn": f1,
                    "plot_name": "f1(z)",
                    "xlim": [0, 10],
                    "ylim": None
                },
                "f2": {
                    "fn_label": "f2",
                    "plot_fn": f2,
                    "plot_name": "f2(z)",
                    "xlim": [0, 10],
                    "ylim": None
                },
                "f1_integrand": {
                    "fn_label": "f1_integrand",
                    "plot_fn": f1_integrand,
                    "plot_name": "f1_integrand(r)",
                    "xlim": [0, 10],
                    "ylim": None
                },
                "f2_integrand": {
                    "fn_label": "f2_integrand",
                    "plot_fn": f2_integrand,
                    "plot_name": "f2_integrand(r)",
                    "xlim": [0, 10],
                    "ylim": None
                }
            }
            plot_bulk_curves(n_component, z, r_plots, fluid.component)
        else:
            st.info("Bulk properties will be available once correlation data is loaded.")

else:
    st.error("Invalid choice of fluid")