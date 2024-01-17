from pathlib import Path
from time import sleep

import numpy as np
import st_redirect as rd
import streamlit as st

from base_streamlit_app import BaseStreamlitApp
from modelling import (
    calc_beta,
    calc_charge,
    calc_charge_pair,
    calc_epsilon,
    calc_f1_integrand,
    calc_f2_integrand,
    calc_hw,
    calc_kappa,
    calc_phiw,
    calc_rho,
    calc_u_lj,
    integral_z_infty_dr_r2_c_short,
    integral_z_infty_dr_r_c_short,
    load_and_interpolate_cr,
    solve_model,
)
from plotting import plot_bulk_curves, plot_convergence, plot_wall_curves


class LennardJonesApp(BaseStreamlitApp):
    def __init__(self, fluid_symbol):
        super().__init__(fluid_symbol)
        self.epsilon_lj = self.other_params[self.fluid_symbol]["epsilon_lj"]
        self.sigma_lj = self.other_params[self.fluid_symbol]["sigma_lj"]

    def run(self):
        # First, call the base class's run method to set up the common interface
        super().run()

        # Fluid-specific logic

        # Set up the model parameters
        self.fluid.beta = calc_beta(self.fluid.temperature)
        self.fluid.epsilon = calc_epsilon(self.fluid.epsilon_r)
        self.fluid.charge = calc_charge(self.fluid.valence)
        self.fluid.charge_pair = calc_charge_pair(
            self.fluid.beta,
            self.fluid.charge,
            self.fluid.epsilon,
            self.n_component,
            self.n_pair,
        )
        self.fluid.rho = calc_rho(self.fluid.concentration)

        # Set up the spatial discretization
        z = np.linspace(0.0, self.z_cutoff, self.n_point)
        z_index = np.arange(0, self.n_point, dtype=int)
        wall_zeros = np.zeros((self.n_point, self.n_component))
        fluid_zeros = np.zeros((self.n_point, self.n_pair))

        # Calculate interaction potentials
        r = z  # r on same discretisation as z
        beta_u = self.fluid.beta * calc_u_lj(
            self.epsilon_lj,
            self.sigma_lj,
            self.n_point,
            self.n_component,
            self.n_pair,
            r,
        )

        # Initialize the model
        model = Model(
            z=z,
            z_index=z_index,
            hw=wall_zeros,
            c_short=fluid_zeros,
            f1=fluid_zeros,
            f2=fluid_zeros,
        )

        # Load and interpolate c(r) data
        CR_PATH = f"{Path.cwd().as_posix()}/data/{self.fluid.cr_filename}"
        try:
            c_short, _ = load_and_interpolate_cr(
                Path(CR_PATH), self.n_point, self.n_pair, z
            )
        except Exception as e:
            st.error(
                f"Please ensure the c(r) data file is available here: {Path(CR_PATH).parent.as_posix()}"
            )
            st.info("Restarting in 10 seconds")
            st.markdown("##")
            st.exception(e)
            sleep(10.0)
            st.experimental_rerun()

        # Calculate integrals
        f1 = integral_z_infty_dr_r_c_short(c_short, self.n_pair, self.n_point, z)
        f2 = integral_z_infty_dr_r2_c_short(c_short, self.n_pair, self.n_point, z)
        model.f1 = f1
        model.f2 = f2

        # ... (rest of the code from the original script)

        # Note: The rest of the code from the original script should be refactored
        # and included here, following the same pattern as above.

        # ... (previous code from the run method)

        # Set up Streamlit interface (in base class run method)
        # st.subheader(f"Fluids near an interface: {self.fluid.name}")
        # tab0, tab1, tab2 = st.tabs(["Bulk properties", "Calculation", "Output graphs"])

        with tab0:
            # Display bulk properties
            self.display_bulk_properties()

        with tab1:
            # Run calculation
            self.perform_calculation(model)

        with tab2:
            # Output graphs
            self.display_output_graphs(model)

    def display_bulk_properties(self):
        st.markdown("#")
        col1, col2, col3 = st.columns(3)
        col1.metric("Temperature (K)", f"{self.fluid.temperature}")
        col2.metric("Concentration (M/dm3)", f"{self.fluid.concentration[0]}")
        col3.metric("Components", f"{', '.join(self.fluid.component)}")
        # ... (additional code to display bulk properties)

    def perform_calculation(self, model):
        if st.button("Run calculation"):
            with st.spinner("Finding optimal solution..."):
                # ... (code to perform the calculation)
                # This will involve calling solve_model and handling the output
                pass

    def display_output_graphs(self, model):
        if "run_calc" in st.session_state and st.session_state.run_calc:
            # ... (code to display output graphs)
            # This will involve calling plot_wall_curves and handling the output
            pass
