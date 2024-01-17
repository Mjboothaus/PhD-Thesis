# src/oo_refactoring/base_streamlit_app.py
from pathlib import Path
from time import sleep

import streamlit as st

from .helper_functions import calc_hw, load_and_interpolate_cr, solve_model
from .numerics import set_num_parameters
from .parameters import fluid_specific_parameters, set_fluid_parameters
from .plotting import plot_bulk_curves, plot_convergence, plot_wall_curves
from .sidebar import create_sidebar


class BaseStreamlitApp:
    def __init__(self, fluid_symbol):
        self.fluid_symbol = fluid_symbol
        self.fluid = set_fluid_parameters(fluid_symbol)
        self.other_params = fluid_specific_parameters(fluid_symbol)
        self.configure_sidebar()
        self.configure_model()

    def configure_sidebar(self):
        if self.fluid is not None:
            (
                self.z_cutoff,
                self.n_point,
                self.psi_0,
                self.tolerance,
                self.max_iteration,
            ) = create_sidebar(self.fluid)
        else:
            st.error("Invalid choice of fluid")

    def configure_model(self):
        self.d = set_num_parameters(
            self.n_point,
            self.z_cutoff,
            self.fluid.n_component,
            self.fluid.n_pair,
            self.tolerance,
            self.max_iteration,
        )
        self.n_component = self.d.n_component
        self.n_pair = self.d.n_pair

    def run(self):
        # Check if fluid parameters are set correctly
        if self.fluid is None:
            st.error("Invalid choice of fluid")
            return

        # Set up the Streamlit interface
        st.subheader(f"{self.fluid.name} near an interface")

        # Set up tabs for different sections of the app
        tab_bulk, tab_calc, tab_results = st.tabs(
            ["Bulk properties", "Calculation", "Results"]
        )

        with tab_bulk:
            self.display_bulk_properties()

        with tab_calc:
            self.perform_calculation()

        with tab_results:
            self.display_results()

    def display_bulk_properties(self):
        # Placeholder for displaying bulk properties
        pass

    def perform_calculation(self):
        # Placeholder for calculation logic
        pass

    def display_results(self):
        # Placeholder for displaying results
        pass


# The placeholders `display_bulk_properties`, `perform_calculation`, and `display_results`
# should be implemented with the specific logic for displaying the bulk properties,
# performing the calculations, and displaying the results, respectively.
