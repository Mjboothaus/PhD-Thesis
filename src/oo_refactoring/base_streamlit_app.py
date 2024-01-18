from oo_refactoring.modelling import ModelManager
from oo_refactoring.numerics import NumericsManager
from oo_refactoring.display import DisplayManager
import streamlit as st


class BaseStreamlitApp:
    def __init__(self, fluid):
        self.fluid = fluid
        self.model = None
        self.results = None

    def run(self):
        # Set up the Streamlit interface
        st.subheader(f"{self.fluid.name} near an interface")

        # Set up tabs for different sections of the app
        tab_bulk, tab_calc, tab_results = st.tabs(
            ["Bulk properties", "Calculation", "Results"]
        )

        with tab_bulk:
            DisplayManager.display_bulk_properties(self.fluid)

        with tab_calc:
            self.model = ModelManager.initialize_model(self.fluid)
            self.results = NumericsManager.perform_calculation(self.model, self.fluid)

        with tab_results:
            DisplayManager.display_results(self.results)
