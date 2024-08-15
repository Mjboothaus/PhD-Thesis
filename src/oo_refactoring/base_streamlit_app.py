from modelling import ModelManager
from numerics import NumericsManager
from display import DisplayManager
import streamlit as st
from sidebar import Sidebar


class BaseStreamlitApp:
    def __init__(self, fluid):
        self.fluid = fluid
        self.model = None  # ModelManager()
        self.numerics = (
            None  # NumericsManager() - these depend on Sidebar values - have defaults
        )
        self.results = None

    def run(self):
        # Set up the Streamlit interface
        st.subheader(f"{self.fluid} near an interface")
        Sidebar.create_sidebar(self.fluid)

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
