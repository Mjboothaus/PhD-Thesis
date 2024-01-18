# src/oo_refactoring/display.py
import streamlit as st


class DisplayManager:
    @staticmethod
    def display_bulk_properties(fluid):
        # Streamlit code to display bulk properties
        st.write("Bulk Properties:")
        st.write(f"Temperature: {fluid.temperature} K")
        st.write(f"Concentration: {fluid.concentration} M")
        # ... additional properties ...

    @staticmethod
    def display_results(results):
        # Streamlit code to display results
        st.write("Calculation Results:")
        # Assuming 'results' is a dictionary with result data
        for key, value in results.items():
            st.write(f"{key}: {value}")
        # ... additional results ...
