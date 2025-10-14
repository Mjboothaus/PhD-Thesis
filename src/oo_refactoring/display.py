import streamlit as st


from optimization_display_component import OptimizationDisplay

class DisplayManager:
    @staticmethod
    def display_bulk_properties(fluid):
        """Display bulk properties in Streamlit."""
        st.write("### Bulk Properties")
        
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Temperature", f"{fluid.temperature:.2f} K")
        with col2:
            st.metric("Concentration", f"{fluid.concentration:.3f} M")
        with col3:
            st.metric("Components", len(fluid.components))
            
        with st.expander("Detailed Properties"):
            st.write("#### System Components")
            components_df = pd.DataFrame({
                'Name': fluid.components,
                'Charge': fluid.charge if hasattr(fluid, 'charge') else ['N/A'] * len(fluid.components),
                'Density': fluid.rho if hasattr(fluid, 'rho') else ['N/A'] * len(fluid.components)
            })
            st.dataframe(components_df)

    @staticmethod
    def display_results(result):
        """Display optimization results with enhanced visualization."""
        if not result:
            st.warning("No results available yet. Run a calculation first.")
            return
            
        # Create tabs for different views
        tab1, tab2, tab3 = st.tabs(["Optimization", "Physical Properties", "Raw Data"])
        
        with tab1:
            OptimizationDisplay.display_optimization_status(result)
            
        with tab2:
            if hasattr(result, 'physical_properties'):
                st.write("### Physical Properties")
                for prop, value in result.physical_properties.items():
                    if isinstance(value, (int, float)):
                        st.metric(
                            prop.replace('_', ' ').title(),
                            f"{value:.6f}"
                        )
            else:
                st.info("No physical properties available.")
                
        with tab3:
            with st.expander("Raw Optimization Results"):
                st.json({
                    "success": bool(result.success),
                    "status": result.status,
                    "message": result.message,
                    "nit": int(result.nit),
                    "nfev": int(result.nfev),
                    "x_shape": list(result.x.shape),
                    "fun_shape": list(result.fun.shape)
                })
