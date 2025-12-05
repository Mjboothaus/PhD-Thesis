"""
This module handles the display of optimization-related information in the Streamlit app.
"""

import streamlit as st
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
import pandas as pd


class OptimizationDisplay:
    @staticmethod
    def display_convergence_metrics(convergence_history: dict) -> None:
        """Display convergence metrics in the Streamlit app."""
        st.subheader("Convergence Metrics")
        
        # Create metrics display
        col1, col2, col3 = st.columns(3)
        col1.metric("Total Iterations", len(convergence_history['iterations']))
        col2.metric("Final Norm", f"{convergence_history['norm_dsqn'][-1]:.2e}")
        col3.metric("Final Mean Change", f"{convergence_history['mean_change'][-1]:.2e}")
        
        # Create convergence plots
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=("Norm DSQN", "Max DSQN", "Mean Change", "All Metrics (Log Scale)"),
            vertical_spacing=0.12
        )
        
        # Individual metric plots
        fig.add_trace(
            go.Scatter(x=convergence_history['iterations'], y=convergence_history['norm_dsqn'],
                      name="Norm DSQN", mode='lines+markers'),
            row=1, col=1
        )
        
        fig.add_trace(
            go.Scatter(x=convergence_history['iterations'], y=convergence_history['max_dsqn'],
                      name="Max DSQN", mode='lines+markers'),
            row=1, col=2
        )
        
        fig.add_trace(
            go.Scatter(x=convergence_history['iterations'], y=convergence_history['mean_change'],
                      name="Mean Change", mode='lines+markers'),
            row=2, col=1
        )
        
        # Combined log plot
        for metric in ['norm_dsqn', 'max_dsqn', 'mean_change']:
            fig.add_trace(
                go.Scatter(x=convergence_history['iterations'], 
                          y=convergence_history[metric],
                          name=metric.replace('_', ' ').title(),
                          mode='lines'),
                row=2, col=2
            )
        
        fig.update_layout(height=800, showlegend=True)
        fig.update_yaxes(type="log", row=2, col=2)
        
        st.plotly_chart(fig, use_container_width=True)
        
        # Display convergence data table
        df = pd.DataFrame(convergence_history)
        df = df.round(8)
        st.dataframe(df, use_container_width=True)

    @staticmethod
    def display_optimization_status(solver_status: dict) -> None:
        """Display optimization status information."""
        st.subheader("Optimization Status")
        
        status_color = "🟢" if solver_status['converged'] else "🔴"
        st.write(f"{status_color} Status: {'Converged' if solver_status['converged'] else 'Not Converged'}")
        
        # Display solver metrics
        metrics_df = pd.DataFrame([solver_status['metrics']])
        st.dataframe(metrics_df, use_container_width=True)
        
        # Display any warnings or messages
        if solver_status.get('warnings'):
            st.warning("\\n".join(solver_status['warnings']))
        
        if solver_status.get('messages'):
            st.info("\\n".join(solver_status['messages']))

    @staticmethod
    def display_solution_properties(properties: dict) -> None:
        """Display solution properties."""
        st.subheader("Solution Properties")
        
        # Create expandable sections for different property types
        with st.expander("Physical Properties", expanded=True):
            for prop, value in properties.get('physical', {}).items():
                st.metric(prop, f"{value:.6f}")
        
        with st.expander("Statistical Properties"):
            for prop, value in properties.get('statistical', {}).items():
                st.metric(prop, f"{value:.6f}")
        
        with st.expander("Raw Properties Data"):
            st.json(properties)