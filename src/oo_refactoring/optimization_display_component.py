"""
A Streamlit component for displaying optimization results.
"""

import streamlit as st
import pandas as pd
import plotly.express as px
import numpy as np

class OptimizationDisplay:
    """
    A class to handle the display of optimization results in the Streamlit app.
    This provides an enhanced view of the optimization process with interactive plots
    and detailed metrics.
    """

    @staticmethod
    def display_optimization_status(result):
        """Display the optimization results with detailed metrics and visualizations."""
        
        # Basic status
        status = "🟢 Success" if result.success else "🔴 Failed"
        st.write(f"### Optimization Status: {status}")
        
        # Display summary metrics in columns
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.metric(
                "Total Iterations",
                result.convergence_summary.get('total_iterations', 'N/A')
            )
        
        with col2:
            st.metric(
                "Final Residual",
                f"{result.convergence_summary.get('final_residual', 'N/A'):.2e}"
            )
        
        with col3:
            st.metric(
                "Final Step Size",
                f"{result.convergence_summary.get('final_step_size', 'N/A'):.2e}"
            )
        
        # Add message if optimization failed
        if not result.success:
            st.error(f"Optimization failed: {result.message}")
        
        # Show convergence plots
        st.write("### Convergence Plots")
        result.convergence_history.plot_convergence()
        
        # Detailed convergence metrics
        with st.expander("Detailed Convergence Metrics"):
            metrics_df = pd.DataFrame({
                'Iteration': result.convergence_history.iterations,
                'Residual': result.convergence_history.residuals,
                'Function Value': result.convergence_history.function_values,
                'Step Size': result.convergence_history.step_sizes,
                **{k: v for k, v in result.convergence_history.convergence_metrics.items()}
            })
            st.dataframe(metrics_df.round(8))
            
            # Add download button for metrics
            csv = metrics_df.to_csv(index=False)
            st.download_button(
                "Download Metrics CSV",
                csv,
                "optimization_metrics.csv",
                "text/csv",
                key='download-optimization-metrics'
            )
        
        # Warning messages and recommendations
        if result.convergence_summary.get('total_iterations', 0) > result.x.size * 2:
            st.warning(
                "⚠️ Optimization took more iterations than expected. Consider adjusting parameters:"
                "\n- Increase tolerance"
                "\n- Try different initial conditions"
                "\n- Check problem scaling"
            )
        
        if result.convergence_summary.get('final_residual', np.inf) > 1e-3:
            st.warning(
                "⚠️ Final residual is relatively high. Solution may not be optimal."
                " Consider:"
                "\n- Increasing max iterations"
                "\n- Adjusting problem formulation"
                "\n- Using a different optimization method"
            )
        
        # Additional information
        with st.expander("Solution Details"):
            st.json({
                "Status": result.status,
                "Success": bool(result.success),
                "Message": result.message,
                "Number of iterations": int(result.nit),
                "Number of function evaluations": int(result.nfev),
                "Final cost value": float(result.fun.max()),
                "Optimization method": "Krylov"
            })