"""Shared calculation component for fluid pages."""
from pathlib import Path
from typing import Optional, Tuple, Dict, Any

import streamlit as st
import st_redirect as rd
from modelling import (
    Model, calc_hw, solve_model, opt_func, 
    load_and_interpolate_cr
)
from utils import get_memory_usage

def update_memory_display(container, initial_mem: Optional[float] = None) -> None:
    """Update memory usage display in the given container."""
    col1, col2 = container.columns(2)
    current_mem = get_memory_usage()
    col1.metric("Current Memory Usage", f"{current_mem:.2f} MB")
    if initial_mem is not None:
        col2.metric("Memory Change", f"{current_mem - initial_mem:+.2f} MB")

def create_calculation_tab(
    fluid: Any,
    model: Model,
    d: Any,
    beta_phiw: float,
    beta_psi_charge: float,
    n_component: int,
    n_pair: int,
    z: Any,
    status: str = "ready"
) -> Tuple[Dict, Any, bool]:
    """Create calculation tab with solver functionality.
    
    Args:
        fluid: Fluid object containing properties
        model: Model object for calculations
        d: Numerical parameters
        beta_phiw: Wall potential parameter
        beta_psi_charge: Charge parameter
        n_component: Number of components
        n_pair: Number of pairs
        z: Distance array
        status: Status of the fluid type ("ready", "in_progress", "not_implemented")
        
    Returns:
        Tuple of (solution dict, hw_solution array)
    """
    solution = None
    hw_solution = None
    tw_initial = model.hw  # Initial guess from model

    st.markdown("### Newton-Krylov Algorithm Solver")
    
    if status != "ready":
        st.warning("⚠️ This fluid type is currently under development")
        return None, None

    # Memory display
    mem_container = st.container()
    update_memory_display(mem_container)

    # Run calculation button
    if run_calc := st.button("Run calculation", key="run_calc_btn"):
        initial_memory = get_memory_usage()
        solver_output_filepath = f"{Path.cwd()}/data/solver_out.txt"

        try:
            with st.spinner("Finding optimal solution..."):
                st.markdown("__Solver progress:__")
                output_area = st.empty()

                with rd.stdout(to=output_area, to_file=solver_output_filepath, format="text", max_buffer=10000):
                    solution = solve_model(
                        opt_func, tw_initial, fluid, model, d,
                        beta_phiw, beta_psi_charge
                    )
                    tw_solution = solution.x

                # Update result display
                result_msg = solution["message"].replace(".", f" after {solution['nit']} iterations.")
                result_msg = result_msg.replace("A s", "S")
                st.success(result_msg)

                # Calculate solution and plot convergence
                hw_solution = calc_hw(tw_solution, n_component, beta_phiw)
                st.markdown("__Convergence Plot:__")
                from plotting import plot_convergence
                plot_convergence(solver_output_filepath)

                # Final memory update
                update_memory_display(mem_container, initial_memory)

        except ValueError as err_message:
            st.error("Solver failed to find a solution")
            st.warning(str(err_message))
            hw_solution = tw_initial
            # Still update memory display
            update_memory_display(mem_container, initial_memory)

    # Return solution, hw_solution, and whether calculation was run
    calculation_run = bool(run_calc)
    return solution, hw_solution, calculation_run
