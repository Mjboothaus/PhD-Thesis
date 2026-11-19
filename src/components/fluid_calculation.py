"""Shared calculation component for fluid pages."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import st_redirect as rd
import streamlit as st
from utils import get_memory_usage

from modelling import Model, calc_hw, opt_func, solve_model


def update_memory_display(container, initial_mem: float | None = None) -> None:
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
    beta_phiw: Any,
    beta_psi_charge: Any,
    n_component: int,
    n_pair: int,
    z: Any,
    status: str = "ready",
) -> tuple[Any, Any, bool]:
    """Create calculation tab with Newton-Krylov solver UI.

    Returns
    -------
    solution, hw_solution, calculation_run
    """
    solution = None
    hw_solution = None
    # Initial guess is wall indirect correlation tw (issue #21)
    tw_initial = model.tw

    st.markdown("### Newton-Krylov Algorithm Solver")

    if status != "ready":
        st.warning("This fluid type is currently under development")
        return None, None, False

    mem_container = st.container()
    update_memory_display(mem_container)

    if run_calc := st.button("Run calculation", key="run_calc_btn"):
        initial_memory = get_memory_usage()
        solver_output_filepath = f"{Path.cwd()}/data/solver_out.txt"

        try:
            with st.spinner("Finding optimal solution..."):
                st.markdown("__Solver progress:__")
                output_area = st.empty()

                with rd.stdout(
                    to=output_area,
                    to_file=solver_output_filepath,
                    format="text",
                    max_buffer=10000,
                ):
                    solution = solve_model(
                        opt_func,
                        tw_initial,
                        fluid,
                        model,
                        d,
                        beta_phiw,
                        beta_psi_charge,
                    )
                    tw_solution = solution.x

                # Persist on shared Model; report iteration count clearly (issue #14)
                model.tw = tw_solution
                hw_solution = calc_hw(tw_solution, n_component, beta_phiw)
                model.hw = hw_solution
                nit = int(solution.get("nit", 0))
                st.success(
                    f"Converged after {nit} iterations. SciPy: {solution.get('message', '')}"
                )
                st.markdown("__Convergence Plot:__")
                from plotting import plot_convergence

                plot_convergence(solver_output_filepath)
                update_memory_display(mem_container, initial_memory)

        except ValueError as err_message:
            st.error("Solver failed to find a solution")
            st.warning(str(err_message))
            hw_solution = model.hw
            update_memory_display(mem_container, initial_memory)

    calculation_run = bool(run_calc)
    return solution, hw_solution, calculation_run
