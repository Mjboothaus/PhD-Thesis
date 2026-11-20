"""Shared calculation component for fluid pages."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import st_redirect as rd
import streamlit as st
from utils import get_memory_usage

from charged_fluids import Model, calc_hw, opt_func, solve_model

# Fluids the wall solver may attempt (same as Marimo / wall_ready_symbols)
_SOLVE_STATUSES = frozenset({"ready", "experimental"})


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
    """Create calculation tab with Newton–Krylov solver UI.

    Returns
    -------
    solution, hw_solution, calculation_run
    """
    solution = None
    hw_solution = None
    tw_initial = model.tw

    st.markdown("### Newton–Krylov algorithm solver")
    st.caption("Uses `charged_fluids.solve_model` (same core as Marimo / tests).")

    if status not in _SOLVE_STATUSES:
        st.warning(f"Fluid status `{status}` — wall solve is not available.")
        return None, None, False

    if status == "experimental":
        st.info(
            "This fluid is **experimental**: the solver will run, but bulk "
            "thermodynamics / residual quality may still be imperfect."
        )

    mem_container = st.container()
    update_memory_display(mem_container)

    if run_calc := st.button("Run calculation", key="run_calc_btn"):
        initial_memory = get_memory_usage()
        solver_output_filepath = str(Path.cwd() / "data" / "solver_out.txt")

        try:
            with st.spinner("Finding optimal solution…"):
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
                        tw_initial.copy(),
                        fluid,
                        model,
                        d,
                        beta_phiw,
                        beta_psi_charge,
                        track_history=True,
                    )
                    tw_solution = solution.x

                model.tw = tw_solution
                hw_solution = calc_hw(tw_solution, n_component, beta_phiw)
                model.hw = hw_solution
                nit = int(getattr(solution, "nit", 0) or 0)
                ok = bool(getattr(solution, "success", False))
                msg = str(getattr(solution, "message", "") or "")
                summary = getattr(solution, "convergence_summary", {}) or {}
                f_final = summary.get("final_residual")
                extra = f" · final |F|₂={f_final:.3e}" if f_final is not None else ""
                if ok:
                    st.success(f"Converged after {nit} iterations{extra}. SciPy: {msg}")
                else:
                    st.warning(
                        f"Finished after {nit} iterations (success={ok}){extra}. SciPy: {msg}"
                    )
                st.markdown("__Convergence plot:__")
                from plotting import plot_convergence

                plot_convergence(
                    solver_output_filepath,
                    history=getattr(solution, "convergence_history", None),
                )
                update_memory_display(mem_container, initial_memory)

        except ValueError as err_message:
            st.error("Solver failed to find a solution")
            st.warning(str(err_message))
            hw_solution = model.hw
            update_memory_display(mem_container, initial_memory)
        except Exception as exc:  # noqa: BLE001 — surface unexpected solver errors in UI
            st.exception(exc)
            hw_solution = model.hw
            update_memory_display(mem_container, initial_memory)

    calculation_run = bool(run_calc)
    return solution, hw_solution, calculation_run
