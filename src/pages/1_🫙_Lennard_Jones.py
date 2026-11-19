"""Lennard-Jones fluid near an interface (Streamlit page)."""

import streamlit as st
from components.fluid_calculation import create_calculation_tab
from plotting import plot_bulk_curves, plot_wall_curves
from sidebar import create_sidebar

from charged_fluids import prepare_fluid_run, set_fluid_parameters

fluid_symbol = "lj1"

_fluid_preview = set_fluid_parameters(fluid_symbol)
if _fluid_preview is None:
    st.error("Invalid choice of fluid")
    st.stop()

z_cutoff, n_point, psi_0, tolerance, max_iteration = create_sidebar(_fluid_preview)

try:
    run = prepare_fluid_run(
        fluid_symbol,
        z_cutoff=float(z_cutoff),
        n_point=int(n_point),
        tolerance=float(tolerance),
        max_iteration=int(max_iteration),
        psi_0_mv=float(psi_0),
    )
except FileNotFoundError as exc:
    st.error(f"Bulk c(r) data missing: {exc}")
    st.stop()
except Exception as exc:
    st.exception(exc)
    st.stop()

fluid = run.fluid
model = run.model
d = run.discrete
z = model.z
n_component = fluid.n_component
n_pair = fluid.n_pair

st.subheader(f"Fluids near an interface: {fluid.name}")
if fluid.status != "ready":
    st.info(f"Status: **{fluid.status}**")

tab_calc, tab_out, tab_bulk = st.tabs(["Calculation", "Output graphs", "Bulk properties"])

with tab_calc:
    solution, hw_solution, calculation_run = create_calculation_tab(
        fluid=fluid,
        model=model,
        d=d,
        beta_phiw=run.beta_phiw,
        beta_psi_charge=run.beta_psi_charge,
        n_component=n_component,
        n_pair=n_pair,
        z=z,
        status="ready" if fluid.status == "ready" else fluid.status,
    )

with tab_out:
    st.markdown("#")
    if calculation_run:
        if solution is not None and hw_solution is not None:
            z_plots = {
                "Solution: g_{wi}(z)": {
                    "fn_label": "g",
                    "plot_fn": hw_solution + 1,
                    "plot_name": "g(z)",
                }
            }
            plot_wall_curves(n_component, z, z_plots, fluid.component)
        else:
            st.info("_Solution not found: no output available._")
    else:
        st.markdown("Select the Calculation tab and press **Run calculation**.")

with tab_bulk:
    st.markdown("#")
    col1, col2, col3 = st.columns(3)
    col1.metric("Temperature (K)", f"{fluid.temperature}")
    col2.metric("Concentration (M/dm3)", f"{fluid.concentration[0]}")
    col3.metric("Components", f"{', '.join(fluid.component)}")
    st.caption(f"Bulk c(r) file: `{run.cr_path.name}`")

    # Issue #3: plot beta*u so the bulk potential is visible on a sensible scale
    r_plots = {
        "c_short": {
            "fn_label": "c_short",
            "plot_fn": model.c_short,
            "plot_name": "c_short(r)",
        },
        "beta*u": {
            "fn_label": "beta u",
            "plot_fn": run.beta_u,
            "plot_name": r"$\beta u(r)$",
            "xlim": [0, 10],
            "ylim": [-1, 10],
        },
        "f1": {
            "fn_label": "f1",
            "plot_fn": model.f1,
            "plot_name": "f1(z)",
            "xlim": [0, 10],
            "ylim": None,
        },
        "f2": {
            "fn_label": "f2",
            "plot_fn": model.f2,
            "plot_name": "f2(z)",
            "xlim": [0, 10],
            "ylim": None,
        },
        "f1_integrand": {
            "fn_label": "f1_integrand",
            "plot_fn": run.f1_integrand,
            "plot_name": "f1_integrand(r)",
            "xlim": [0, 10],
            "ylim": None,
        },
        "f2_integrand": {
            "fn_label": "f2_integrand",
            "plot_fn": run.f2_integrand,
            "plot_name": "f2_integrand(r)",
            "xlim": [0, 10],
            "ylim": None,
        },
    }
    plot_bulk_curves(n_component, z, r_plots, fluid.component)
