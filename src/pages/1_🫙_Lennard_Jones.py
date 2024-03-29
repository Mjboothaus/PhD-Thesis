from pathlib import Path
from time import sleep

import pandas as pd
import st_redirect as rd
import streamlit as st
from modelling import *
from numerics import set_num_parameters
from parameters import fluid_specific_parameters, set_fluid_parameters
from plotting import plot_bulk_curves, plot_convergence, plot_wall_curves
from sidebar import create_sidebar

# Initialise fluid and numerical parameters

fluid_symbol = "lj2"

fluid = set_fluid_parameters(fluid_symbol)
if fluid is not None:
    z_cutoff, n_point, psi_0, tolerance, max_iteration = create_sidebar(fluid)
else:
    st.error("Invalid choice of fluid")

other_params = fluid_specific_parameters(fluid_symbol)

epsilon_lj = other_params[fluid_symbol]["epsilon_lj"]
sigma_lj = other_params[fluid_symbol]["sigma_lj"]

d = set_num_parameters(n_point, z_cutoff, fluid.n_component,
                       fluid.n_pair, tolerance, max_iteration)

n_component = d.n_component
n_pair = d.n_pair

z = np.linspace(0.0, z_cutoff, n_point)
z_index = np.arange(0, n_point, dtype=int)
wall_zeros = np.zeros((n_point, n_component))
fluid_zeros = np.zeros((n_point, n_pair))

fluid.beta = calc_beta(fluid.temperature)
fluid.epsilon = calc_epsilon(fluid.epsilon_r)
fluid.charge = calc_charge(fluid.valence)
fluid.charge_pair = calc_charge_pair(
    fluid.beta, fluid.charge, fluid.epsilon, n_component, n_pair)
fluid.rho = calc_rho(fluid.concentration)

kappa = calc_kappa(fluid.beta, fluid.charge, fluid.rho, fluid.epsilon)

if fluid.symbol not in ["lj1", "lj2"]:  # kappa only relevant for charged fluids
    st.sidebar.markdown(f"__kappa: {kappa}__")

# st.sidebar.text(fluid.cr_filename)

r = z    # r on same discretisation as z
beta_u = fluid.beta * \
    calc_u_lj(epsilon_lj, sigma_lj, n_point, n_component, n_pair, r)

ur = calc_u_lj(epsilon_lj, sigma_lj, n_point, n_component, n_pair, r)


beta_phiw = fluid.beta * calc_phiw(z, n_point, n_component)
beta_psi = 0.0
beta_psi_charge = -beta_psi * fluid.charge

model = Model(z=z, z_index=z_index, hw=wall_zeros,
              c_short=fluid_zeros, f1=fluid_zeros, f2=fluid_zeros)

# Bulk-fluid inputs (direct correlation function

# TODO: Run bulk calculation on the fly (or from stored cr)

CR_PATH = f"{Path.cwd().as_posix()}/data/{fluid.cr_filename}"

try:
    c_short, _ = load_and_interpolate_cr(Path(CR_PATH), n_point, n_pair, z)
except Exception as e:
    st.error(
        f"Please ensure the c(r) data file is available here: {Path(CR_PATH).parent.as_posix()}")
    st.info("Restarting in 10 seconds")
    st.markdown("##")
    st.exception(e)
    sleep(10.0)
    st.experimental_rerun()

f1 = integral_z_infty_dr_r_c_short(c_short, n_pair, n_point, z)
f2 = integral_z_infty_dr_r2_c_short(c_short, n_pair, n_point, z)
model.f1 = f1
model.f2 = f2

f1_integrand = calc_f1_integrand(c_short, n_pair, z, n_point)
f2_integrand = calc_f2_integrand(c_short, n_pair, z, n_point)

# initial guess of zero - maybe should be \beta \phi

tw_initial = np.zeros((n_point, n_component))
hw_initial = calc_hw(tw_initial, n_component, beta_phiw)

if fluid.symbol in ["lj1", "lj2"]:
    st.subheader(f"Fluids near an interface: {fluid.name}")
else:
    st.subheader(f"Charged fluids near an interface: {fluid.name}")

tab0, tab1, tab2 = st.tabs(["Bulk properties", "Calculation", "Output graphs"])

with tab0:
    st.markdown("#")
    col1, col2, col3 = st.columns(3)
    col1.metric("Temperature (K)", f"{fluid.temperature}")
    col2.metric("Concentration (M/dm3)", f"{fluid.concentration[0]}")
    col3.metric("Components", f"{', '.join(fluid.component)}")

    r_plots = dict({"c_short": dict({"fn_label": "c_short", "plot_fn": c_short,
                                    "plot_name": "c_short(r)"})})
    r_plots["u(r)"] = dict({"fn_label": "u(r)", "plot_fn": ur, "plot_name": r'u(r)',
                            "xlim": [0, 10], "ylim": [-1, 10]})
    r_plots["f1"] = dict({"fn_label": "f1", "plot_fn": f1,
                        "plot_name": "f1(z)",  "xlim": [0, 10], "ylim": None})
    r_plots["f2"] = dict({"fn_label": "f2", "plot_fn": f2,
                        "plot_name": "f2(z)",  "xlim": [0, 10], "ylim": None})

    r_plots["f1_integrand"] = dict({"fn_label": "f1_integrand", "plot_fn": f1_integrand,
                                "plot_name": "f1_integrand(r)",  "xlim": [0, 10], "ylim": None})
    r_plots["f2_integrand"] = dict({"fn_label": "f2_integrand", "plot_fn": f2_integrand,
                                "plot_name": "f2_integrand(r)",  "xlim": [0, 10], "ylim": None})

    plot_bulk_curves(n_component, r, r_plots, fluid.component)

with tab1:
    # Run calculation
    st.markdown("#")
    if run_calc := st.button("Run calculation"):
        with st.spinner("__Finding optimal solution:__"):
            st.markdown("Solver output (Newton-Krylov algorithm)")
            to_out = st.empty()
            solver_output_filepath = f"{Path.cwd()}/data/solver_out.txt"
            with rd.stdout(to=to_out, to_file=solver_output_filepath, format="text", max_buffer=10000):
                    # Solve non-linear equation
                    try:
                        solution = solve_model(opt_func, tw_initial, fluid, model, d,
                                            beta_phiw, beta_psi_charge)
                        tw_solution = solution.x
                    except ValueError as err_message:
                        solution = None
                        st.info(err_message)
        if solution is not None:
            st.write(solution["message"].replace(".", " after " + str(solution["nit"]) + " iterations.").replace("A s", "S"))
            hw_solution = calc_hw(tw_solution, n_component, beta_phiw)
            plot_convergence(solver_output_filepath)
        else:
            st.error("Solver failed to find a solution: see error message above.")
            hw_solution = hw_initial

with tab2:
        st.markdown("#")
        if run_calc:
            if solution is not None:
                z_plots = dict({"Solution: g_{wi}(z)": dict({"fn_label": "g", 
                                                    "plot_fn": hw_solution+1,
                                                    "plot_name": "Solution: g(z)"})})
                plot_wall_curves(n_component, z, z_plots, fluid.component)
            else:
                st.info("_Solution not found: no output available._")
        else:
            st.markdown("Select the Calculation tab and press the __[Run calculation]__ button.")