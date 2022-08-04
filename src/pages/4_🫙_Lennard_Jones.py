import sys
from pathlib import Path
from time import sleep

import pandas as pd
#import st_redirect as rd
import st_redirect_new as rd
import streamlit as st
from modelling import *
from numerics import create_sidebar, set_num_parameters
from plotting import plotly_line

# from helper_functions import read_render_markdown_file

# Local plotting functions

def plot_selected_wall_curves(n_component, z, z_plots):
    for plot_title, details in z_plots.items():
        y_col_labels = ["z"]
        fn_label = details["fn_label"]
        y_col_labels.extend(f"{fn_label}_{i}" for i in range(n_component))
        try:
            xliml=details["xlim"]
        except Exception:
            xliml=None
        try:
            yliml=details["ylim"]
        except Exception:
            yliml=None

        fig = plotly_line(z, details["plot_fn"], y_col_labels, y_label=details["plot_name"], legend_label="",
                             xliml=xliml, yliml=yliml, title=plot_title)
        st.plotly_chart(fig)


def plot_selected_bulk_curves(n_component, r, r_plots):
    for plot_title, details in r_plots.items():
        y_col_labels = ["r"]
        fn_label = details["fn_label"]
        for i in range(n_component):
            y_col_labels.extend(f"{fn_label}_{i}{j}" for j in range(i, n_component))
        try:
            xliml=details["xlim"]
        except Exception:
            xliml=None
        try:
            yliml=details["ylim"]
        except Exception:
            yliml=None

        fig = plotly_line(r, details["plot_fn"], y_col_labels, y_label=details["plot_name"], legend_label="",
                            xliml=xliml, yliml=yliml, title=plot_title)
        st.plotly_chart(fig)


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

d = set_num_parameters(n_point, z_cutoff, fluid.n_component, fluid.n_pair, tolerance, max_iteration)

n_component = d.n_component
n_pair = d.n_pair

z = np.linspace(0.0, z_cutoff, n_point)
z_index = np.arange(0, n_point, dtype=int)
wall_zeros = np.zeros((n_point, n_component))
fluid_zeros = np.zeros((n_point, n_pair))

fluid.beta = calc_beta(fluid.temperature)
fluid.epsilon = calc_epsilon(fluid.epsilon_r)
fluid.charge = calc_charge(fluid.valence)
fluid.charge_pair = calc_charge_pair(fluid.beta, fluid.charge, fluid.epsilon, n_component, n_pair)
fluid.rho = calc_rho(fluid.concentration)
kappa = calc_kappa(fluid.beta, fluid.charge, fluid.rho, fluid.epsilon)
st.sidebar.text(f"kappa: {kappa}")


r = z    # r on same discretisation as z
beta_u = fluid.beta * calc_u_lj(epsilon_lj, sigma_lj, n_point, n_component, n_pair, r)

beta_phiw = fluid.beta * calc_phiw(z, n_point, n_component)
beta_psi =  0.0
beta_psi_charge = -beta_psi * fluid.charge

model = Model(z=z, z_index=z_index, hw=wall_zeros, 
                c_short=fluid_zeros, f1=fluid_zeros, f2=fluid_zeros)

# Bulk-fluid inputs (direct correlation function
#TODO: Run bulk calculation on the fly (or from stored cr)

CR_PATH = f"{Path.cwd().as_posix()}/data/{fluid.cr_filename}"

st.sidebar.text(fluid.cr_filename)

# Run calculation

if run_calc := st.button("Run calculation"):
    try:
        c_short, r_short = load_and_interpolate_cr(Path(CR_PATH), n_point, n_pair, z)
    except Exception as e:
        st.error(f"Please ensure the c(r) data file is available here: {Path(CR_PATH).parent.as_posix()}")
        st.info("Restarting in 10 seconds")
        # st.exception(e)
        sleep(10.0)
        st.experimental_rerun()

    f1 = integral_z_infty_dr_r_c_short(c_short, n_pair, z, model.f1)
    f2 = integral_z_infty_dr_r2_c_short(c_short, n_pair, z, model.f2)

    f1_integrand = calc_f1_integrand(c_short, n_pair, z, n_point)
    f2_integrand = calc_f2_integrand(c_short, n_pair, z, n_point)

    # initial guess of zero - maybe should be \beta \phi

    tw_initial = np.zeros((n_point, n_component))
    hw_initial = calc_hw(tw_initial, n_component, beta_phiw)

    model.f1 = f1
    model.f2 = f2

    # Output to main page

    #TODO: Save/write solution & params to disk (on a "continuous" basis e.g. after every 10 iterations)
    #TODO: Handle numerical (e.g. Jacobian) exceptions gracefully
    #TODO: Plot |F(x)| convergence - pull values out of st_redirect

    col1, col2 = st.columns([1, 2])

    # Solve equation

    out_filepath = f"{Path.cwd()}/output/optim-solver-out.txt"


    with col1:
        with st.spinner("Finding optimal solution:"):
            st.markdown("Solver output (Newton-Krylov)")
            to_out = st.empty()
            with rd.stdout(to=to_out, format='text', max_buffer=20): # , rd.stdout(format='text', to=to_out_file):
                solution = solve_model(opt_func, tw_initial, fluid, model, d, 
                                            beta_phiw, beta_psi_charge)
        tw_solution = solution.x
        hw_solution = calc_hw(tw_solution, n_component, beta_phiw)
        st.write(solution)



    with col2:
        z_plots = dict({"hw_solution": dict({"fn_label": "hw", "plot_fn": hw_solution+1, 
            "plot_name": "Solution: gw(z)"})})

        plot_selected_wall_curves(n_component, z, z_plots)

        r_plots = dict({"c_short": dict({"fn_label": "c_short", "plot_fn": c_short, 
            "plot_name": "c_short(r)"})})
        r_plots["beta*u"] = dict({"fn_label": "beta u", "plot_fn": beta_u, "plot_name": "beta u(r)",
            "xlim": [0, 10], "ylim": [-10, 100]})
        r_plots["f1"] = dict({"fn_label": "f1", "plot_fn": f1, "plot_name": "f1(r)",  "xlim": [0, 10], "ylim": None})
        r_plots["f2"] = dict({"fn_label": "f2", "plot_fn": f2, "plot_name": "f2(r)",  "xlim": [0, 10], "ylim": None})

        r_plots["f1_integrand"] = dict({"fn_label": "f1_integrand", "plot_fn": f1_integrand, "plot_name": "f1_integrand(r)",  "xlim": [0, 10], "ylim": None})
        r_plots["f2_integrand"] = dict({"fn_label": "f2_integrand", "plot_fn": f2_integrand, "plot_name": "f2_integrand(r)",  "xlim": [0, 10], "ylim": None})

        plot_selected_bulk_curves(n_component, r, r_plots)
