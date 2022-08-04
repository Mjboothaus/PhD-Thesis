import sys
from pathlib import Path
from time import sleep

import pandas as pd
#import st_redirect as rd
import st_redirect_new as rd
import streamlit as st
from modelling import *
from numerics import set_num_parameters
from plotting import plot_bulk_curves, plot_wall_curves
from sidebar import create_sidebar


# Initialise fluid and numerical parameters

fluid_symbol = "kcl"

fluid = set_fluid_parameters(fluid_symbol)
if fluid is not None:
    z_cutoff, n_point, psi_0, tolerance, max_iteration = create_sidebar(fluid)
else:
    st.error("Invalid choice of fluid")

other_params = fluid_specific_parameters(fluid_symbol)

#TODO: Work out how to make the other parameters less custom
n_outer_shell = other_params["kcl"]["n_outer_shell"]
b = other_params["kcl"]["b"]
alpha = other_params["kcl"]["alpha"]
sigma = other_params["kcl"]["sigma"]
cap_c = other_params["kcl"]["cap_c"]
cap_d = other_params["kcl"]["cap_d"]

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

#TODO: Move these prelim calcs inside the potential (in modelling)
beta_pauling = calc_beta_pauling(fluid.valence, n_outer_shell, n_component, n_pair)
cap_b = calc_cap_b(beta_pauling, b, alpha, sigma, n_component, n_pair)

r = z    # r on same discretisation as z
beta_u = fluid.beta * calc_u(fluid.charge, cap_b, alpha, cap_c, cap_d,
           n_point, n_component, n_pair, fluid.epsilon, r)

beta_phiw = fluid.beta * calc_phiw(z, n_point, n_component)
beta_psi =  fluid.beta * psi_0 * 1.0e-3  # Convert to mV (in Volts)
beta_psi_charge = -beta_psi * fluid.charge

model = Model(z=z, z_index=z_index, hw=wall_zeros, 
                c_short=fluid_zeros, f1=fluid_zeros, f2=fluid_zeros)

# Bulk-fluid inputs (direct correlation function


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
    normal_stdout = sys.stdout
    
    to_out = st.empty()
    with col1:
        with st.spinner("Finding optimal solution:"):
            st.markdown("Solver output (Newton-Krylov)")

            #to_out_file = ""
            #out_file = open(out_filepath, 'w')
            with rd.stdout(to=to_out, format='markdown', max_buffer=20): # , rd.stdout(format='text', to=to_out_file):
                solution = solve_model(opt_func, tw_initial, fluid, model, d, 
                                            beta_phiw, beta_psi_charge)
        tw_solution = solution.x
        hw_solution = calc_hw(tw_solution, n_component, beta_phiw)
        st.write(solution)

    sys.stdout = normal_stdout

    with open(out_filepath, "w") as out_file:
        out_file.writelines(str(sys.stdout))


    with col2:
        z_plots = dict({"hw_solution": dict({"fn_label": "hw", "plot_fn": hw_solution+1, 
            "plot_name": "Solution: gw(z)"})})

        plot_wall_curves(n_component, z, z_plots)

        r_plots = dict({"c_short": dict({"fn_label": "c_short", "plot_fn": c_short, 
            "plot_name": "c_short(r)"})})
        r_plots["beta*u"] = dict({"fn_label": "beta u", "plot_fn": beta_u, "plot_name": "beta u(r)",
            "xlim": [0, 10], "ylim": [-100, 200]})
        r_plots["f1"] = dict({"fn_label": "f1", "plot_fn": f1, "plot_name": "f1(r)",  "xlim": [0, 10], "ylim": [-10, 10]})
        r_plots["f2"] = dict({"fn_label": "f2", "plot_fn": f2, "plot_name": "f2(r)",  "xlim": [0, 10], "ylim": [-10, 10]})

        r_plots["f1_integrand"] = dict({"fn_label": "f1_integrand", "plot_fn": f1_integrand, "plot_name": "f1_integrand(r)",  "xlim": [0, 10], "ylim": [-10, 10]})
        r_plots["f2_integrand"] = dict({"fn_label": "f2_integrand", "plot_fn": f2_integrand, "plot_name": "f2_integrand(r)",  "xlim": [0, 10], "ylim": [-10, 10]})

        plot_bulk_curves(n_component, r, r_plots)


        # beta_u, beta_phiw, c_short, f1, f2, f1_integrand, f2_integrand, hw_initial, hw_solution
