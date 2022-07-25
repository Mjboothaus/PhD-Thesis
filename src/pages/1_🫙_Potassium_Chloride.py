import streamlit as st
from model_parameters import *
from numerical_parameters import create_sidebar, set_num_parameters
from plotting import plotly_line
import pandas as pd
from pathlib import Path

# from helper_functions import read_render_markdown_file

fluid_symbol = "kcl"

fluid = set_fluid_parameters(fluid_symbol)

st.write(fluid.charge_pair)

if fluid is not None:
    z_cutoff, n_point, psi_0, tolerance, max_iteration = create_sidebar(fluid)
else:
    st.error("Invalid choice of fluid")


other_params = fluid_specific_parameters(fluid_symbol)

n_outer_shell = other_params["kcl"]["n_outer_shell"]
b = other_params["kcl"]["b"]
alpha = other_params["kcl"]["alpha"]
sigma = other_params["kcl"]["sigma"]
cap_c = other_params["kcl"]["cap_c"]
cap_d = other_params["kcl"]["cap_d"]

d = set_num_parameters(n_point, z_cutoff, fluid.n_component, fluid.n_pair)

n_component = d.n_component
n_pair = d.n_pair

fluid.beta = calc_beta(fluid.temperature)
fluid.epsilon = calc_epsilon(fluid.epsilon_r)

beta_pauling = calc_beta_pauling(
    fluid.valence, n_outer_shell, n_component, n_pair)

cap_b = calc_cap_b(beta_pauling, b, alpha, sigma, n_component, n_pair)

charge = calc_charge(fluid.valence)
fluid.charge_pair = calc_charge_pair(fluid.beta, charge, fluid.epsilon, n_component, n_pair)

fluid.rho = calc_rho(fluid.concentration)

z = np.linspace(0.0, z_cutoff, n_point)
z_index = np.arange(0, n_point, dtype=int)
wall_zeros = np.zeros((n_point, n_component))
fluid_zeros = np.zeros((n_point, n_pair))

r = z    # r on same discretisation as z
beta_u = fluid.beta * calc_u(charge, cap_b, alpha, cap_c, cap_d,
           n_point, n_component, n_pair, fluid.epsilon, r)

model = Model(z=z, z_index=z_index, hw=wall_zeros, tw=wall_zeros, phiw=wall_zeros,
    c_short=fluid_zeros, f1=fluid_zeros, f2=fluid_zeros,
    integral_0_z=wall_zeros, integral_z_infty=wall_zeros)






kappa = calc_kappa(fluid.beta, charge, fluid.rho, fluid.epsilon)
st.sidebar.text(f"kappa: {kappa}")

beta_phiw = fluid.beta * calc_phiw(z, n_component)
beta_psi =  fluid.beta * psi_0 * 1.0e-3  # 100 mV (in Volts)
beta_psi_charge = -beta_psi * charge

# Bulk-fluid inputs (direct correlation function

# cr_path = "/Users/mjboothaus/code/github/mjboothaus/PhD-Thesis/pyOZ_bulk_fluid/tests/lj/nrcg-cr.dat.orig"

cr_path = "data/pyOZ_bulk_fluid/tests/lj/nrcg-cr.dat.orig"
c_short, r_short = load_and_interpolate_cr(Path(cr_path), n_point, n_pair, z)


f1 = integral_z_infty_dr_r_c_short(c_short, n_pair, z, d.f1)
f2 = integral_z_infty_dr_r2_c_short(c_short, n_pair, z, d.f2)

# initial guess of zero - maybe should be \beta \phi

tw_initial = np.zeros((n_point, n_component))
hw_initial = calc_hw(tw_initial, n_component, beta_phiw)


# tw = calc_tw(tw_initial, fluid, model, d)
# hw = calc_hw(tw, n_component, beta_phiw)

# Solve equation

solution = solve_model(opt_func, tw_initial, tolerance, max_iteration)




# Output to main page




fig = plotly_line(r, beta_u, ["r", "u0", "u1", "u2"], y_label="beta * u", legend_label="",
                  xliml=[0, 10], yliml=[-100, 200], title="Dimensionless ion-ion potential")
st.plotly_chart(fig)


fig = plotly_line(z, beta_phiw, ["z", "phi0", "phi1"], y_label="beta * phiw", legend_label="",
                  xliml=[0, 10], yliml=[-20, 40], title="Dimensionless short-range wall-ion potential")
st.plotly_chart(fig)


fig = plotly_line(r, c_short, ["r", "c0", "c1", "c2"], y_label="c_ij(r)", legend_label="",
                  xliml=[0, 10], yliml=[-2, 2], title="'Short-range' direct correlation function")
st.plotly_chart(fig)


fig = plotly_line(z, f1, ["z", "f1_0", "f1_1", "f1_2"], y_label="f1_ij(r)", legend_label="",
                  xliml=[0, 10], yliml=[-15, 5], title="f1 function")
st.plotly_chart(fig)


fig = plotly_line(z, f2, ["z", "f2_0", "f2_1", "f2_2"], y_label="f2_ij(r)", legend_label="",
                  xliml=[0, 10], yliml=[-40, 5], title="f2 function")
st.plotly_chart(fig)


fig = plotly_line(z, hw_initial, ["z", "hw0", "hw1"], y_label="hw", legend_label="",
                  xliml=[0, 10], yliml=[-2, 2], title="hw_initial for tw = tw_initial (zero guess)")
st.plotly_chart(fig)


fig = plotly_line(z, tw, ["z", "tw0", "tw1"], y_label="tw", legend_label="",
                  xliml=[0, 10], title="tw")
st.plotly_chart(fig)


fig = plotly_line(z, hw, ["z", "hw0", "hw1"], y_label="hw", legend_label="",
                  xliml=[0, 10], title="hw")
st.plotly_chart(fig)