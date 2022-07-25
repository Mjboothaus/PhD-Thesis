from ctypes import c_short
from unicodedata import numeric
import streamlit as st
from model_parameters import *
from numerical_parameters import create_sidebar, set_num_parameters
from plotting import make_simple_plot, plotly_line
import pandas as pd
from pathlib import Path

# from helper_functions import read_render_markdown_file

fluid_symbol = "kcl"

fluid = set_fluid_parameters(fluid_symbol)


if fluid is not None:
    z_cutoff, n_point, psi_0 = create_sidebar(fluid)
else:
    st.error("Invalid choice of fluid")


other_params = fluid_specific_parameters(fluid_symbol)

n_outer_shell = other_params["kcl"]["n_outer_shell"]
b = other_params["kcl"]["b"]
alpha = other_params["kcl"]["alpha"]
sigma = other_params["kcl"]["sigma"]
cap_c = other_params["kcl"]["cap_c"]
cap_d = other_params["kcl"]["cap_d"]

discrete = set_num_parameters(
    n_point, z_cutoff, fluid.n_component, fluid.n_pair)

n_component = discrete.n_component
n_pair = discrete.n_pair

beta = calc_beta(fluid.temperature)
epsilon = calc_epsilon(fluid.epsilon_r)

beta_pauling = calc_beta_pauling(
    fluid.valence, n_outer_shell, n_component, n_pair)

cap_b = calc_cap_b(beta_pauling, b, alpha, sigma, n_component, n_pair)

charge = calc_charge(fluid.valence)
charge_pair = calc_charge_pair(beta, charge, epsilon, n_component, n_pair)

z = discrete.z
r = discrete.z    # r on same discretisation as z

beta_u = beta * calc_u(charge, cap_b, alpha, cap_c, cap_d,
           n_point, n_component, n_pair, epsilon, r)

rho = calc_rho(fluid.concentration)

kappa = calc_kappa(beta, charge, rho, epsilon)
st.sidebar.text(f"kappa: {kappa}")

beta_phiw = beta * calc_phiw(z, n_component, discrete.phiw)
beta_psi =  beta * psi_0 * 1.0e-3  # 100 mV (in Volts)
beta_psi_charge = -beta_psi * charge

# Bulk-fluid inputs (direct correlation function

# cr_path = "/Users/mjboothaus/code/github/mjboothaus/PhD-Thesis/pyOZ_bulk_fluid/tests/lj/nrcg-cr.dat.orig"

cr_path = "data/pyOZ_bulk_fluid/tests/lj/nrcg-cr.dat.orig"
c_short, r_short = load_and_intepolate_cr(Path(cr_path), n_point, n_pair, z)


f1 = integral_z_infty_dr_r_c_short(c_short, n_pair, z, discrete.f1)
f2 = integral_z_infty_dr_r2_c_short(c_short, n_pair, z, discrete.f2)

# initial guess of zero - maybe should be \beta \phi
tw_initial = np.zeros((n_point, n_component))

hw_initial = calc_hw(tw_initial, n_component, beta_phiw, discrete.hw)


tw = calc_tw(tw_initial, discrete.hw, discrete.tw, beta_phiw, beta_psi_charge, charge_pair, 
            rho, f1, f2, z, n_component, n_point, discrete.z_index, 
            discrete.integral_z_infty, discrete.integral_0_z)





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