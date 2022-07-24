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

st.markdown(
    "Molten potassium chloride at a charged interface: Now in progress...")


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

z = discrete.z
r = discrete.z    # r on same discretisation as z

u = calc_u(charge, cap_b, alpha, cap_c, cap_d,
           n_point, n_component, n_pair, epsilon, r)
u = beta * u


rho = calc_rho(fluid.concentration)

kappa = calc_kappa(beta, charge, rho, epsilon)
st.sidebar.text(f"kappa: {kappa}")

phiw = calc_phiw(z, n_component, discrete.phiw)
phiw = beta * phiw

psi_0 = psi_0 * 1e-3     # 100 mV (in Volts)

# Bulk-fluid inputs (direct correlation function

cr_path = "/Users/mjboothaus/code/github/mjboothaus/PhD-Thesis/pyOZ_bulk_fluid/tests/lj/nrcg-cr.dat.orig"
c_short, r_short = load_and_intepolate_cr(Path(cr_path), n_point, n_pair, z)

# Output to main page

fig = plotly_line(r, u, ["r", "u0", "u1", "u2"], y_label="beta * u", legend_label="", 
        xliml=[0, 10], yliml=[-100, 200], title="Dimensionless ion-ion potential")
st.plotly_chart(fig)

fig = plotly_line(z, phiw, ["z", "phi0", "phi1"], y_label="beta * phiw", legend_label="", 
        xliml=[0, 10], yliml=[-20, 40], title="Dimensionless short-range wall-ion potential")
st.plotly_chart(fig)

fig = plotly_line(r, c_short, ["r", "c0", "c1", "c2"], y_label="c_ij(r)", legend_label="", 
        xliml=[0, 10], yliml=[-2, 2], title="'Short-range' direct correlation function")
st.plotly_chart(fig)