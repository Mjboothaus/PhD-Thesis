import streamlit as st
from model_parameters import *
from numerical_parameters import create_sidebar, set_num_parameters
from plotting import make_simple_plot, plot_plotly_line
import pandas as pd

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

r = discrete.z

u = calc_u(charge, cap_b, alpha, cap_c, cap_d,
           n_point, n_component, n_pair, epsilon, r)


# fig = make_simple_plot(r, beta * u, "r", "u", "Dimensionless ion-ion potential", xliml=[0, 10], yliml=[-100, 200])
# st.pyplot(fig)

fig = plot_plotly_line(r, beta*u, ["r", "u0", "u1", "u2"], y_label="beta * u", legend_label="", 
        xliml=[0, 10], yliml=[-100, 200], title="Dimensionless ion-ion potential")
st.plotly_chart(fig)


psi_0 = psi_0 * 1e-3     # 100 mV (in Volts)
