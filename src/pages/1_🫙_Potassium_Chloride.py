import streamlit as st
from model_parameters import Fluid, calc_beta, calc_epsilon, set_fluid_parameters, fluid_specific_parameters
from numerical_parameters import create_sidebar, set_num_parameters

# from helper_functions import read_render_markdown_file

fluid_symbol = "kcl"

fluid = set_fluid_parameters(fluid_symbol)

if fluid is not None:
    z_cutoff, n_point, psi_0 = create_sidebar(fluid)
    st.write(fluid)
else:
    st.error("Invalid choice of fluid")

st.markdown("Molten potassium chloride at a charged interface: Now in progress...")


other_params = fluid_specific_parameters(fluid_symbol)

# e.g. other_params["kcl"]["n_outer_shell"]

discrete = set_num_parameters(n_point, z_cutoff, fluid.n_pair)

beta = calc_beta(fluid.temperature)

epsilon = calc_epsilon(fluid.epsilon_r)



