"""Liquid water page - not yet fully implemented."""

import streamlit as st
from sidebar import create_sidebar

from charged_fluids import set_fluid_parameters

fluid_symbol = "h2o"
fluid = set_fluid_parameters(fluid_symbol)
if fluid is None:
    st.error("Invalid choice of fluid")
    st.stop()

_ = create_sidebar(fluid)
st.subheader(f"Charged fluids near an interface: {fluid.name}")
st.warning(
    "Liquid water is **not implemented** yet: bulk c(r) and the molecular model "
    "are still outstanding (see `data/fluid_parameters.toml` status)."
)
st.info("Status from TOML: **" + fluid.status + "**")
