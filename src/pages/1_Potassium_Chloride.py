
from dacite import from_dict
from dataclasses import dataclass
import numpy as np
import streamlit as st
# from helper_functions import read_render_markdown_file
from model_parameters import get_parameter_variables


@dataclass
class Fluid:
    name: str
    component: list[str]
    valence: np.ndarray
    temperature: float
    concentration: np.ndarray
    epsilon_r: float
    n_component: int
    n_pair: int


def create_sidebar(fluid):
    st.sidebar.text(fluid.name)
    st.sidebar.text("Inputs:")
    st.sidebar.number_input(label=r"$\psi_0$ (mV):",
                            value=0, min_value=-1000, max_value=1000)
    z_cutoff = st.sidebar.number_input(
        "z_cutoff (A)", value=25, min_value=20, max_value=50)
    nPoint = st.sidebar.number_input("nPoint:", value=2001)
    grid_size = z_cutoff / (nPoint - 1)
    st.sidebar.text(f"Grid size (A): {grid_size}")
    return None


def set_fluid_parameters():
    parameters = dict({"name": "KCl", "component": ["K", "Cl"], "valence": np.array([
                    1.0, -1.0]), "temperature": 1075.0, "concentration": np.array([19.265, 19.265]),
        "epsilon_r": 1.0})

    n_component = parameters["n_component"] = len(parameters["component"])
    parameters["n_pair"] = int((n_component+1) * (n_component) / 2)
    return from_dict(data_class=Fluid, data=parameters)

fluid = set_fluid_parameters()

create_sidebar(fluid)
