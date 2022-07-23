
from dataclasses import dataclass
from scipy.constants import epsilon_0, elementary_charge, Boltzmann, Avogadro
import numpy as np
from dacite import from_dict

# from helper_functions import read_render_markdown_file


@dataclass
class Fluid:
    name: str
    symbol: str
    component: list[str]
    valence: np.ndarray
    temperature: float
    concentration: np.ndarray
    epsilon_r: float
    n_component: int
    n_pair: int
    index: int


fluid_parameters = dict({"kcl": dict({"name": "Potassium Chloride", "component": ["K", "Cl"], "valence": np.array([
    1.0, -1.0]), "temperature": 1075.0, "concentration": np.array([19.265, 19.265]),
    "epsilon_r": 1.0, "index": 1})})

fluid_parameters["h2o"] = dict({"name": "Liquid water", "component": ["H", "2O"], "valence": np.array([
    1.0, -1.0]), "temperature": 298.0, "concentration": np.array([1.0, 1.0]),
    "epsilon_r": 1.0, "index": 2})

fluid_parameters["2_2"] = dict({"name": "2-2 Aqueous electrolyte", "component": ["+2", "-2"], "valence": np.array([
    2.0, -2.0]), "temperature": 298.0, "concentration": np.array([1.0, 1.0]),
    "epsilon_r": 78.3, "index": 3})


def set_fluid_parameters(symbol):
    symbol = symbol.lower()
    if symbol not in fluid_parameters:
        return None
    parameters = fluid_parameters[symbol]
    parameters["symbol"] = "".join(parameters["component"])
    n_component = parameters["n_component"] = len(parameters["component"])
    parameters["n_pair"] = int((n_component+1) * (n_component) / 2)
    return from_dict(data_class=Fluid, data=parameters)


def fluid_specific_parameters(symbol):
    if symbol == "kcl":
        other_params = dict({"kcl": dict({"n_outer_shell": np.array([8., 8.])})})
        
    return other_params


def calc_beta(temperature):
    return 1.0 / (Boltzmann * temperature)


def calc_epsilon(epsilon_r):
     # units same as $\epsilon_0$ (need to allow for distances in angstrom)
    return 4.0 * np.pi * epsilon_r * epsilon_0 


def calc_l_index(i, j):
    return i + j    


def calc_beta_pauling(valence, n_outer_shell, n_component, n_pair):
    beta_pauling = np.zeros(n_pair)
    for i in range(n_component):
        for j in range(i, n_component):
            l = calc_l_index(i, j)
            beta_pauling[l] = 1.0 + valence[i] / n_outer_shell[i] + valence[j] / n_outer_shell[j]
    return beta_pauling