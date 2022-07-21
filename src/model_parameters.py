
import numpy as np
import pandas as pd

parameters = dict({"kcl": dict({"component": ["K", "Cl"], "valence": np.array([
                  1.0, -1.0]), "temperature": 1075.0, "concentration": [19.265, 19.265],
    "eps_r": 1.0})})

#TODO: Generalise for other fluids

def get_parameter_variables(fluid_name):
    for param_name in parameters[fluid_name].keys():
        locals().__setitem__(param_name, parameters[fluid_name][param_name])
        print(
            f"Created variable - {param_name}: {str(parameters[fluid_name][param_name])}")
            
    n_component = len(component)
    n_pair = int((n_component+1) * (n_component) / 2)
    return component, valence, temperature, concentration, eps_r, n_component, n_pair
