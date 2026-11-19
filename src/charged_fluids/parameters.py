"""Fluid parameter definitions and factories.

Note: charge and rho placeholders are filled by callers after construction.
A later phase will load a single TOML source of truth instead of the in-module
dict used here (characterisation-compatible with the historical app).
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass

import numpy as np


@dataclass
class Fluid:
    name: str
    symbol: str
    component: list[str]
    valence: np.ndarray
    charge: np.ndarray
    temperature: float
    concentration: np.ndarray
    epsilon_r: float
    n_component: int
    n_pair: int
    index: int
    charge_pair: np.ndarray
    rho: np.ndarray
    beta: float
    epsilon: float
    cr_filename: str = ""


# Placeholder arrays keep typing stable until derived fields are filled.
_fluid_parameters: dict[str, dict] = {
    "kcl": {
        "name": "Potassium Chloride",
        "component": ["K", "Cl"],
        "valence": np.array([1.0, -1.0]),
        "charge": np.array([0.0]),
        "temperature": 1075.0,
        "concentration": np.array([19.5, 19.5]),
        "epsilon_r": 1.0,
        "index": 1,
        "charge_pair": np.array([0.0]),
        "rho": np.array([0.0]),
        "beta": 0.0,
        "epsilon": 0.0,
        "cr_filename": "nrcg-cr-approx-maybe.dat",
    },
    "lj1": {
        "name": "Lennard-Jones liquid (1-comp)",
        "component": ["LJ"],
        "valence": np.array([0.0]),
        "charge": np.array([0.0]),
        "temperature": 298.15,
        "concentration": np.array([0.5]),
        "epsilon_r": 1.0,
        "index": 2,
        "charge_pair": np.array([0.0]),
        "rho": np.array([0.0]),
        "beta": 0.0,
        "epsilon": 0.0,
        "cr_filename": "lj1-cr.data",
    },
    "lj2": {
        "name": "Lennard-Jones liquid (2-comp)",
        "component": ["1", "2"],
        "valence": np.array([0.0, 0.0]),
        "charge": np.array([0.0]),
        "temperature": 298.15,
        "concentration": np.array([0.5, 0.5]),
        "epsilon_r": 1.0,
        "index": 5,
        "charge_pair": np.array([0.0]),
        "rho": np.array([0.0]),
        "beta": 0.0,
        "epsilon": 0.0,
        "cr_filename": "pyoz-cr-lj-equal-2-comp.data",
    },
    "h2o": {
        "name": "Liquid water",
        "component": ["H", "2O"],
        "valence": np.array([1.0, -1.0]),
        "charge": np.array([0.0]),
        "temperature": 298.0,
        "concentration": np.array([1.0, 1.0]),
        "epsilon_r": 1.0,
        "index": 2,
        "charge_pair": np.array([0.0]),
        "rho": np.array([0.0]),
        "beta": 0.0,
        "epsilon": 0.0,
        "cr_filename": "TO_BE_DEFINED.dat",
    },
    "2_2": {
        "name": "2-2 Aqueous electrolyte",
        "component": ["+2", "-2"],
        "valence": np.array([2.0, -2.0]),
        "charge": np.array([0.0]),
        "temperature": 298.15,
        "concentration": np.array([1.0, 1.0]),
        "epsilon_r": 78.3,
        "index": 5,
        "charge_pair": np.array([0.0]),
        "rho": np.array([0.0]),
        "beta": 0.0,
        "epsilon": 0.0,
        "cr_filename": "TO_BE_DEFINED.dat",
    },
}

# Public alias retained for callers that inspect the table directly.
fluid_parameters = _fluid_parameters


def set_fluid_parameters(symbol: str) -> Fluid | None:
    """Return a :class:`Fluid` for ``symbol``, or ``None`` if unknown."""
    symbol = symbol.lower()
    if symbol not in _fluid_parameters:
        return None

    parameters = deepcopy(_fluid_parameters[symbol])
    parameters["symbol"] = symbol
    n_component = len(parameters["component"])
    parameters["n_component"] = n_component
    parameters["n_pair"] = int((n_component + 1) * n_component / 2)

    return Fluid(
        name=parameters["name"],
        symbol=parameters["symbol"],
        component=parameters["component"],
        valence=parameters["valence"],
        charge=parameters["charge"],
        temperature=parameters["temperature"],
        concentration=parameters["concentration"],
        epsilon_r=parameters["epsilon_r"],
        n_component=parameters["n_component"],
        n_pair=parameters["n_pair"],
        index=parameters["index"],
        charge_pair=parameters["charge_pair"],
        rho=parameters["rho"],
        beta=parameters["beta"],
        epsilon=parameters["epsilon"],
        cr_filename=parameters["cr_filename"],
    )


def fluid_specific_parameters(symbol: str) -> dict:
    """Potential-model extras keyed by fluid symbol.

    Raises
    ------
    KeyError
        If the symbol has no fluid-specific table (e.g. incomplete water/2-2).
    """
    symbol = symbol.lower()
    if symbol == "kcl":
        return {
            "kcl": {
                "n_outer_shell": np.array([8.0, 8.0]),
                "alpha": 1.0 / 0.337,
                "b": 0.338e-19,
                "sigma": [1.463, 1.585],
                "cap_c": np.array([24.3, 48.0, 124.5]) * 1e-19,
                "cap_d": np.array([24.0, 73.0, 250.0]) * 1e-19,
            }
        }
    if symbol == "lj1":
        return {"lj1": {"epsilon_lj": np.array([0.1]), "sigma_lj": np.array([0.5])}}
    if symbol == "lj2":
        return {
            "lj2": {
                "epsilon_lj": np.array([0.1, 0.1, 0.1]),
                "sigma_lj": np.array([0.5, 0.5, 0.5]),
            }
        }
    raise KeyError(f"No fluid-specific parameters defined for symbol={symbol!r}")
