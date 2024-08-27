# File: pyoz_potential.py

import numpy as np
from typing import Dict, Tuple, Any


class PotentialCalculator:
    def __init__(
        self, config: Dict[str, Any], system: Dict[str, Any], parameters: Dict[str, Any], constants: Any
    ):
        self.config = config
        self.system = system
        self.parameters = parameters
        self.constants = constants

    def calculate_potentials(
        self, dft: Any, r: np.ndarray, k: np.ndarray
    ) -> Tuple[
        np.ndarray, Dict[str, np.ndarray], Dict[str, np.ndarray], np.ndarray, Dict[str, np.ndarray]
    ]:
        U_ij = np.zeros((self.system["ncomponents"], self.system["ncomponents"], self.config["npoints"]))
        U_ij_individual = {}
        dU_ij_individual = {}
        U_discontinuity = np.zeros(
            (self.system["ncomponents"], self.system["ncomponents"], self.config["npoints"])
        )
        U_erf_ij = {"real": np.zeros_like(U_ij), "fourier": np.zeros_like(U_ij)}

        for potential_type in self.system["potential"]:
            if potential_type == "hs":
                self._add_hard_sphere_potential(
                    U_ij, U_ij_individual, dU_ij_individual, U_discontinuity, r
                )
            elif potential_type == "lj":
                self._add_lennard_jones_potential(U_ij, U_ij_individual, dU_ij_individual, r)
            elif potential_type == "coulomb":
                self._add_coulomb_potential(U_ij, U_ij_individual, dU_ij_individual, U_erf_ij, dft, r, k)
            # Add other potential types here

        return U_ij, U_ij_individual, dU_ij_individual, U_discontinuity, U_erf_ij

    def _add_hard_sphere_potential(
        self,
        U_ij: np.ndarray,
        U_ij_individual: Dict[str, np.ndarray],
        dU_ij_individual: Dict[str, np.ndarray],
        U_discontinuity: np.ndarray,
        r: np.ndarray,
    ):
        U_ij_individual["hs"] = np.zeros_like(U_ij)
        dU_ij_individual["hs"] = np.zeros_like(U_ij)

        for i in range(self.system["ncomponents"]):
            for j in range(i, self.system["ncomponents"]):
                sigma = 0.5 * (self.parameters["sigma"][i] + self.parameters["sigma"][j])
                U_ij_individual["hs"][i, j] = np.where(r < sigma, np.inf, 0)
                U_ij_individual["hs"][j, i] = U_ij_individual["hs"][i, j]
                U_discontinuity[i, j] = np.where(r < sigma, np.inf, 0)
                U_discontinuity[j, i] = U_discontinuity[i, j]

        U_ij += U_ij_individual["hs"]

    def _add_lennard_jones_potential(
        self,
        U_ij: np.ndarray,
        U_ij_individual: Dict[str, np.ndarray],
        dU_ij_individual: Dict[str, np.ndarray],
        r: np.ndarray,
    ):
        U_ij_individual["lj"] = np.zeros_like(U_ij)
        dU_ij_individual["lj"] = np.zeros_like(U_ij)

        for i in range(self.system["ncomponents"]):
            for j in range(i, self.system["ncomponents"]):
                epsilon = np.sqrt(self.parameters["epsilon"][i] * self.parameters["epsilon"][j])
                sigma = 0.5 * (self.parameters["sigma"][i] + self.parameters["sigma"][j])

                U_lj = 4 * epsilon * ((sigma / r) ** 12 - (sigma / r) ** 6)
                U_ij_individual["lj"][i, j] = U_lj
                U_ij_individual["lj"][j, i] = U_lj

                dU_lj = 4 * epsilon * (-12 * (sigma**12 / r**13) + 6 * (sigma**6 / r**7))
                dU_ij_individual["lj"][i, j] = dU_lj
                dU_ij_individual["lj"][j, i] = dU_lj

        U_ij += U_ij_individual["lj"]

    def _add_coulomb_potential(
        self,
        U_ij: np.ndarray,
        U_ij_individual: Dict[str, np.ndarray],
        dU_ij_individual: Dict[str, np.ndarray],
        U_erf_ij: Dict[str, np.ndarray],
        dft: Any,
        r: np.ndarray,
        k: np.ndarray,
    ):
        U_ij_individual["coulomb"] = np.zeros_like(U_ij)
        dU_ij_individual["coulomb"] = np.zeros_like(U_ij)

        for i in range(self.system["ncomponents"]):
            for j in range(i, self.system["ncomponents"]):
                q_ij = self.parameters["charge"][i] * self.parameters["charge"][j]
                U_coulomb = q_ij / r
                U_ij_individual["coulomb"][i, j] = U_coulomb
                U_ij_individual["coulomb"][j, i] = U_coulomb

                dU_coulomb = -q_ij / r**2
                dU_ij_individual["coulomb"][i, j] = dU_coulomb
                dU_ij_individual["coulomb"][j, i] = dU_coulomb

                U_erf = q_ij * np.erf(self.config["alpha"] * r) / r
                U_erf_ij["real"][i, j] = U_erf
                U_erf_ij["real"][j, i] = U_erf

                U_erf_fourier = (
                    4 * np.pi * q_ij * np.exp(-((k / (2 * self.config["alpha"])) ** 2)) / k**2
                )
                U_erf_ij["fourier"][i, j] = U_erf_fourier
                U_erf_ij["fourier"][j, i] = U_erf_fourier

        U_ij += U_ij_individual["coulomb"]

    def def_modMayerFunc(
        self,
        U_ij: np.ndarray,
        U_ij_individual: Dict[str, np.ndarray],
        U_discontinuity: np.ndarray,
        U_erf_ij_real: np.ndarray,
    ) -> Dict[str, np.ndarray]:
        modMayerFunc = {
            "u_ij": np.exp(-self.constants.beta * U_ij),
            "u_hs": np.exp(-self.constants.beta * U_ij_individual["hs"]),
            "u_lj": np.exp(-self.constants.beta * U_ij_individual["lj"]),
            "u_coulomb": np.exp(-self.constants.beta * U_ij_individual["coulomb"]),
            "u_discontinuity": np.exp(-self.constants.beta * U_discontinuity),
            "u_erf": np.exp(self.constants.beta * U_erf_ij_real),
        }
        return modMayerFunc


# Usage:
# potential_calculator = PotentialCalculator(config, system, parameters, constants)
# U_ij, U_ij_individual, dU_ij_individual, U_discontinuity, U_erf_ij = potential_calculator.calculate_potentials(dft, r, k)
# modMayerFunc = potential_calculator.def_modMayerFunc(U_ij, U_ij_individual, U_discontinuity, U_erf_ij['real'])
