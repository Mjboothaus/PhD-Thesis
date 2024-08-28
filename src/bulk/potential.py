import numpy as np
from scipy.special import erf
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
        U_discontinuity = np.zeros_like(U_ij)
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
        # Convert sigma to a NumPy array if it's not already
        sigma_array = np.array(self.parameters["sigma"])

        # Vectorized approach: Use broadcasting to calculate sigma for all pairs
        sigma = 0.5 * (sigma_array[:, None] + sigma_array[None, :])

        # Create a mask for r < sigma
        mask = r[None, None, :] < sigma[:, :, None]

        # Use the mask to set values
        U_ij_individual["hs"] = np.where(mask, np.inf, 0)
        U_discontinuity[:] = U_ij_individual["hs"]  # Use [:] to modify in-place

        U_ij += U_ij_individual["hs"]

    def _add_lennard_jones_potential(
        self,
        U_ij: np.ndarray,
        U_ij_individual: Dict[str, np.ndarray],
        dU_ij_individual: Dict[str, np.ndarray],
        r: np.ndarray,
    ):
        # Convert parameters to NumPy arrays
        epsilon_array = np.array(self.parameters["epsilon"])
        sigma_array = np.array(self.parameters["sigma"])

        # Vectorized approach: Use broadcasting to calculate epsilon and sigma for all pairs
        epsilon = np.sqrt(epsilon_array[:, None] * epsilon_array[None, :])
        sigma = 0.5 * (sigma_array[:, None] + sigma_array[None, :])

        # Calculate LJ potential for all pairs and distances at once
        sigma_r = sigma[:, :, None] / r[None, None, :]
        U_lj = 4 * epsilon[:, :, None] * (sigma_r**12 - sigma_r**6)
        dU_lj = 4 * epsilon[:, :, None] * (-12 * sigma_r**12 + 6 * sigma_r**6) / r[None, None, :]

        U_ij_individual["lj"] = U_lj
        dU_ij_individual["lj"] = dU_lj
        U_ij += U_lj

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
        # Convert charge to a NumPy array if it's not already
        charge_array = np.array(self.parameters["charge"])

        # Vectorized approach: Calculate q_ij for all pairs
        q_ij = charge_array[:, None] * charge_array[None, :]

        # Calculate Coulomb potential for all pairs and distances at once
        U_coulomb = q_ij[:, :, None] / r[None, None, :]
        dU_coulomb = -q_ij[:, :, None] / r[None, None, :] ** 2

        U_ij_individual["coulomb"] = U_coulomb
        dU_ij_individual["coulomb"] = dU_coulomb

        # Calculate erf-corrected potentials
        alpha = self.config["alpha"]
        U_erf = q_ij[:, :, None] * erf(alpha * r[None, None, :]) / r[None, None, :]
        U_erf_ij["real"] = U_erf

        U_erf_fourier = (
            4
            * np.pi
            * q_ij[:, :, None]
            * np.exp(-((k[None, None, :] / (2 * alpha)) ** 2))
            / k[None, None, :] ** 2
        )
        U_erf_ij["fourier"] = U_erf_fourier

        U_ij += U_coulomb

    def def_modMayerFunc(
        self,
        U_ij: np.ndarray,
        U_ij_individual: Dict[str, np.ndarray],
        U_discontinuity: np.ndarray,
        U_erf_ij_real: np.ndarray,
    ) -> Dict[str, np.ndarray]:
        # Vectorized approach: Calculate all modified Mayer functions at once
        modMayerFunc = {
            "u_ij": np.exp(-self.constants.beta * U_ij),
            "u_discontinuity": np.exp(-self.constants.beta * U_discontinuity),
            "u_erf": np.exp(self.constants.beta * U_erf_ij_real),
        }

        # Add individual potential contributions
        for potential_type in ["hs", "lj", "coulomb"]:
            if potential_type in U_ij_individual:
                modMayerFunc[f"u_{potential_type}"] = np.exp(
                    -self.constants.beta * U_ij_individual[potential_type]
                )

        return modMayerFunc


# Usage:
# potential_calculator = PotentialCalculator(config, system, parameters, constants)
# U_ij, U_ij_individual, dU_ij_individual, U_discontinuity, U_erf_ij = potential_calculator.calculate_potentials(dft, r, k)
# modMayerFunc = potential_calculator.def_modMayerFunc(U_ij, U_ij_individual, U_discontinuity, U_erf_ij['real'])
