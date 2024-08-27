# File: pyoz_property.py

import numpy as np
from typing import Dict, Any


class PropertyCalculator:
    def __init__(self, config: Dict[str, Any], system: Dict[str, Any], constants: Any):
        self.config = config
        self.system = system
        self.constants = constants

    def calculate_properties(
        self, r: np.ndarray, g_r_ij: np.ndarray, S_k_ij: np.ndarray
    ) -> Dict[str, float]:
        properties = {}
        properties.update(self._calculate_thermodynamic_properties(r, g_r_ij))
        properties.update(self._calculate_structural_properties(S_k_ij))
        return properties

    def _calculate_thermodynamic_properties(self, r: np.ndarray, g_r_ij: np.ndarray) -> Dict[str, float]:
        properties = {}
        properties["pressure"] = self._calculate_pressure(r, g_r_ij)
        properties["internal_energy"] = self._calculate_internal_energy(r, g_r_ij)
        properties["isothermal_compressibility"] = self._calculate_isothermal_compressibility(g_r_ij)
        return properties

    def _calculate_structural_properties(self, S_k_ij: np.ndarray) -> Dict[str, float]:
        properties = {}
        properties["S_0"] = self._calculate_structure_factor_at_zero(S_k_ij)
        properties["correlation_length"] = self._calculate_correlation_length(S_k_ij)
        return properties

    def _calculate_pressure(self, r: np.ndarray, g_r_ij: np.ndarray) -> float:
        """
        Calculate pressure using the virial equation:
        $$P = \rho k_B T - \frac{2\pi\rho^2}{3} \int_0^\infty r^3 \frac{dU(r)}{dr} g(r) dr$$
        """
        rho = self.system["density"]
        kT = self.constants.kT
        dr = r[1] - r[0]

        integral = np.sum(r[1:] ** 3 * (g_r_ij[1:] - 1) * dr)
        return rho * kT * (1 - 2 * np.pi * rho * integral / 3)

    def _calculate_internal_energy(self, r: np.ndarray, g_r_ij: np.ndarray) -> float:
        """
        Calculate internal energy:
        $$U = \frac{3}{2}Nk_BT + 2\pi\rho N \int_0^\infty r^2 U(r) g(r) dr$$
        """
        N = self.system["number_of_particles"]
        rho = self.system["density"]
        kT = self.constants.kT
        dr = r[1] - r[0]

        U_r = self.system["potential"](r)  # Assuming a method to get potential
        integral = np.sum(r**2 * U_r * g_r_ij * dr)
        return 1.5 * N * kT + 2 * np.pi * rho * N * integral

    def _calculate_isothermal_compressibility(self, g_r_ij: np.ndarray) -> float:
        """
        Calculate isothermal compressibility:
        $$\kappa_T = \frac{1}{\rho k_B T} \left(1 + 4\pi\rho \int_0^\infty r^2 [g(r) - 1] dr\right)$$
        """
        rho = self.system["density"]
        kT = self.constants.kT
        dr = self.config["deltar"]

        integral = np.sum((g_r_ij - 1) * dr)
        return (1 + 4 * np.pi * rho * integral) / (rho * kT)

    def _calculate_structure_factor_at_zero(self, S_k_ij: np.ndarray) -> float:
        """
        Calculate structure factor at k=0:
        $$S(0) = \lim_{k \to 0} S(k)$$
        """
        return S_k_ij[0]

    def _calculate_correlation_length(self, S_k_ij: np.ndarray) -> float:
        """
        Calculate correlation length:
        xi = \sqrt{\frac{1}{6\pi^2\rho} \int_0^\infty k^2 [S(k) - 1] dk}
        """
        rho = self.system["density"]
        dk = self.config["deltak"]
        k = np.arange(len(S_k_ij)) * dk

        integral = np.sum(k**2 * (S_k_ij - 1) * dk)
        return np.sqrt(integral / (6 * np.pi**2 * rho))


# Usage:
# property_calculator = PropertyCalculator(config, system, constants)
# properties = property_calculator.calculate_properties(r, g_r_ij, S_k_ij)
