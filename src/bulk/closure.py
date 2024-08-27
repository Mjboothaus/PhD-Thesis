# File: pyoz_closure.py

import numpy as np
from typing import Dict, Tuple, Any


class ClosureCalculator:
    def __init__(self, system: Dict[str, Any]):
        self.system = system

    def calculate_gamma_term(
        self,
        r: np.ndarray,
        mod_mayer_func: Dict[str, np.ndarray],
        u_discontinuity: np.ndarray,
        g_r_ij: np.ndarray,
    ) -> Tuple[np.ndarray, np.ndarray]:
        cs_r_ij = np.zeros_like(g_r_ij)

        for i in range(self.system["ncomponents"]):
            for j in range(self.system["ncomponents"]):
                cs_r_ij[i, j] = self._calculate_component_gamma(
                    r, mod_mayer_func, u_discontinuity, g_r_ij, i, j
                )

        return cs_r_ij, g_r_ij

    def _calculate_component_gamma(
        self,
        r: np.ndarray,
        mod_mayer_func: Dict[str, np.ndarray],
        u_discontinuity: np.ndarray,
        g_r_ij: np.ndarray,
        i: int,
        j: int,
    ) -> np.ndarray:
        if self.system["closure"] == "hnc":
            return self._hnc_closure(mod_mayer_func, g_r_ij, i, j)
        elif self.system["closure"] == "py":
            return self._py_closure(mod_mayer_func, g_r_ij, i, j)
        elif self.system["closure"] == "msa":
            return self._msa_closure(r, mod_mayer_func, u_discontinuity, g_r_ij, i, j)
        else:
            raise ValueError(f"Unknown closure: {self.system['closure']}")

    def _hnc_closure(
        self, mod_mayer_func: Dict[str, np.ndarray], g_r_ij: np.ndarray, i: int, j: int
    ) -> np.ndarray:
        """
        Hypernetted-chain (HNC) closure
        $$c(r) = e^{-\beta u(r) + \gamma(r)} - 1 - \gamma(r)$$
        """
        return (mod_mayer_func["u_ij"][i, j] - 1.0) * mod_mayer_func["u_erf"][i, j] + np.log(
            g_r_ij[i, j] * mod_mayer_func["u_erf"][i, j]
        )

    def _py_closure(
        self, mod_mayer_func: Dict[str, np.ndarray], g_r_ij: np.ndarray, i: int, j: int
    ) -> np.ndarray:
        """
        Percus-Yevick (PY) closure
        $$c(r) = (1 - e^{\beta u(r)})(1 + \gamma(r))$$
        """
        return (mod_mayer_func["u_ij"][i, j] - 1.0) * mod_mayer_func["u_erf"][i, j] * g_r_ij[i, j]

    def _msa_closure(
        self,
        r: np.ndarray,
        mod_mayer_func: Dict[str, np.ndarray],
        u_discontinuity: np.ndarray,
        g_r_ij: np.ndarray,
        i: int,
        j: int,
    ) -> np.ndarray:
        """
        Mean Spherical Approximation (MSA) closure
        $$c(r) = -\beta u(r) \quad \text{for } r > \sigma$$
        $$g(r) = 0 \quad \text{for } r < \sigma$$
        where $\sigma$ is the particle diameter.
        """
        cs_r = np.zeros_like(r)
        mask = r < self.system["sigma"][i, j]
        cs_r[mask] = -u_discontinuity[i, j, mask]
        cs_r[~mask] = (mod_mayer_func["u_ij"][i, j, ~mask] - 1.0) * mod_mayer_func["u_erf"][i, j, ~mask]
        return cs_r


# Usage:
# closure_calculator = ClosureCalculator(system)
# cs_r_ij, g_r_ij = closure_calculator.calculate_gamma_term(r, mod_mayer_func, u_discontinuity, g_r_ij)
