# File: pyoz_solver.py

import numpy as np
from typing import Callable, Union


class Solver:
    def __init__(self, ncomponents: int):
        self.ncomponents = ncomponents
        self.solver_function = self._select_solver()

    def _select_solver(self) -> Callable:
        if self.ncomponents == 1:
            return self._solver_1
        elif self.ncomponents == 2:
            return self._solver_2
        else:
            return self._solver_n

    def solve(self, matrix: np.ndarray, vector: np.ndarray, npoints: int) -> np.ndarray:
        return self.solver_function(matrix, vector, npoints)

    @staticmethod
    def _solver_1(matrix: np.ndarray, vector: np.ndarray, npoints: int) -> np.ndarray:
        """
        Optimized solver for 1 component
        """
        return vector / (1 - matrix * vector)

    @staticmethod
    def _solver_2(matrix: np.ndarray, vector: np.ndarray, npoints: int) -> np.ndarray:
        """
        Optimized solver for 2 components
        """
        det = (
            1
            - (matrix[0, 0] * vector[0, 0] + matrix[1, 1] * vector[1, 1])
            + (matrix[0, 0] * matrix[1, 1] - matrix[0, 1] * matrix[1, 0]) * vector[0, 0] * vector[1, 1]
        )

        result = np.empty_like(vector)
        result[0, 0] = (
            vector[0, 0] * (1 - matrix[1, 1] * vector[1, 1]) + matrix[0, 1] * vector[0, 1] * vector[1, 1]
        ) / det
        result[0, 1] = (
            vector[0, 1] * (1 - matrix[0, 0] * vector[0, 0]) + matrix[1, 0] * vector[0, 0] * vector[1, 1]
        ) / det
        result[1, 0] = result[0, 1]
        result[1, 1] = (
            vector[1, 1] * (1 - matrix[0, 0] * vector[0, 0]) + matrix[1, 0] * vector[0, 1] * vector[0, 0]
        ) / det

        return result

    @staticmethod
    def _solver_n(matrix: np.ndarray, vector: np.ndarray, npoints: int) -> np.ndarray:
        """
        General solver for n components using numpy's linalg.solve
        """
        return np.linalg.solve(np.eye(matrix.shape[0]) - matrix, vector)


# Usage:
# solver = Solver(ncomponents)
# H_f_ij = solver.solve(E_ij - dft.ft_convolution_factor * C_f_ij, C_f_ij, ctrl['npoints'])
