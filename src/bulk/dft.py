# File: pyoz_dft.py

import numpy as np
from scipy.special import jn
from typing import Tuple


class DFT:
    def __init__(self, npoints: int, deltar: float, deltak: float, r: np.ndarray, k: np.ndarray):
        self.npoints = npoints
        self.deltar = deltar
        self.deltak = deltak
        self.r = r
        self.k = k
        self.ft_convolution_factor = 2 * np.pi * deltar * deltak
        self._initialize_transform_matrices()

    def _initialize_transform_matrices(self):
        self.fbt_matrix = np.zeros((self.npoints, self.npoints))
        self.ifbt_matrix = np.zeros((self.npoints, self.npoints))

        for i in range(self.npoints):
            for j in range(self.npoints):
                self.fbt_matrix[i, j] = self.r[j] * jn(0, self.k[i] * self.r[j])
                self.ifbt_matrix[i, j] = self.k[j] * jn(0, self.k[j] * self.r[i])

        self.fbt_matrix *= self.deltar
        self.ifbt_matrix *= self.deltak

    def dfbt(self, func: np.ndarray, norm: float, corr: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Discrete Fourier-Bessel Transform
        \tilde{f}(k_i) = \Delta r \sum_{j=1}^{N} r_j f(r_j) J_0(k_i r_j)

        \text{Where:}
        \begin{itemize}
            \item $f(r)$ is the function in real space
            \item $\tilde{f}(k)$ is the transformed function in reciprocal space
            \item $J_0$ is the Bessel function of the first kind of order zero
            \item $\Delta r$ and $\Delta k$ are the step sizes in real and reciprocal space respectively
            \item $N$ is the number of points
            \item $r_i$ and $k_i$ are the discrete points in real and reciprocal space
        \end{itemize}
        """
        cs_f = np.dot(self.fbt_matrix, func * self.r)
        c_f = cs_f + corr
        return cs_f * norm, c_f * norm

    def idfbt(self, func: np.ndarray, norm: float, corr: np.ndarray) -> np.ndarray:
        """
        Inverse Discrete Fourier-Bessel Transform
        f(r_i) = \Delta k \sum_{j=1}^{N} k_j \tilde{f}(k_j) J_0(k_j r_i)
        See details above.
        """
        return np.dot(self.ifbt_matrix, func * self.k) / norm + corr

    def print_status(self):
        print("DFT initialised:")
        print(f"  Number of points: {self.npoints}")
        print(f"  Delta r: {self.deltar:.6f}")
        print(f"  Delta k: {self.deltak:.6f}")
        print(f"  FT convolution factor: {self.ft_convolution_factor:.6f}")


# Usage example:
# dft = DFT(npoints, deltar, deltak, r, k)
# cs_f, c_f = dft.dfbt(func, norm, corr)
# g_r = dft.idfbt(func, norm, corr)
