"""Sine-transform helpers for isotropic bulk OZ Fourier space."""

from __future__ import annotations

import numpy as np


def radial_grids(n_points: int, delta_r: float) -> tuple[np.ndarray, np.ndarray, float]:
    """Return r, k grids and delta_k for a sine-transform pair.

    Uses r_i = (i+1) * dr, k_i = (i+1) * dk with dk = pi / (n * dr) so that
    the discrete sine transform is unitary up to known prefactors.
    """
    r = (np.arange(n_points) + 1.0) * delta_r
    delta_k = np.pi / (n_points * delta_r)
    k = (np.arange(n_points) + 1.0) * delta_k
    return r, k, delta_k


def forward_sine(f_r: np.ndarray, r: np.ndarray, k: np.ndarray, delta_r: float) -> np.ndarray:
    """4π/k ∫_0^∞ r f(r) sin(kr) dr  (discrete)."""
    # f_r, r, k shape (n,)
    integrand = r * f_r  # length n
    # DST-II style via FFT of odd extension is overkill; use matrix sin for clarity/n small
    # For n~2048 this is fine and exact for tests
    kr = np.outer(k, r)
    f_k = 4.0 * np.pi * (delta_r / k) * (np.sin(kr) @ integrand)
    return f_k


def inverse_sine(f_k: np.ndarray, r: np.ndarray, k: np.ndarray, delta_k: float) -> np.ndarray:
    """1/(2π² r) ∫_0^∞ k f(k) sin(kr) dk  (discrete)."""
    integrand = k * f_k
    kr = np.outer(r, k)
    f_r = (delta_k / (2.0 * np.pi**2 * r)) * (np.sin(kr) @ integrand)
    return f_r
