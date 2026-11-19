"""Sine-transform helpers for isotropic bulk OZ Fourier space.

Default path uses an O(n log n) FFT of the odd extension, which matches the
historical O(n²) ``sin(kr) @ v`` matrix product exactly on the standard grid
``r_i = (i+1) dr``, ``k_i = (i+1) dk``, ``dk = π / (n dr)``.
"""

from __future__ import annotations

from typing import Literal

import numpy as np

TransformBackend = Literal["fft", "direct"]

# Below this size the dense product is cheap and keeps tests simple to debug.
_FFT_MIN_POINTS = 32


def radial_grids(n_points: int, delta_r: float) -> tuple[np.ndarray, np.ndarray, float]:
    """Return r, k grids and delta_k for a sine-transform pair.

    Uses r_i = (i+1) * dr, k_i = (i+1) * dk with dk = pi / (n * dr) so that
    the discrete sine transform is unitary up to known prefactors.
    """
    r = (np.arange(n_points) + 1.0) * delta_r
    delta_k = np.pi / (n_points * delta_r)
    k = (np.arange(n_points) + 1.0) * delta_k
    return r, k, delta_k


def _sin_matmul_direct(left: np.ndarray, right: np.ndarray, vec: np.ndarray) -> np.ndarray:
    """Dense ``sin(outer(left, right)) @ vec``."""
    return np.sin(np.outer(left, right)) @ vec


def _sin_matmul_fft(vec: np.ndarray) -> np.ndarray:
    """Compute ``sum_i sin(π (j+1)(i+1)/n) v_i`` via FFT of the odd extension.

    Exactly matches the dense product on the standard ``(i+1)`` grid (factor 1/2
    from the full complex FFT of a pure-odd real sequence).
    """
    n = vec.shape[0]
    y = np.zeros(2 * n, dtype=float)
    y[1 : n + 1] = vec
    # y[n] contributes 0 to every mode (sin(j π) = 0); mirror the first n-1 samples.
    if n > 1:
        y[n + 1 :] = -np.flip(vec[: n - 1])
    spectrum = np.fft.fft(y)
    return -0.5 * spectrum.imag[1 : n + 1]


def sin_matmul(
    left: np.ndarray,
    right: np.ndarray,
    vec: np.ndarray,
    *,
    backend: TransformBackend = "fft",
) -> np.ndarray:
    """Matrix-free or dense sine product used by forward/inverse transforms."""
    n = vec.shape[0]
    use_fft = backend == "fft" and n >= _FFT_MIN_POINTS
    # FFT path assumes the standard commensurate grid from radial_grids.
    if use_fft:
        return _sin_matmul_fft(vec)
    return _sin_matmul_direct(left, right, vec)


def forward_sine(
    f_r: np.ndarray,
    r: np.ndarray,
    k: np.ndarray,
    delta_r: float,
    *,
    backend: TransformBackend = "fft",
) -> np.ndarray:
    """4π/k ∫_0^∞ r f(r) sin(kr) dr  (discrete)."""
    integrand = r * f_r
    sine = sin_matmul(k, r, integrand, backend=backend)
    return 4.0 * np.pi * (delta_r / k) * sine


def inverse_sine(
    f_k: np.ndarray,
    r: np.ndarray,
    k: np.ndarray,
    delta_k: float,
    *,
    backend: TransformBackend = "fft",
) -> np.ndarray:
    """1/(2π² r) ∫_0^∞ k f(k) sin(kr) dk  (discrete)."""
    integrand = k * f_k
    sine = sin_matmul(r, k, integrand, backend=backend)
    return (delta_k / (2.0 * np.pi**2 * r)) * sine
