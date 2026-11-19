"""Bulk Ornstein-Zernike solver (one- and two-component).

Implements Picard iteration with HNC or PY closures and optional Ng-style
long-range Coulomb separation so short-range direct correlations vanish at
large r (thesis c_sr = c + beta q_i q_j / (eps r) for bare Coulomb).

Closures live in :mod:`charged_fluids.bulk.closures`; sine transforms in
:mod:`charged_fluids.bulk.transforms`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from scipy.constants import Boltzmann

from charged_fluids.bulk.closures import apply_closure as _apply_closure
from charged_fluids.bulk.transforms import forward_sine, inverse_sine, radial_grids
from charged_fluids.modelling import calc_l_index, n_pairs

Closure = Literal["hnc", "py"]


@dataclass
class BulkResult:
    r: np.ndarray
    c_full: np.ndarray  # shape (n_pair, n_r) upper-triangle order (0,0), (0,1), (1,1), ...
    c_short: np.ndarray
    g: np.ndarray
    gamma: np.ndarray
    n_iter: int
    converged: bool
    residual: float
    beta: float
    rho: np.ndarray
    meta: dict


def pair_index(i: int, j: int) -> int:
    """Upper-triangle pair index matching wall-fluid :func:`calc_l_index`."""
    return calc_l_index(i, j)


def beta_from_temperature(temperature: float) -> float:
    return 1.0 / (Boltzmann * temperature)


def lj_potential(r: np.ndarray, epsilon: float, sigma: float) -> np.ndarray:
    ratio = np.where(r > 0, sigma / r, np.inf)
    u = 4.0 * epsilon * (ratio**12 - ratio**6)
    u = np.where(np.isfinite(u), u, np.inf)
    return u


def coulomb_potential(r: np.ndarray, qi: float, qj: float, epsilon: float) -> np.ndarray:
    """Coulomb pair energy matching wall-fluid calc_u units (J)."""
    out = np.zeros_like(r, dtype=float)
    mask = r > 0
    out[mask] = qi * qj * 1.0e-10 / (epsilon * r[mask])
    return out


def apply_closure(
    closure: Closure | str,
    beta_u: np.ndarray,
    gamma: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return c(r), g(r) from beta*u and gamma = h - c (registry-backed)."""
    return _apply_closure(closure, beta_u, gamma)


def oz_one_component_fourier(c_k: np.ndarray, rho: float) -> tuple[np.ndarray, np.ndarray]:
    """h_hat = c_hat / (1 - rho c_hat); gamma_hat = h_hat - c_hat."""
    denom = 1.0 - rho * c_k
    # Avoid division blow-ups
    denom = np.where(np.abs(denom) < 1e-14, np.sign(denom) * 1e-14 + 1e-14, denom)
    h_k = c_k / denom
    gamma_k = h_k - c_k
    return h_k, gamma_k


def solve_bulk_lj(
    *,
    temperature: float,
    density: float,
    epsilon: float,
    sigma: float,
    n_points: int = 2048,
    delta_r: float = 0.02,
    closure: Closure = "hnc",
    mix: float = 0.5,
    max_iter: int = 500,
    tol: float = 1e-8,
) -> BulkResult:
    """Solve bulk OZ for a one-component Lennard-Jones fluid.

    Parameters
    ----------
    density:
        Number density in particles per angstrom^3.
    epsilon, sigma:
        LJ parameters in energy (J) and angstrom. For reduced units set
        temperature so that beta*epsilon is the desired reduced epsilon.
    """
    beta = beta_from_temperature(temperature)
    r, k, delta_k = radial_grids(n_points, delta_r)
    u = lj_potential(r, epsilon, sigma)
    beta_u = beta * u

    gamma = np.zeros_like(r)
    converged = False
    residual = np.inf
    n_iter = 0
    c = np.zeros_like(r)
    g = np.zeros_like(r)

    for iteration in range(1, max_iter + 1):
        n_iter = iteration
        c, g = apply_closure(closure, beta_u, gamma)
        c_k = forward_sine(c, r, k, delta_r)
        _, gamma_k = oz_one_component_fourier(c_k, density)
        gamma_new = inverse_sine(gamma_k, r, k, delta_k)
        residual = float(np.sqrt(np.mean((gamma_new - gamma) ** 2)))
        gamma = (1.0 - mix) * gamma + mix * gamma_new
        if residual < tol:
            converged = True
            break

    c, g = apply_closure(closure, beta_u, gamma)
    # Neutral LJ: short-range is full c
    c_pair = c.reshape(1, -1)
    g_pair = g.reshape(1, -1)
    gamma_pair = gamma.reshape(1, -1)

    return BulkResult(
        r=r,
        c_full=c_pair,
        c_short=c_pair.copy(),
        g=g_pair,
        gamma=gamma_pair,
        n_iter=n_iter,
        converged=converged,
        residual=residual,
        beta=beta,
        rho=np.array([density]),
        meta={
            "fluid": "lj1",
            "closure": closure,
            "epsilon": epsilon,
            "sigma": sigma,
            "temperature": temperature,
            "density": density,
            "n_points": n_points,
            "delta_r": delta_r,
        },
    )


def dimensionless_coulomb_from_charge_pair(r: np.ndarray, charge_pair: np.ndarray) -> np.ndarray:
    """Dimensionless beta*u_coulomb matching wall-fluid ``charge_pair`` scale.

    ``charge_pair[l] = 2e10 * beta * q_i * q_j / epsilon`` from
    :func:`charged_fluids.modelling.calc_charge_pair`, so
    ``beta u_c(r) = charge_pair[l] / (2 r)``.
    """
    charge_pair = np.asarray(charge_pair, dtype=float).ravel()
    out = np.zeros((charge_pair.size, r.size))
    mask = r > 0
    # broadcast: (n_pair, 1) / (1, n_r)
    out[:, mask] = charge_pair[:, None] / (2.0 * r[None, mask])
    return out


def solve_bulk_binary_coulomb_lj(
    *,
    temperature: float,
    densities: tuple[float, float],
    charge_pair: np.ndarray,
    sigmas: tuple[float, float],
    epsilons_lj: tuple[float, float, float],
    n_points: int = 2048,
    delta_r: float = 0.025,
    closure: Closure = "hnc",
    mix: float = 0.15,
    max_iter: int = 1500,
    tol: float = 1e-6,
    ng_alpha: float = 1.2,
) -> BulkResult:
    """Two-component LJ + Coulomb bulk OZ with Ng long-range renormalisation.

    Pair order: (0,0), (0,1), (1,1) via :func:`pair_index` for n=2.
    Coulomb uses wall-fluid ``charge_pair`` units. Ng splits
    ``beta u_c = beta u_long + beta u_c_short`` with
    ``u_long ~ erf(alpha r)/r`` so the iterated potential is short-ranged;
    final ``c_short = c_full + beta u_c`` (thesis) vanishes at large r when
    ``c_full -> -beta u_c``.
    """
    from scipy.special import erf

    beta = beta_from_temperature(temperature)
    r, k, delta_k = radial_grids(n_points, delta_r)
    n_comp = 2
    n_pair = n_pairs(n_comp)
    rho = np.asarray(densities, dtype=float)
    charge_pair = np.asarray(charge_pair, dtype=float).ravel()
    if charge_pair.size != n_pair:
        raise ValueError(f"charge_pair must have length {n_pair}, got {charge_pair.size}")

    pair_list = [(0, 0), (0, 1), (1, 1)]
    eps_lj = {
        (0, 0): epsilons_lj[0],
        (0, 1): epsilons_lj[1],
        (1, 1): epsilons_lj[2],
    }

    beta_u_c = dimensionless_coulomb_from_charge_pair(r, charge_pair)
    # Ng long-range piece (screened Coulomb)
    beta_u_long = np.zeros_like(beta_u_c)
    beta_u_short_c = np.zeros_like(beta_u_c)
    erf_term = np.where(r > 0, erf(ng_alpha * r) / (2.0 * r), ng_alpha / np.sqrt(np.pi))
    beta_u_long = charge_pair[:, None] * erf_term[None, :]
    beta_u_short_c = beta_u_c - beta_u_long

    beta_u_s = np.zeros((n_pair, n_points))
    for i, j in pair_list:
        l = pair_index(i, j)
        sig = 0.5 * (sigmas[i] + sigmas[j])
        u_lj = lj_potential(r, eps_lj[(i, j)], sig)
        beta_u_s[l] = beta * u_lj + beta_u_short_c[l]

    # Initial gamma: cancel long-range so first closure is better conditioned
    gamma = -beta_u_long.copy()
    c = np.zeros((n_pair, n_points))
    g = np.zeros((n_pair, n_points))
    converged = False
    residual = np.inf
    n_iter = 0

    for iteration in range(1, max_iter + 1):
        n_iter = iteration
        for l in range(n_pair):
            # HNC/PY on short-range potential + short-range gamma
            c[l], g[l] = apply_closure(closure, beta_u_s[l], gamma[l])

        c_k = np.zeros((n_comp, n_comp, n_points))
        for i, j in pair_list:
            l = pair_index(i, j)
            ck = forward_sine(c[l], r, k, delta_r)
            c_k[i, j] = ck
            c_k[j, i] = ck

        gamma_k = np.zeros((n_comp, n_comp, n_points))
        for ik in range(n_points):
            C = c_k[:, :, ik]
            R = np.diag(rho)
            A = np.eye(n_comp) - C @ R
            try:
                H = np.linalg.solve(A, C)
            except np.linalg.LinAlgError:
                H = C
            gamma_k[:, :, ik] = H - C

        gamma_new = np.zeros_like(gamma)
        for i, j in pair_list:
            l = pair_index(i, j)
            gamma_new[l] = inverse_sine(gamma_k[i, j], r, k, delta_k)

        residual = float(np.sqrt(np.mean((gamma_new - gamma) ** 2)))
        # Adaptive mix when residual grows
        step = mix if residual < 1.0 else mix * 0.25
        gamma = (1.0 - step) * gamma + step * gamma_new
        if residual < tol:
            converged = True
            break

    for l in range(n_pair):
        c[l], g[l] = apply_closure(closure, beta_u_s[l], gamma[l])

    # Full direct correlation approx: short-range c minus long-range beta u
    # (Ng: c = c_s - beta u_long). Thesis wall input wants c_short = c + beta u_c
    # with beta u_c full Coulomb, equivalent to c_s + beta u_short_c when
    # c = c_s - beta u_long and beta u_c = beta u_long + beta u_short_c.
    c_full = c - beta_u_long
    c_short = c + beta_u_short_c  # = c_full + beta_u_c

    return BulkResult(
        r=r,
        c_full=c_full,
        c_short=c_short,
        g=g,
        gamma=gamma,
        n_iter=n_iter,
        converged=converged,
        residual=residual,
        beta=beta,
        rho=rho,
        meta={
            "fluid": "binary_coulomb_lj_ng",
            "closure": closure,
            "temperature": temperature,
            "densities": list(map(float, densities)),
            "charge_pair": charge_pair.tolist(),
            "sigmas": list(sigmas),
            "epsilons_lj": list(epsilons_lj),
            "n_points": n_points,
            "delta_r": delta_r,
            "ng_alpha": ng_alpha,
            "c_short_definition": "c_s + beta_u_c_short (= c_full + beta_u_coulomb)",
        },
    )
