"""Bulk Ornstein-Zernike solver (one- and two-component).

Implements Picard iteration with HNC or PY closures and optional Ng-style
long-range Coulomb separation so short-range direct correlations vanish at
large r (thesis c_sr = c + beta q_i q_j / (eps r) for bare Coulomb).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from scipy.constants import Boltzmann

from charged_fluids.bulk.transforms import forward_sine, inverse_sine, radial_grids

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
    """Upper-triangle pair index matching wall-fluid calc_l_index for n<=2."""
    return i + j


def n_pairs(n_component: int) -> int:
    return n_component * (n_component + 1) // 2


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
    closure: Closure,
    beta_u: np.ndarray,
    gamma: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return c(r), g(r) from beta*u and gamma = h - c."""
    # Soften infinite cores
    bu = np.where(np.isfinite(beta_u), beta_u, 1.0e6)
    if closure == "hnc":
        # c = exp(-bu + gamma) - 1 - gamma
        t = -bu + gamma
        t = np.clip(t, -100.0, 50.0)
        e = np.exp(t)
        c = e - 1.0 - gamma
        g = e
    elif closure == "py":
        # c = (e^{-bu}-1)(1+gamma)
        t = np.clip(-bu, -100.0, 50.0)
        f = np.exp(t) - 1.0
        c = f * (1.0 + gamma)
        g = (1.0 + f) * (1.0 + gamma)
    else:
        raise ValueError(f"Unknown closure: {closure}")
    return c, g


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


def solve_bulk_binary_coulomb_lj(
    *,
    temperature: float,
    densities: tuple[float, float],
    charges: tuple[float, float],
    epsilon: float,
    sigmas: tuple[float, float],
    epsilons_lj: tuple[float, float, float],
    n_points: int = 2048,
    delta_r: float = 0.025,
    closure: Closure = "hnc",
    mix: float = 0.25,
    max_iter: int = 800,
    tol: float = 1e-7,
) -> BulkResult:
    """Two-component charged LJ + Coulomb with Ng short-range c extraction.

    Pair order: (0,0), (0,1), (1,1) → indices 0,1,2 via i+j for n=2.
    ``epsilons_lj`` is (e00, e01, e11).
    """
    beta = beta_from_temperature(temperature)
    r, k, delta_k = radial_grids(n_points, delta_r)
    n_comp = 2
    n_pair = 3
    rho = np.asarray(densities, dtype=float)

    # Pair potentials
    beta_u = np.zeros((n_pair, n_points))
    beta_u_coul = np.zeros((n_pair, n_points))
    pair_list = [(0, 0), (0, 1), (1, 1)]
    eps_lj = {
        (0, 0): epsilons_lj[0],
        (0, 1): epsilons_lj[1],
        (1, 1): epsilons_lj[2],
    }
    for i, j in pair_list:
        l = pair_index(i, j)
        sig = 0.5 * (sigmas[i] + sigmas[j])
        u_lj = lj_potential(r, eps_lj[(i, j)], sig)
        u_c = coulomb_potential(r, charges[i], charges[j], epsilon)
        beta_u[l] = beta * (u_lj + u_c)
        beta_u_coul[l] = beta * u_c

    # Picard on gamma for each pair (approximate multi-comp OZ via diagonal density weighting)
    # Full multi-comp OZ matrix invert per k-point:
    gamma = np.zeros((n_pair, n_points))
    c = np.zeros((n_pair, n_points))
    g = np.zeros((n_pair, n_points))
    converged = False
    residual = np.inf
    n_iter = 0

    for iteration in range(1, max_iter + 1):
        n_iter = iteration
        for l in range(n_pair):
            c[l], g[l] = apply_closure(closure, beta_u[l], gamma[l])

        # Fourier c
        c_k = np.zeros((n_comp, n_comp, n_points))
        for i, j in pair_list:
            l = pair_index(i, j)
            ck = forward_sine(c[l], r, k, delta_r)
            c_k[i, j] = ck
            c_k[j, i] = ck

        gamma_k = np.zeros((n_comp, n_comp, n_points))
        for ik in range(n_points):
            # C matrix density-scaled: solve h = c + c rho h  => h = (I - c rho)^{-1} c
            C = c_k[:, :, ik]
            # rho matrix diagonal
            R = np.diag(rho)
            A = np.eye(n_comp) - C @ R
            try:
                H = np.linalg.solve(A, C)
            except np.linalg.LinAlgError:
                H = C
            G = H - C
            gamma_k[:, :, ik] = G

        gamma_new = np.zeros_like(gamma)
        for i, j in pair_list:
            l = pair_index(i, j)
            gamma_new[l] = inverse_sine(gamma_k[i, j], r, k, delta_k)

        residual = float(np.sqrt(np.mean((gamma_new - gamma) ** 2)))
        gamma = (1.0 - mix) * gamma + mix * gamma_new
        if residual < tol:
            converged = True
            break

    for l in range(n_pair):
        c[l], g[l] = apply_closure(closure, beta_u[l], gamma[l])

    # Thesis short-range: c_sr = c + beta * q_i q_j / (eps r) = c + beta_u_coul
    # (beta_u_coul already includes 1/r Coulomb)
    c_short = c + beta_u_coul

    return BulkResult(
        r=r,
        c_full=c,
        c_short=c_short,
        g=g,
        gamma=gamma,
        n_iter=n_iter,
        converged=converged,
        residual=residual,
        beta=beta,
        rho=rho,
        meta={
            "fluid": "binary_coulomb_lj",
            "closure": closure,
            "temperature": temperature,
            "densities": list(densities),
            "charges": list(charges),
            "epsilon": epsilon,
            "sigmas": list(sigmas),
            "epsilons_lj": list(epsilons_lj),
            "n_points": n_points,
            "delta_r": delta_r,
            "c_short_definition": "c + beta * u_coulomb",
        },
    )
