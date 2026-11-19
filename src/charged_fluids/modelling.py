"""Integral-equation modelling helpers for fluids near interfaces.

This module is UI-free. Streamlit caching belongs at the app edge, not here.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import scipy.optimize as optim
from scipy import interpolate
from scipy.constants import Avogadro, Boltzmann, elementary_charge, epsilon_0
from scipy.integrate import trapezoid

from charged_fluids.memory import monitor_memory


@dataclass
class Model:
    """Discrete wall-fluid fields on a common z-grid.

    All pair/wall arrays that the solver and UI share live here (issue #21):
    ``hw``, ``tw``, ``c_short``, ``f1``, ``f2`` - not mixed between free variables
    and ad-hoc locals.
    """

    z: np.ndarray
    z_index: np.ndarray
    hw: np.ndarray
    tw: np.ndarray
    c_short: np.ndarray
    f1: np.ndarray
    f2: np.ndarray


def calc_beta(temperature: float) -> float:
    return 1.0 / (Boltzmann * temperature)


def calc_epsilon(epsilon_r: float) -> float:
    # Units match epsilon_0; distances in the model are in angstrom.
    return 4.0 * np.pi * epsilon_r * epsilon_0


def calc_l_index(i: int, j: int) -> int:
    """Pair index used by the historical code (i + j).

    Adequate for one- and two-component upper-triangle layouts used today;
    multicomponent generality will be revisited under tests.
    """
    return i + j


def calc_beta_pauling(
    valence: np.ndarray,
    n_outer_shell: np.ndarray,
    n_component: int,
    n_pair: int,
) -> np.ndarray:
    beta_pauling = np.zeros(n_pair)
    for i in range(n_component):
        for j in range(i, n_component):
            l = calc_l_index(i, j)
            beta_pauling[l] = 1.0 + valence[i] / n_outer_shell[i] + valence[j] / n_outer_shell[j]
    return beta_pauling


def calc_cap_b(
    beta_pauling: np.ndarray,
    b: float,
    alpha: float,
    sigma: list[float] | np.ndarray,
    n_component: int,
    n_pair: int,
) -> np.ndarray:
    cap_b = np.zeros(n_pair)
    for i in range(n_component):
        for j in range(i, n_component):
            l = calc_l_index(i, j)
            cap_b[l] = beta_pauling[l] * b * np.exp(alpha * (sigma[i] + sigma[j]))
    return cap_b


def calc_charge(valence: np.ndarray) -> np.ndarray:
    return valence * elementary_charge


def calc_charge_pair(
    beta: float,
    charge: np.ndarray,
    epsilon: float,
    n_component: int,
    n_pair: int,
) -> np.ndarray:
    charge_pair = np.zeros(n_pair)
    for i in range(n_component):
        for j in range(i, n_component):
            l = calc_l_index(i, j)
            charge_pair[l] = 2.0e10 * beta * charge[i] * charge[j] / epsilon
    return charge_pair


@monitor_memory
def calc_u(
    charge: np.ndarray,
    cap_b: np.ndarray,
    alpha: float,
    cap_c: np.ndarray,
    cap_d: np.ndarray,
    n_point: int,
    n_component: int,
    n_pair: int,
    epsilon: float,
    r: np.ndarray,
) -> np.ndarray:
    u = np.zeros((n_point, n_pair))
    r_inv = np.zeros_like(r)
    r_inv[1:] = 1.0 / r[1:]
    r_inv_e10 = r_inv * 1e-10
    r_inv_6 = r_inv**6
    r_inv_8 = r_inv**8
    exp_term = np.exp(-alpha * r)

    for i in range(n_component):
        for j in range(i, n_component):
            l = calc_l_index(i, j)
            charge_term = (charge[i] * charge[j]) * r_inv_e10[1:] / epsilon
            repulsive_term = cap_c[l] * r_inv_6[1:] + cap_d[l] * r_inv_8[1:]
            exp_contribution = cap_b[l] * exp_term[1:]
            u[1:, l] = charge_term + repulsive_term + exp_contribution
            u[0, l] = u[1, l]
    return u


def calc_u_lj(
    epsilon_lj: np.ndarray,
    sigma_lj: np.ndarray,
    n_point: int,
    n_component: int,
    n_pair: int,
    r: np.ndarray,
) -> np.ndarray:
    """Bulk Lennard-Jones pair potential on the radial grid ``r``."""
    u = np.zeros((n_point, n_pair))
    u[0, :] = np.inf
    for ir in range(1, n_point):
        for i in range(n_component):
            for j in range(i, n_component):
                l = calc_l_index(i, j)
                ratio = sigma_lj[l] / r[ir]
                u[ir, l] = 4.0 * epsilon_lj[l] * (ratio**12 - ratio**6)
    return u


def calc_rho(concentration: np.ndarray | list[float]) -> np.ndarray:
    """Convert mol/dm^3 to number density per angstrom^3."""
    return np.asarray(concentration, dtype=float) * Avogadro / 1.0e27


def calc_kappa(beta: float, charge: np.ndarray, rho: np.ndarray, epsilon: float) -> float:
    return float(np.sqrt(4.0 * np.pi * beta / epsilon * 1e10 * sum(np.multiply(charge**2, rho))))


def calc_phiw(z: np.ndarray, n_point: int, n_component: int) -> np.ndarray:
    phiw = np.zeros((n_point, n_component))
    capital_a = 16.274e-19  # joules
    wall_d = 2.97  # inverse angstrom
    for i in range(n_component):
        phiw[:, i] = np.exp(-wall_d * z) * capital_a * (wall_d * z + 2)
    return phiw


def interpolate_cr(
    r_in: np.ndarray,
    cr_in: np.ndarray,
    n_point: int,
    n_pair: int,
    z: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    cr = np.zeros((n_point, n_pair))
    # Normalise 1-D bulk files to shape (n_r, n_pair).
    if cr_in.ndim == 1:
        cr_in = cr_in.reshape(-1, 1)
    for l in range(n_pair):
        f = interpolate.interp1d(r_in, cr_in[:, l], bounds_error=False, fill_value="extrapolate")
        cr[:, l] = f(z)
    return cr, z


def load_and_interpolate_cr(
    cr_path: str | Path,
    n_point: int,
    n_pair: int,
    z: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    cr_df = pd.read_csv(cr_path, header=None, sep=r"\s+")
    cr_df = cr_df.set_index(0)
    r = cr_df.index.to_numpy()
    cr = cr_df.to_numpy()
    return interpolate_cr(r, cr, n_point, n_pair, z)


def calc_f1_integrand(c_short: np.ndarray, n_pair: int, z: np.ndarray, n_point: int) -> np.ndarray:
    f1_integrand = np.zeros((n_point, n_pair))
    for ij in range(n_pair):
        f1_integrand[:, ij] = z * c_short[:, ij]
    return f1_integrand


def calc_f2_integrand(c_short: np.ndarray, n_pair: int, z: np.ndarray, n_point: int) -> np.ndarray:
    f2_integrand = np.zeros((n_point, n_pair))
    for ij in range(n_pair):
        f2_integrand[:, ij] = z * z * c_short[:, ij]
    return f2_integrand


@monitor_memory
def integral_z_infty_dr_r_c_short(
    c_short: np.ndarray, n_pair: int, n_point: int, z: np.ndarray
) -> np.ndarray:
    f1 = np.zeros((n_point, n_pair))
    for ij in range(n_pair):
        zc = z * c_short[:, ij]
        for k in range(n_point):
            if k < n_point - 1:
                f1[k, ij] = trapezoid(zc[k:], z[k:])
    return f1


@monitor_memory
def integral_z_infty_dr_r2_c_short(
    c_short: np.ndarray, n_pair: int, n_point: int, z: np.ndarray
) -> np.ndarray:
    f2 = np.zeros((n_point, n_pair))
    z2 = z * z
    for ij in range(n_pair):
        z2c = z2 * c_short[:, ij]
        for k in range(n_point):
            if k < n_point - 1:
                f2[k, ij] = trapezoid(z2c[k:], z[k:])
    return f2


def calc_hw(tw: np.ndarray, n_component: int, beta_phiw: np.ndarray) -> np.ndarray:
    hw = np.zeros((len(tw), n_component))
    for i in range(n_component):
        hw[:, i] = np.exp(tw[:, i] - beta_phiw[:, i]) - 1.0
    return hw


@monitor_memory
def calc_tw(
    tw_in: np.ndarray,
    beta_phiw: np.ndarray,
    beta_psi_charge: np.ndarray,
    charge_pair: np.ndarray,
    rho: np.ndarray,
    f1: np.ndarray,
    f2: np.ndarray,
    z: np.ndarray,
    n_component: int,
    n_point: int,
    z_index: np.ndarray,
) -> np.ndarray:
    two_pi = 2.0 * np.pi
    tw = np.zeros((n_point, n_component))
    hw = np.exp(tw_in - beta_phiw) - 1.0

    def calc_integrals(hw_arr: np.ndarray, z_arr: np.ndarray, k: int, i: int):
        if k > 0:
            integral_0_z = trapezoid(hw_arr[:k, i], z_arr[:k])
        else:
            integral_0_z = 0.0
        if k < n_point:
            z_slice = z_arr[k:]
            integral_z_infty = trapezoid(z_slice * hw_arr[k:, i], z_slice)
        else:
            integral_z_infty = 0.0
        return integral_0_z, integral_z_infty

    for i in range(n_component):
        for k in range(n_point):
            integral_0_z_k, integral_z_infty_k = calc_integrals(hw, z, k, i)
            for j in range(i, n_component):
                l = calc_l_index(i, j)
                tw[k, i] = beta_psi_charge[i]
                main_terms = z[k] * f1[k, l] - f2[k, l]
                charge_terms = charge_pair[l] * (integral_z_infty_k + z[k] * integral_0_z_k)
                if k > 0:
                    trap_term1 = trapezoid(hw[:k, j] * f1[np.flip(z_index[:k]), l])
                else:
                    trap_term1 = 0.0
                if k < n_point:
                    trap_term2 = trapezoid(hw[k:, j] * f1[z_index[k:] - k, l])
                else:
                    trap_term2 = 0.0
                tw[k, i] += two_pi * rho[j] * (main_terms + charge_terms + trap_term1 + trap_term2)
    return tw


def opt_func(
    tw_in: np.ndarray,
    beta_phiw: np.ndarray,
    beta_psi_charge: np.ndarray,
    charge_pair: np.ndarray,
    rho: np.ndarray,
    f1: np.ndarray,
    f2: np.ndarray,
    z: np.ndarray,
    n_component: int,
    n_point: int,
    z_index: np.ndarray,
) -> np.ndarray:
    tw = calc_tw(
        tw_in,
        beta_phiw,
        beta_psi_charge,
        charge_pair,
        rho,
        f1,
        f2,
        z,
        n_component,
        n_point,
        z_index,
    )
    return tw_in - tw


@monitor_memory
def solve_model(
    residual: Any,
    tw_initial: np.ndarray,
    fluid: Any,
    model: Model,
    discrete: Any,
    beta_phiw: np.ndarray,
    beta_psi_charge: np.ndarray,
):
    """Solve the singlet OZ residual with SciPy's Newton-Krylov method."""
    charge_pair = fluid.charge_pair
    n_component = fluid.n_component
    rho = fluid.rho
    f1 = model.f1
    f2 = model.f2
    z = model.z
    z_index = model.z_index
    n_point = discrete.n_point
    tolerance = discrete.tolerance
    max_iteration = discrete.max_iteration
    tw_args = (
        beta_phiw,
        beta_psi_charge,
        charge_pair,
        rho,
        f1,
        f2,
        z,
        n_component,
        n_point,
        z_index,
    )
    return optim.root(
        residual,
        tw_initial,
        args=tw_args,
        method="krylov",
        jac=None,
        callback=None,
        options={"disp": True, "maxiter": max_iteration, "fatol": tolerance},
    )
