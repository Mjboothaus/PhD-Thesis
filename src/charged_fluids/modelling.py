"""Integral-equation modelling helpers for fluids near interfaces.

This module is UI-free. Streamlit caching belongs at the app edge, not here.

Unit conventions (Å grids, joule potentials, charge_pair scale) are documented
in :mod:`charged_fluids.units`.
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
from scipy.integrate import cumulative_trapezoid, trapezoid

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
    """Upper-triangle pair index including the diagonal (symmetric in i, j).

    Uses ``l = j*(j+1)//2 + i`` with ``i <= j``. For one- and two-component
    fluids this matches the historical ``i + j`` layout of thesis data files.
    For ``n > 2`` it avoids collisions that ``i + j`` would produce.
    """
    if i > j:
        i, j = j, i
    return j * (j + 1) // 2 + i


def n_pairs(n_component: int) -> int:
    """Number of unique unordered pairs including like-like terms."""
    return n_component * (n_component + 1) // 2


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
    if n_point < 2:
        return u
    r_tail = r[1:, None]  # (n_point-1, 1)
    sigma = np.asarray(sigma_lj, dtype=float).reshape(1, -1)
    eps = np.asarray(epsilon_lj, dtype=float).reshape(1, -1)
    ratio = sigma / r_tail
    u[1:, :] = 4.0 * eps * (ratio**12 - ratio**6)
    # Only fill pair columns that belong to the component layout.
    if n_pair != n_pairs(n_component):
        # Defensive: fall back to historical loop if sizes disagree.
        u = np.zeros((n_point, n_pair))
        u[0, :] = np.inf
        for ir in range(1, n_point):
            for i in range(n_component):
                for j in range(i, n_component):
                    l = calc_l_index(i, j)
                    ratio_ij = sigma_lj[l] / r[ir]
                    u[ir, l] = 4.0 * epsilon_lj[l] * (ratio_ij**12 - ratio_ij**6)
    return u


def calc_rho(concentration: np.ndarray | list[float]) -> np.ndarray:
    """Convert mol/dm^3 to number density per angstrom^3."""
    return np.asarray(concentration, dtype=float) * Avogadro / 1.0e27


def calc_kappa(beta: float, charge: np.ndarray, rho: np.ndarray, epsilon: float) -> float:
    return float(np.sqrt(4.0 * np.pi * beta / epsilon * 1e10 * sum(np.multiply(charge**2, rho))))


def calc_phiw(z: np.ndarray, n_point: int, n_component: int) -> np.ndarray:
    capital_a = 16.274e-19  # joules
    wall_d = 2.97  # inverse angstrom
    # Same wall potential on every component (historical thesis form).
    profile = np.exp(-wall_d * z) * capital_a * (wall_d * z + 2)
    return np.column_stack([profile] * n_component) if n_component else np.zeros((n_point, 0))


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
    del n_pair, n_point  # shapes come from arrays
    return z[:, None] * c_short


def calc_f2_integrand(c_short: np.ndarray, n_pair: int, z: np.ndarray, n_point: int) -> np.ndarray:
    del n_pair, n_point
    return (z * z)[:, None] * c_short


def _integral_z_to_end(values: np.ndarray, z: np.ndarray) -> np.ndarray:
    """For each k, trapezoid(values[k:], z[k:]) via one cumulative pass.

    ``values`` and result have shape (n_point,) or (n_point, n_cols).
    """
    values = np.asarray(values, dtype=float)
    z = np.asarray(z, dtype=float)
    n_point = z.size
    if n_point < 2:
        return np.zeros_like(values, dtype=float)

    if values.ndim == 1:
        cum = np.zeros(n_point, dtype=float)
        cum[1:] = cumulative_trapezoid(values, z)
        out = np.zeros(n_point, dtype=float)
        out[:-1] = cum[-1] - cum[:-1]
        return out

    # Multi-column: one cumulative_trapezoid call via ravel is awkward; use axis loop
    # but avoid per-k Python trapezoid. cumulative_trapezoid supports axis=.
    cum = np.zeros_like(values, dtype=float)
    cum[1:, :] = cumulative_trapezoid(values, z, axis=0)
    out = np.zeros_like(values, dtype=float)
    out[:-1, :] = cum[-1:, :] - cum[:-1, :]
    return out


def _trapz_unit_dx(y: np.ndarray) -> float:
    """Trapezoid with unit sample spacing (historical wall-kernel convention).

    The original thesis discretisation called ``trapezoid(y)`` *without* the ``z``
    abscissa for the hwxf1 convolution pieces, i.e. dx = 1 in index space.
    Charge-path integrals still use physical ``z``.
    """
    if y.size < 2:
        return 0.0
    return float(trapezoid(y))


def _wall_pair_kernels(
    hw_j: np.ndarray,
    f1_l: np.ndarray,
    z: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Precompute trap_term1[k] and trap_term2[k] for one (j, pair) channel.

    Matches the historical index-space trapezoid (unit dx) used by the thesis code::

        trap_term1[k] ~ trapz(hw[:k] * flip(f1[:k]))
        trap_term2[k] ~ trapz(hw[k:] * f1[:n-k])
    """
    del z  # physical z is not used for these kernels (unit dx)
    n = hw_j.size
    t1 = np.zeros(n, dtype=float)
    t2 = np.zeros(n, dtype=float)
    for k in range(n):
        if k > 0:
            t1[k] = _trapz_unit_dx(hw_j[:k] * np.flip(f1_l[:k]))
        # k == n-1 -> single sample -> 0; still call for parity with historical branch
        t2[k] = _trapz_unit_dx(hw_j[k:] * f1_l[: n - k])
    return t1, t2


@monitor_memory
def integral_z_infty_dr_r_c_short(
    c_short: np.ndarray, n_pair: int, n_point: int, z: np.ndarray
) -> np.ndarray:
    del n_pair, n_point
    return _integral_z_to_end(z[:, None] * c_short, z)


@monitor_memory
def integral_z_infty_dr_r2_c_short(
    c_short: np.ndarray, n_pair: int, n_point: int, z: np.ndarray
) -> np.ndarray:
    del n_pair, n_point
    return _integral_z_to_end((z * z)[:, None] * c_short, z)


def calc_hw(tw: np.ndarray, n_component: int, beta_phiw: np.ndarray) -> np.ndarray:
    del n_component
    return np.exp(tw - beta_phiw) - 1.0


def _prefix_trapezoid(y: np.ndarray, z: np.ndarray) -> np.ndarray:
    """prefix[k] = trapezoid(y[:k], z[:k]) with prefix[0] = 0."""
    n = z.size
    out = np.zeros(n, dtype=float)
    if n < 2:
        return out
    cum = np.zeros(n, dtype=float)
    cum[1:] = cumulative_trapezoid(y, z)
    # trapezoid(y[:k], z[:k]) == cum[k-1] for k>=1 when endpoints match...
    # cumulative_trapezoid: cum[k] = int from z[0] to z[k] for k>=1, cum[0]=0.
    # trapezoid(y[:k], z[:k]) integrates points 0..k-1, i.e. up to z[k-1] = cum[k-1].
    out[1:] = cum[:-1]
    # For k=1: trapezoid(y[:1], z[:1]) = 0 (single point) - cum[0]=0, out[1]=0 OK
    # For k=2: trapezoid(y[:2], z[:2]) = cum[1], out[2]=cum[1] OK
    return out


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
    """Singlet wall indirect correlation from bulk short-range kernels.

    For species ``i``::

        tw_i(z) = beta*psi_i
                  + sum_j 2pi rho_j [ z f1_ij - f2_ij
                    + charge_pair_ij * (int _z^∞ z' hw_i dz' + z int _0^z hw_i dz')
                    + int _0^z hw_j(z') f1_ij(z-z') dz'
                    + int _z^∞ hw_j(z') f1_ij(z'-z) dz' ]

    Pair index ``l = calc_l_index(i, j)`` is symmetric. The sum runs over **all**
    partner species ``j`` (correct multi-component form).
    """
    del z_index  # kept in signature for call-site compatibility
    two_pi = 2.0 * np.pi
    tw = np.zeros((n_point, n_component))
    hw = np.exp(tw_in - beta_phiw) - 1.0

    # Charge-path integrals per species: int _0^z hw and int _z^∞ z hw
    integral_0_z = np.zeros((n_point, n_component))
    integral_z_infty = np.zeros((n_point, n_component))
    for i in range(n_component):
        integral_0_z[:, i] = _prefix_trapezoid(hw[:, i], z)
        integral_z_infty[:, i] = _integral_z_to_end(z * hw[:, i], z)

    # Precompute wall pair kernels for every (j, pair channel)
    # Cache by (j, l) since l depends on (i,j)
    kernel_cache: dict[tuple[int, int], tuple[np.ndarray, np.ndarray]] = {}

    for i in range(n_component):
        tw[:, i] = beta_psi_charge[i]
        for j in range(n_component):
            l = calc_l_index(i, j)
            key = (j, l)
            if key not in kernel_cache:
                kernel_cache[key] = _wall_pair_kernels(hw[:, j], f1[:, l], z)
            t1, t2 = kernel_cache[key]
            main = z * f1[:, l] - f2[:, l]
            charge = charge_pair[l] * (integral_z_infty[:, i] + z * integral_0_z[:, i])
            tw[:, i] += two_pi * rho[j] * (main + charge + t1 + t2)
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
