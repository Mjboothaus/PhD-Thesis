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
    anneal_steps: int = 1,
    use_anderson: bool = True,
    t_hot: float | None = None,
) -> BulkResult:
    """Two-component LJ + Coulomb bulk OZ with Ng long-range renormalisation.

    Pair order: (0,0), (0,1), (1,1) via :func:`pair_index` for n=2.
    Coulomb uses wall-fluid ``charge_pair`` units. Ng splits
    ``beta u_c = beta u_long + beta u_c_short`` with
    ``u_long ~ erf(alpha r)/r`` so the iterated potential is short-ranged;
    final ``c_short = c_full + beta u_c`` (thesis) vanishes at large r when
    ``c_full -> -beta u_c``.

    For dense melts (e.g. KCl), set ``anneal_steps > 1`` to cool from a hot
    start with continuous ``gamma`` (rescaled by β). Optional Anderson mixing
    stabilises the late Picard iterations.
    """
    from scipy.special import erf

    t_final = float(temperature)
    beta_final = beta_from_temperature(t_final)
    r, k, delta_k = radial_grids(n_points, delta_r)
    n_comp = 2
    n_pair = n_pairs(n_comp)
    rho = np.asarray(densities, dtype=float)
    charge_pair_final = np.asarray(charge_pair, dtype=float).ravel()
    if charge_pair_final.size != n_pair:
        raise ValueError(f"charge_pair must have length {n_pair}, got {charge_pair_final.size}")

    pair_list = [(0, 0), (0, 1), (1, 1)]
    eps_lj = {
        (0, 0): epsilons_lj[0],
        (0, 1): epsilons_lj[1],
        (1, 1): epsilons_lj[2],
    }

    if anneal_steps <= 1:
        temps = [t_final]
    else:
        hot = float(t_hot) if t_hot is not None else max(t_final * 2.5, t_final + 800.0)
        temps = list(np.geomspace(hot, t_final, int(anneal_steps)))

    gamma = np.zeros((n_pair, n_points))
    c = np.zeros((n_pair, n_points))
    g = np.zeros((n_pair, n_points))
    converged = False
    residual = np.inf
    n_iter_total = 0
    prev_t: float | None = None

    for i_t, t_step in enumerate(temps):
        t_step = float(t_step)
        beta = beta_from_temperature(t_step)
        scale = beta / beta_final
        cp = charge_pair_final * scale
        erf_term = np.where(r > 0, erf(ng_alpha * r) / (2.0 * r), ng_alpha / np.sqrt(np.pi))
        beta_u_c = dimensionless_coulomb_from_charge_pair(r, cp)
        beta_u_long = cp[:, None] * erf_term[None, :]
        beta_u_short_c = beta_u_c - beta_u_long

        beta_u_s = np.zeros((n_pair, n_points))
        for i, j in pair_list:
            l = pair_index(i, j)
            sig = 0.5 * (sigmas[i] + sigmas[j])
            u_lj = lj_potential(r, eps_lj[(i, j)], sig)
            beta_u_s[l] = beta * u_lj + beta_u_short_c[l]

        if prev_t is None:
            # γ=0 is more stable for the first hot stage than −βu_long at melt density
            gamma = np.zeros((n_pair, n_points))
        else:
            gamma = gamma * (beta / beta_from_temperature(prev_t))
        prev_t = t_step

        is_final = i_t == len(temps) - 1
        if is_final:
            mix_stage, stage_iters, stage_tol = mix, max_iter, tol
        elif t_step >= max(t_final * 1.4, 1500.0):
            # Hot stages: larger Picard steps (melt KCl needs ~0.2 here)
            mix_stage = min(0.35, max(mix * 5.0, 0.2))
            stage_iters = min(max_iter, 800)
            stage_tol = max(tol * 100.0, 1e-3)
        else:
            mix_stage = min(0.15, max(mix * 2.0, 0.08))
            stage_iters = min(max_iter, 2000)
            stage_tol = max(tol * 20.0, 2e-4)

        hist_g: list[np.ndarray] = []
        hist_f: list[np.ndarray] = []
        best: tuple[float, np.ndarray, np.ndarray, np.ndarray] | None = None

        for _iteration in range(1, stage_iters + 1):
            n_iter_total += 1
            for l in range(n_pair):
                c[l], g[l] = apply_closure(closure, beta_u_s[l], gamma[l])
                g[l] = np.maximum(g[l], 0.0)

            c_k = np.zeros((n_comp, n_comp, n_points))
            for i, j in pair_list:
                l = pair_index(i, j)
                ck = forward_sine(c[l], r, k, delta_r)
                c_k[i, j] = ck
                c_k[j, i] = ck

            gamma_k = _binary_oz_gamma_k(c_k, rho)
            gamma_new = np.zeros_like(gamma)
            for i, j in pair_list:
                l = pair_index(i, j)
                gamma_new[l] = inverse_sine(gamma_k[i, j], r, k, delta_k)

            if not np.isfinite(gamma_new).all():
                if best is not None:
                    residual, gamma, c, g = best[0], best[1].copy(), best[2].copy(), best[3].copy()
                converged = False
                break

            residual = float(np.sqrt(np.mean((gamma_new - gamma) ** 2)))
            if best is None or residual < best[0]:
                best = (residual, gamma.copy(), c.copy(), g.copy())

            if use_anderson and stage_tol < residual < 1.0:
                gamma = _anderson_mix(
                    gamma, gamma_new, hist_g, hist_f, m_hist=6, damp=mix_stage
                )
            else:
                if residual > 10.0:
                    step = mix_stage * 0.1
                elif residual > 1.0:
                    step = mix_stage * 0.4
                else:
                    step = mix_stage
                gamma = (1.0 - step) * gamma + step * gamma_new
                hist_g.clear()
                hist_f.clear()

            if residual < stage_tol:
                converged = is_final
                break
        else:
            converged = False
            if best is not None:
                residual, gamma, c, g = best[0], best[1].copy(), best[2].copy(), best[3].copy()

    # Final closure at target T
    beta = beta_final
    erf_term = np.where(r > 0, erf(ng_alpha * r) / (2.0 * r), ng_alpha / np.sqrt(np.pi))
    beta_u_c = dimensionless_coulomb_from_charge_pair(r, charge_pair_final)
    beta_u_long = charge_pair_final[:, None] * erf_term[None, :]
    beta_u_short_c = beta_u_c - beta_u_long
    beta_u_s = np.zeros((n_pair, n_points))
    for i, j in pair_list:
        l = pair_index(i, j)
        sig = 0.5 * (sigmas[i] + sigmas[j])
        u_lj = lj_potential(r, eps_lj[(i, j)], sig)
        beta_u_s[l] = beta * u_lj + beta_u_short_c[l]
    for l in range(n_pair):
        c[l], g[l] = apply_closure(closure, beta_u_s[l], gamma[l])
        g[l] = np.maximum(g[l], 0.0)

    c_full = c - beta_u_long
    c_short = c + beta_u_short_c  # = c_full + beta_u_c
    tail = float(np.max(np.abs(c_short[:, -max(5, n_points // 20) :])))

    return BulkResult(
        r=r,
        c_full=c_full,
        c_short=c_short,
        g=g,
        gamma=gamma,
        n_iter=n_iter_total,
        converged=converged,
        residual=residual,
        beta=beta,
        rho=rho,
        meta={
            "fluid": "binary_coulomb_lj_ng",
            "closure": closure,
            "temperature": t_final,
            "densities": list(map(float, densities)),
            "charge_pair": charge_pair_final.tolist(),
            "sigmas": list(sigmas),
            "epsilons_lj": list(epsilons_lj),
            "n_points": n_points,
            "delta_r": delta_r,
            "ng_alpha": ng_alpha,
            "anneal_steps": len(temps),
            "use_anderson": use_anderson,
            "tail_abs_max_c_short": tail,
            "g_peak": [float(np.max(g[l])) for l in range(n_pair)],
            "c_short_definition": "c_s + beta_u_c_short (= c_full + beta_u_coulomb)",
        },
    )


def _binary_oz_gamma_k(
    c_k: np.ndarray,
    rho: np.ndarray,
) -> np.ndarray:
    """gamma_k = H - C for binary OZ in k-space; c_k shape (2,2,n_k)."""
    r0, r1 = float(rho[0]), float(rho[1])
    c00, c01 = c_k[0, 0], c_k[0, 1]
    c10, c11 = c_k[1, 0], c_k[1, 1]
    a00 = 1.0 - c00 * r0
    a01 = -c01 * r1
    a10 = -c10 * r0
    a11 = 1.0 - c11 * r1
    det = a00 * a11 - a01 * a10
    det = np.where(np.abs(det) < 1e-14, np.sign(det) * 1e-14 + 1e-14, det)
    inv00 = a11 / det
    inv01 = -a01 / det
    inv10 = -a10 / det
    inv11 = a00 / det
    h00 = inv00 * c00 + inv01 * c10
    h01 = inv00 * c01 + inv01 * c11
    h10 = inv10 * c00 + inv11 * c10
    h11 = inv10 * c01 + inv11 * c11
    gamma_k = np.empty_like(c_k)
    gamma_k[0, 0] = h00 - c00
    gamma_k[0, 1] = h01 - c01
    gamma_k[1, 0] = h10 - c10
    gamma_k[1, 1] = h11 - c11
    return gamma_k


def _anderson_mix(
    gamma: np.ndarray,
    gamma_new: np.ndarray,
    hist_g: list[np.ndarray],
    hist_f: list[np.ndarray],
    *,
    m_hist: int = 5,
    damp: float = 0.15,
) -> np.ndarray:
    """Simple Anderson (DIIS-like) mixing on flattened pair fields."""
    f = gamma_new - gamma
    hist_g.append(gamma.copy())
    hist_f.append(f.copy())
    if len(hist_g) > m_hist:
        hist_g.pop(0)
        hist_f.pop(0)
    m = len(hist_f)
    if m < 2:
        return (1.0 - damp) * gamma + damp * gamma_new
    # Solve min || sum c_i F_i || with sum c_i = 1
    dim = f.size
    F = np.column_stack([fi.ravel() for fi in hist_f])  # (dim, m)
    # Regularised least squares for coefficients with constraint via Lagrange
    A = F.T @ F
    A = A + 1e-8 * np.eye(m)
    ones = np.ones(m)
    try:
        # Block system: [A 1; 1^T 0] [c; λ] = [0; 1]
        top = np.hstack([A, ones[:, None]])
        bot = np.hstack([ones, np.zeros(1)])
        mat = np.vstack([top, bot])
        rhs = np.zeros(m + 1)
        rhs[-1] = 1.0
        sol = np.linalg.solve(mat, rhs)
        coeffs = sol[:-1]
    except np.linalg.LinAlgError:
        return (1.0 - damp) * gamma + damp * gamma_new
    g_stack = np.column_stack([gi.ravel() for gi in hist_g])
    f_stack = F
    mixed = (g_stack + damp * f_stack) @ coeffs
    return mixed.reshape(gamma.shape)


def _cf1_temperature_ladder(t_final: float, anneal_steps: int) -> list[float]:
    """Geometric cool-down from a hot start into the association regime."""
    t_final = float(t_final)
    if anneal_steps <= 1:
        return [t_final]
    # Deep CF1 OH well needs T ≳ 5e3 K before HNC Mayer f-bonds are mild.
    t_hot = max(t_final * 20.0, 6000.0)
    if anneal_steps <= 8:
        return list(np.geomspace(t_hot, t_final, anneal_steps))
    # Two-phase ladder: coarse cool to ~4 T_final, then fine steps to target.
    t_mid = max(4.0 * t_final, 1200.0)
    n_coarse = max(6, anneal_steps // 3)
    n_fine = max(anneal_steps - n_coarse, 8)
    coarse = list(np.geomspace(t_hot, t_mid, n_coarse))
    fine = list(np.geomspace(t_mid, t_final, n_fine + 1))[1:]
    return coarse + fine


def _cf1_stage_controls(temperature: float, mix: float, max_iter: int, tol: float) -> tuple[float, int, float, float]:
    """Return mix, max_iter, tol, exp_clip for one annealing temperature."""
    t = float(temperature)
    if t >= 2000.0:
        return min(0.35, max(mix, 0.2)), min(max_iter, 400), max(tol, 5e-4), 20.0
    if t >= 1000.0:
        return min(0.15, max(mix, 0.1)), min(max_iter, 800), max(tol, 5e-4), 28.0
    if t >= 600.0:
        return min(0.08, max(mix * 1.2, 0.05)), min(max_iter, 1500), max(tol, 2e-4), 35.0
    if t >= 400.0:
        return min(0.05, max(mix, 0.03)), min(max_iter, 2500), max(tol, 1e-4), 45.0
    return max(0.02, min(mix, 0.04)), max_iter, tol, 50.0


def solve_bulk_cf1(
    *,
    temperature: float,
    densities: tuple[float, float],
    charge_pair: np.ndarray,
    n_points: int = 1024,
    delta_r: float = 0.025,
    closure: Closure = "hnc",
    mix: float = 0.05,
    max_iter: int = 4000,
    tol: float = 5e-5,
    ng_alpha: float = 1.5,
    anneal_steps: int = 24,
    use_anderson: bool = True,
    t_hot: float | None = None,
) -> BulkResult:
    """Bulk OZ for CF1 water (O, H) with Ng Coulomb split + temperature annealing.

    Pair order OO, OH, HH. CF1 already contains Coulomb; long-range Coulomb is
    taken from wall-fluid ``charge_pair`` (matches CF1 1/r coeffs to ~0.03%) and
    removed from the iterated potential so Ng renormalisation is not double-counted.

    Direct Picard at 298 K fails because the intramolecular OH Mayer f-bond is
    enormous. A geometric temperature cool-down from ~6000 K with continuous
    ``gamma`` (rescaled by β) converges HNC/Ng to a wall-usable ``c_short`` with
    vanishing large-r tail. Bare HNC under-associates vs experiment/Duh–Haymet
    bridge IE (coord_OH ≲ 2); bridge functions remain future work.
    """
    from scipy.special import erf

    from charged_fluids.potentials.cf1 import calc_u_cf1

    if closure not in ("hnc", "py"):
        raise ValueError(f"CF1 bulk solver supports hnc|py, got {closure!r}")

    rho = np.asarray(densities, dtype=float)
    if rho.size != 2:
        raise ValueError("CF1 densities must be (rho_O, rho_H)")
    charge_pair_final = np.asarray(charge_pair, dtype=float).ravel()
    n_comp = 2
    n_pair = n_pairs(n_comp)
    if charge_pair_final.size != n_pair:
        raise ValueError(f"charge_pair must have length {n_pair}, got {charge_pair_final.size}")

    pair_list = [(0, 0), (0, 1), (1, 1)]
    r, k, delta_k = radial_grids(n_points, delta_r)
    u_full = calc_u_cf1(n_points, r).T  # (n_pair, n_r) joules, OO/OH/HH

    t_final = float(temperature)
    beta_final = beta_from_temperature(t_final)
    temps = _cf1_temperature_ladder(t_final, int(anneal_steps))
    if t_hot is not None and len(temps) > 1:
        temps = list(np.geomspace(float(t_hot), t_final, len(temps)))

    gamma = np.zeros((n_pair, n_points))
    c = np.zeros((n_pair, n_points))
    g = np.zeros((n_pair, n_points))
    converged = False
    residual = np.inf
    n_iter_total = 0
    prev_t: float | None = None

    def _hnc_like(beta_u_s_l: np.ndarray, gamma_l: np.ndarray, exp_clip: float) -> tuple[np.ndarray, np.ndarray]:
        if closure == "hnc":
            bu = np.where(np.isfinite(beta_u_s_l), beta_u_s_l, 1.0e6)
            t = np.clip(-bu + gamma_l, -100.0, exp_clip)
            e = np.exp(t)
            return e - 1.0 - gamma_l, e
        # PY with the same exp clip for association stability
        bu = np.where(np.isfinite(beta_u_s_l), beta_u_s_l, 1.0e6)
        t = np.clip(-bu, -100.0, exp_clip)
        fbond = np.exp(t) - 1.0
        c_l = fbond * (1.0 + gamma_l)
        g_l = (1.0 + fbond) * (1.0 + gamma_l)
        return c_l, g_l

    for i_t, t_step in enumerate(temps):
        t_step = float(t_step)
        beta = beta_from_temperature(t_step)
        scale = beta / beta_final
        cp = charge_pair_final * scale
        erf_term = np.where(
            r > 0, erf(ng_alpha * r) / (2.0 * r), ng_alpha / np.sqrt(np.pi)
        )
        beta_u_long = cp[:, None] * erf_term[None, :]
        # Short-range potential for the Picard unknown (Ng): full CF1 minus long Coulomb
        beta_u_s = beta * u_full - beta_u_long

        if prev_t is not None:
            # Continuity in β when cooling; structural cores + Coulomb γ both track β roughly
            gamma = gamma * (beta / beta_from_temperature(prev_t))
        prev_t = t_step

        mix_stage, stage_iters, stage_tol, exp_clip = _cf1_stage_controls(
            t_step, mix, max_iter, tol
        )
        if i_t == len(temps) - 1:
            stage_iters = max(stage_iters, max_iter)
            stage_tol = tol
            mix_stage = mix

        hist_g: list[np.ndarray] = []
        hist_f: list[np.ndarray] = []
        best: tuple[float, np.ndarray, np.ndarray, np.ndarray] | None = None

        for _iteration in range(1, stage_iters + 1):
            n_iter_total += 1
            for l in range(n_pair):
                c[l], g[l] = _hnc_like(beta_u_s[l], gamma[l], exp_clip)

            c_k = np.zeros((n_comp, n_comp, n_points))
            for i, j in pair_list:
                l = pair_index(i, j)
                ck = forward_sine(c[l], r, k, delta_r)
                c_k[i, j] = ck
                c_k[j, i] = ck

            gamma_k = _binary_oz_gamma_k(c_k, rho)
            gamma_new = np.zeros_like(gamma)
            for i, j in pair_list:
                l = pair_index(i, j)
                gamma_new[l] = inverse_sine(gamma_k[i, j], r, k, delta_k)

            if not np.isfinite(gamma_new).all():
                if best is not None:
                    residual, gamma, c, g = best[0], best[1].copy(), best[2].copy(), best[3].copy()
                converged = False
                break

            residual = float(np.sqrt(np.mean((gamma_new - gamma) ** 2)))
            if best is None or residual < best[0]:
                best = (residual, gamma.copy(), c.copy(), g.copy())

            if use_anderson and tol < residual < 1.0:
                gamma = _anderson_mix(
                    gamma, gamma_new, hist_g, hist_f, m_hist=5, damp=mix_stage
                )
            else:
                if residual > 1.0e3:
                    step = mix_stage * 0.05
                elif residual > 10.0:
                    step = mix_stage * 0.2
                elif residual > 1.0:
                    step = mix_stage * 0.5
                else:
                    step = mix_stage
                gamma = (1.0 - step) * gamma + step * gamma_new
                hist_g.clear()
                hist_f.clear()

            if residual < stage_tol:
                converged = i_t == len(temps) - 1
                break
        else:
            converged = False
            if best is not None:
                residual, gamma, c, g = best[0], best[1].copy(), best[2].copy(), best[3].copy()

        if not np.isfinite(gamma).all():
            converged = False
            break

    # Final closure / c_short at the target temperature
    beta = beta_final
    erf_term = np.where(r > 0, erf(ng_alpha * r) / (2.0 * r), ng_alpha / np.sqrt(np.pi))
    beta_u_c = dimensionless_coulomb_from_charge_pair(r, charge_pair_final)
    beta_u_long = charge_pair_final[:, None] * erf_term[None, :]
    beta_u_short_c = beta_u_c - beta_u_long
    beta_u_s = beta * u_full - beta_u_long
    for l in range(n_pair):
        c[l], g[l] = _hnc_like(beta_u_s[l], gamma[l], 50.0)
        g[l] = np.maximum(g[l], 0.0)

    c_full = c - beta_u_long
    c_short = c + beta_u_short_c  # = c_full + beta_u_c

    def _coord(g_pair: np.ndarray, density: float, r_cut: float) -> float:
        mask = r <= r_cut
        return float(4.0 * np.pi * density * np.trapezoid(g_pair[mask] * r[mask] ** 2, r[mask]))

    n_oh = _coord(g[pair_index(0, 1)], float(rho[1]), 1.3)
    n_hh = _coord(g[pair_index(1, 1)], float(rho[1]), 1.9)
    tail = float(np.max(np.abs(c_short[:, -max(5, n_points // 20) :])))
    i_oh = int(np.argmax(g[pair_index(0, 1)]))

    return BulkResult(
        r=r,
        c_full=c_full,
        c_short=c_short,
        g=g,
        gamma=gamma,
        n_iter=n_iter_total,
        converged=converged,
        residual=residual,
        beta=beta,
        rho=rho,
        meta={
            "fluid": "cf1_ng_anneal",
            "closure": closure,
            "temperature": t_final,
            "densities": [float(rho[0]), float(rho[1])],
            "charge_pair": charge_pair_final.tolist(),
            "n_points": n_points,
            "delta_r": delta_r,
            "ng_alpha": ng_alpha,
            "anneal_steps": len(temps),
            "temperature_ladder": [float(t) for t in temps],
            "use_anderson": use_anderson,
            "coord_OH_1.3A": n_oh,
            "coord_HH_1.9A": n_hh,
            "g_OH_peak": float(g[pair_index(0, 1), i_oh]),
            "g_OH_peak_r": float(r[i_oh]),
            "tail_abs_max_c_short": tail,
            "c_short_definition": "c_s + beta_u_c_short (= c_full + beta_u_coulomb)",
            "note": (
                "CF1 bulk OZ with Ng Coulomb split and temperature annealing; "
                "HNC/PY without bridge (under-associates vs Duh–Haymet). "
                "Aimed at wall-usable c_short with small large-r tail."
            ),
        },
    )
