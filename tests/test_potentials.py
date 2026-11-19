"""Tests for pair potentials and the wall potential."""

from __future__ import annotations

import numpy as np
import pytest

from charged_fluids import (
    calc_beta_pauling,
    calc_cap_b,
    calc_phiw,
    calc_u,
    calc_u_lj,
    fluid_specific_parameters,
)


def test_lj_potential_minimum_near_sigma_times_2_to_1_6(lj1_params, small_grid):
    _, _, r, _ = small_grid
    n_point = len(r)
    u = calc_u_lj(
        lj1_params["epsilon_lj"],
        lj1_params["sigma_lj"],
        n_point,
        n_component=1,
        n_pair=1,
        r=r,
    )
    assert np.isinf(u[0, 0])
    # Finite for r > 0
    assert np.isfinite(u[1:, 0]).all()
    # Classic LJ minimum at r = 2^(1/6) sigma
    sigma = float(lj1_params["sigma_lj"][0])
    eps = float(lj1_params["epsilon_lj"][0])
    r_min = sigma * 2.0 ** (1.0 / 6.0)
    # Evaluate on a fine local grid for the analytic check
    r_fine = np.linspace(0.2, 2.0, 400)
    u_fine = calc_u_lj(np.array([eps]), np.array([sigma]), len(r_fine), 1, 1, r_fine)
    imin = int(np.nanargmin(u_fine[1:, 0])) + 1
    assert r_fine[imin] == pytest.approx(r_min, rel=2e-2)
    assert u_fine[imin, 0] == pytest.approx(-eps, rel=2e-2)


def test_phiw_positive_and_decays(small_grid):
    n_point, _, z, _ = small_grid
    phiw = calc_phiw(z, n_point, n_component=1)
    assert phiw.shape == (n_point, 1)
    assert (phiw[:, 0] > 0).all()
    assert phiw[-1, 0] < phiw[0, 0]


def test_kcl_pauling_and_cap_b_shapes():
    params = fluid_specific_parameters("kcl")["kcl"]
    valence = np.array([1.0, -1.0])
    n_component = 2
    n_pair = 3
    beta_p = calc_beta_pauling(valence, params["n_outer_shell"], n_component, n_pair)
    assert beta_p.shape == (n_pair,)
    assert beta_p[0] == pytest.approx(1.0 + 1.0 / 8.0 + 1.0 / 8.0)
    cap_b = calc_cap_b(beta_p, params["b"], params["alpha"], params["sigma"], n_component, n_pair)
    assert cap_b.shape == (n_pair,)
    assert (cap_b > 0).all()


def test_ionic_potential_finite_away_from_origin(kcl_fluid, small_grid):
    params = fluid_specific_parameters("kcl")["kcl"]
    n_point, _, r, _ = small_grid
    n_component = kcl_fluid.n_component
    n_pair = kcl_fluid.n_pair
    beta_p = calc_beta_pauling(kcl_fluid.valence, params["n_outer_shell"], n_component, n_pair)
    cap_b = calc_cap_b(beta_p, params["b"], params["alpha"], params["sigma"], n_component, n_pair)
    u = calc_u(
        kcl_fluid.charge,
        cap_b,
        params["alpha"],
        params["cap_c"],
        params["cap_d"],
        n_point,
        n_component,
        n_pair,
        kcl_fluid.epsilon,
        r,
    )
    assert u.shape == (n_point, n_pair)
    assert np.isfinite(u[1:]).all()


def test_fluid_specific_unknown_raises():
    with pytest.raises(KeyError):
        fluid_specific_parameters("h2o")
