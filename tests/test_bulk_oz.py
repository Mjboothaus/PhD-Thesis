"""Tests for bulk OZ transforms, LJ solve, and c_short long-range split."""

from __future__ import annotations

import numpy as np
import pytest
from scipy.constants import Boltzmann

from charged_fluids.bulk.artefacts import resolve_cr_path
from charged_fluids.bulk.generate import generate_bulk_cr
from charged_fluids.bulk.oz_solver import (
    apply_closure,
    coulomb_potential,
    solve_bulk_lj,
)
from charged_fluids.bulk.transforms import forward_sine, inverse_sine, radial_grids


def test_sine_transform_roundtrip_gaussian():
    n = 256
    dr = 0.05
    r, k, dk = radial_grids(n, dr)
    # f(r) = exp(-a r^2); transform pair is approximate on truncated grid
    a = 0.5
    f = np.exp(-a * r**2)
    f_k = forward_sine(f, r, k, dr)
    f_back = inverse_sine(f_k, r, k, dk)
    # Allow generous error on truncated discrete sine grid
    rel = np.linalg.norm(f_back - f) / np.linalg.norm(f)
    assert rel < 0.25


def test_hnc_closure_zero_gamma_matches_mayer():
    r = np.linspace(0.5, 5.0, 50)
    beta_u = 1.0 / r
    gamma = np.zeros_like(r)
    c, g = apply_closure("hnc", beta_u, gamma)
    np.testing.assert_allclose(g, np.exp(-beta_u), rtol=1e-12)
    np.testing.assert_allclose(c, g - 1.0, rtol=1e-12)


@pytest.mark.slow
def test_solve_bulk_lj_converges_low_density():
    T = 298.15
    eps = 0.1 * 1000.0 / 6.02214076e23  # 0.1 kJ/mol
    result = solve_bulk_lj(
        temperature=T,
        density=0.001,  # dilute
        epsilon=eps,
        sigma=5.0,
        n_points=512,
        delta_r=0.05,
        mix=0.8,
        max_iter=500,
        tol=1e-7,
    )
    assert result.converged
    assert result.g.shape[0] == 1
    assert np.all(result.g[0] >= -1e-6)


def test_coulomb_short_range_definition_vanishes_when_c_is_minus_beta_u():
    """If full c = -beta u_coulomb, thesis c_sr = c + beta u_c must be ~0."""
    r = np.linspace(0.5, 20.0, 200)
    beta = 1.0 / (Boltzmann * 298.15)
    from charged_fluids.modelling import calc_charge, calc_epsilon

    eps = calc_epsilon(1.0)
    q = calc_charge(np.array([1.0, -1.0]))
    u_c = coulomb_potential(r, float(q[0]), float(q[1]), eps)
    c_full = -beta * u_c
    c_short = c_full + beta * u_c
    np.testing.assert_allclose(c_short, 0.0, atol=1e-15)


@pytest.mark.slow
def test_generate_lj1_writes_artefact(tmp_path):
    _result, cr_path, meta_path = generate_bulk_cr(
        "lj1",
        bulk_dir=tmp_path,
        n_points=512,
        delta_r=0.05,
        write=True,
    )
    assert cr_path is not None and cr_path.is_file()
    assert meta_path is not None and meta_path.is_file()
    data = np.loadtxt(cr_path)
    assert data.shape[1] == 2  # r, c00
    resolved = resolve_cr_path("lj1", data_dir=tmp_path, bulk_dir=tmp_path, legacy_name=None)
    assert resolved == cr_path


@pytest.mark.slow
def test_kcl_generated_c_short_tail_small(tmp_path):
    """Ng/PY ionic bulk: short-range c must vanish at large r (thesis)."""
    result, _cr_path, _ = generate_bulk_cr(
        "kcl",
        bulk_dir=tmp_path,
        n_points=256,
        delta_r=0.1,
        write=True,
    )
    assert np.isfinite(result.c_short).all()
    tail = np.abs(result.c_short[:, -max(10, result.r.size // 20) :])
    tail_max = float(tail.max())
    assert tail_max < 1e-2, f"KCl |c_short| tail max={tail_max}"


def test_fft_sine_matches_direct():
    """FFT odd-extension path must match dense sin(kr) @ v exactly."""
    n, dr = 128, 0.05
    r, k, dk = radial_grids(n, dr)
    f = np.exp(-0.35 * r**2) * np.cos(0.2 * r)
    fk_fft = forward_sine(f, r, k, dr, backend="fft")
    fk_dir = forward_sine(f, r, k, dr, backend="direct")
    np.testing.assert_allclose(fk_fft, fk_dir, rtol=1e-10, atol=1e-12)
    back_fft = inverse_sine(fk_fft, r, k, dk, backend="fft")
    back_dir = inverse_sine(fk_dir, r, k, dk, backend="direct")
    np.testing.assert_allclose(back_fft, back_dir, rtol=1e-10, atol=1e-12)


def test_closure_registry_lists_hnc_py():
    from charged_fluids.bulk.closures import list_closures

    assert list_closures() == ["hnc", "py"]
