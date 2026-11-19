"""Tests for bulk c(r) I/O and short-range integral transforms."""

from __future__ import annotations

import numpy as np
import pytest

from charged_fluids import (
    calc_f1_integrand,
    calc_f2_integrand,
    calc_hw,
    integral_z_infty_dr_r2_c_short,
    integral_z_infty_dr_r_c_short,
    interpolate_cr,
    load_and_interpolate_cr,
)


def test_load_lj1_cr_shape_and_finite(lj1_cr_path, small_grid):
    n_point, _, z, _ = small_grid
    c_short, z_out = load_and_interpolate_cr(lj1_cr_path, n_point, n_pair=1, z=z)
    assert z_out.shape == (n_point,)
    assert c_short.shape == (n_point, 1)
    assert np.isfinite(c_short).all()


def test_interpolate_cr_constant_field(small_grid):
    n_point, _, z, _ = small_grid
    r_in = np.linspace(0.0, 10.0, 21)
    cr_in = np.full((21, 1), -0.5)
    cr, _ = interpolate_cr(r_in, cr_in, n_point, n_pair=1, z=z)
    np.testing.assert_allclose(cr, -0.5)


def test_f_integrands_match_definition(small_grid):
    n_point, _, z, _ = small_grid
    c_short = np.ones((n_point, 1))
    f1 = calc_f1_integrand(c_short, 1, z, n_point)
    f2 = calc_f2_integrand(c_short, 1, z, n_point)
    np.testing.assert_allclose(f1[:, 0], z)
    np.testing.assert_allclose(f2[:, 0], z * z)


def test_integral_f1_zero_for_zero_c_short(small_grid):
    n_point, _, z, _ = small_grid
    c_short = np.zeros((n_point, 1))
    f1 = integral_z_infty_dr_r_c_short(c_short, 1, n_point, z)
    f2 = integral_z_infty_dr_r2_c_short(c_short, 1, n_point, z)
    np.testing.assert_allclose(f1, 0.0)
    np.testing.assert_allclose(f2, 0.0)


def test_integral_f1_analytic_constant_c(small_grid):
    """For c=const and last point truncated, f1[k] ≈ c * ∫_z^zmax z' dz'."""
    n_point, z_cutoff, z, _ = small_grid
    c_val = 2.0
    c_short = np.full((n_point, 1), c_val)
    f1 = integral_z_infty_dr_r_c_short(c_short, 1, n_point, z)
    # Compare interior points away from the open endpoint
    for k in (0, n_point // 4, n_point // 2):
        expected = c_val * 0.5 * (z_cutoff**2 - z[k] ** 2)
        assert f1[k, 0] == pytest.approx(expected, rel=5e-2, abs=1e-3)


def test_calc_hw_zero_tw_equals_boltzmann_factor(small_grid):
    _n_point, _, z, _ = small_grid
    beta_phiw = np.column_stack([0.1 * z, 0.2 * z])
    tw = np.zeros_like(beta_phiw)
    hw = calc_hw(tw, n_component=2, beta_phiw=beta_phiw)
    np.testing.assert_allclose(hw, np.exp(-beta_phiw) - 1.0)


def test_kcl_bulk_artefact_short_range_vanishes_at_large_r(kcl_cr_path):
    """Shipped ``data/bulk/kcl.cr.dat`` must have a near-zero short-range tail."""
    data = np.loadtxt(kcl_cr_path)
    assert data.ndim == 2 and data.shape[1] >= 2
    r = data[:, 0]
    cr = data[:, 1:]
    assert np.isfinite(cr).all()
    tail = cr[int(0.95 * len(r)) :]
    tail_max = float(np.max(np.abs(tail)))
    assert tail_max < 1e-2, (
        f"KCl short-range |c(r)| tail max={tail_max:.4g}; expected near zero"
    )
