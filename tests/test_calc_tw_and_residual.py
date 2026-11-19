"""Characterisation tests for calc_tw and the Newton residual."""

from __future__ import annotations

import numpy as np

from charged_fluids import calc_hw, calc_tw, opt_func
from charged_fluids.memory import get_memory_usage
from charged_fluids.modelling import interpolate_cr
from charged_fluids.parameters import fluid_specific_parameters, set_fluid_parameters


def _neutral_one_component_setup(n_point: int = 16, z_cutoff: float = 4.0):
    z = np.linspace(0.0, z_cutoff, n_point)
    z_index = np.arange(n_point, dtype=int)
    n_component = 1
    n_pair = 1
    # Mild wall field so exp(tw - beta_phiw) stays well behaved
    beta_phiw = np.column_stack([0.05 * z])
    beta_psi_charge = np.zeros(n_component)
    charge_pair = np.zeros(n_pair)
    rho = np.array([0.01])
    f1 = np.zeros((n_point, n_pair))
    f2 = np.zeros((n_point, n_pair))
    return {
        "z": z,
        "z_index": z_index,
        "n_point": n_point,
        "n_component": n_component,
        "beta_phiw": beta_phiw,
        "beta_psi_charge": beta_psi_charge,
        "charge_pair": charge_pair,
        "rho": rho,
        "f1": f1,
        "f2": f2,
    }


def test_calc_tw_zero_inputs_stays_finite():
    s = _neutral_one_component_setup()
    tw_in = np.zeros((s["n_point"], s["n_component"]))
    tw = calc_tw(
        tw_in,
        s["beta_phiw"],
        s["beta_psi_charge"],
        s["charge_pair"],
        s["rho"],
        s["f1"],
        s["f2"],
        s["z"],
        s["n_component"],
        s["n_point"],
        s["z_index"],
    )
    assert tw.shape == tw_in.shape
    assert np.isfinite(tw).all()


def test_opt_func_is_tw_in_minus_calc_tw():
    s = _neutral_one_component_setup()
    tw_in = np.full((s["n_point"], s["n_component"]), 0.1)
    residual = opt_func(
        tw_in,
        s["beta_phiw"],
        s["beta_psi_charge"],
        s["charge_pair"],
        s["rho"],
        s["f1"],
        s["f2"],
        s["z"],
        s["n_component"],
        s["n_point"],
        s["z_index"],
    )
    tw = calc_tw(
        tw_in,
        s["beta_phiw"],
        s["beta_psi_charge"],
        s["charge_pair"],
        s["rho"],
        s["f1"],
        s["f2"],
        s["z"],
        s["n_component"],
        s["n_point"],
        s["z_index"],
    )
    np.testing.assert_allclose(residual, tw_in - tw)


def test_calc_tw_with_nonzero_f_kernels_changes_output():
    s = _neutral_one_component_setup()
    tw_in = np.zeros((s["n_point"], s["n_component"]))
    baseline = calc_tw(
        tw_in,
        s["beta_phiw"],
        s["beta_psi_charge"],
        s["charge_pair"],
        s["rho"],
        s["f1"],
        s["f2"],
        s["z"],
        s["n_component"],
        s["n_point"],
        s["z_index"],
    )
    f1 = np.full_like(s["f1"], 0.2)
    f2 = np.full_like(s["f2"], 0.1)
    perturbed = calc_tw(
        tw_in,
        s["beta_phiw"],
        s["beta_psi_charge"],
        s["charge_pair"],
        s["rho"],
        f1,
        f2,
        s["z"],
        s["n_component"],
        s["n_point"],
        s["z_index"],
    )
    assert not np.allclose(baseline, perturbed)
    assert np.isfinite(perturbed).all()


def test_calc_hw_matches_definition_for_nonzero_tw():
    n_point = 8
    z = np.linspace(0.0, 1.0, n_point)
    beta_phiw = np.column_stack([z, 2 * z])
    tw = np.column_stack([0.5 * z, -0.25 * z])
    hw = calc_hw(tw, n_component=2, beta_phiw=beta_phiw)
    np.testing.assert_allclose(hw, np.exp(tw - beta_phiw) - 1.0)


def test_interpolate_cr_accepts_1d_cr_column():
    z = np.linspace(0.0, 5.0, 11)
    r_in = np.linspace(0.0, 5.0, 6)
    cr_1d = -0.25 * np.ones(6)  # shape (n_r,) not (n_r, 1)
    cr, z_out = interpolate_cr(r_in, cr_1d, n_point=len(z), n_pair=1, z=z)
    assert z_out.shape == z.shape
    np.testing.assert_allclose(cr, -0.25)


def test_get_memory_usage_positive():
    mem = get_memory_usage()
    assert mem > 0.0


def test_lj2_fluid_specific_parameters_shape():
    params = fluid_specific_parameters("lj2")["lj2"]
    assert params["epsilon_lj"].shape == (3,)
    assert params["sigma_lj"].shape == (3,)
    fluid = set_fluid_parameters("lj2")
    assert fluid is not None
    assert fluid.n_component == 2
    assert fluid.n_pair == 3


def test_calc_tw_two_component_neutral_runs():
    n_point = 12
    z = np.linspace(0.0, 3.0, n_point)
    z_index = np.arange(n_point, dtype=int)
    n_component = 2
    n_pair = 3
    beta_phiw = np.column_stack([0.02 * z, 0.03 * z])
    tw_in = np.zeros((n_point, n_component))
    tw = calc_tw(
        tw_in,
        beta_phiw,
        np.zeros(n_component),
        np.zeros(n_pair),
        np.array([0.005, 0.005]),
        np.zeros((n_point, n_pair)),
        np.zeros((n_point, n_pair)),
        z,
        n_component,
        n_point,
        z_index,
    )
    assert tw.shape == (n_point, n_component)
    assert np.isfinite(tw).all()
