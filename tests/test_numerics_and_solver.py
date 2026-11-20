"""Discretisation helpers and optional LJ solver smoke test."""

from __future__ import annotations

import numpy as np
import pytest

from charged_fluids import (
    Model,
    calc_hw,
    calc_phiw,
    calc_u_lj,
    integral_z_infty_dr_r2_c_short,
    integral_z_infty_dr_r_c_short,
    load_and_interpolate_cr,
    opt_func,
    set_num_parameters,
    solve_model,
)


def test_set_num_parameters_grid_size():
    d = set_num_parameters(
        n_point=101,
        z_cutoff=10.0,
        n_component=1,
        n_pair=1,
        tolerance=1e-6,
        max_iteration=100,
    )
    assert d.grid_size == pytest.approx(0.1)
    assert d.n_point == 101


@pytest.mark.slow
def test_lj1_solver_smoke_converges(lj1_fluid, lj1_params, lj1_cr_path, small_grid, lj1_discrete):
    """Small-grid LJ wall solve should converge with relaxed tolerance."""
    n_point, _, z, z_index = small_grid
    n_component = lj1_fluid.n_component
    n_pair = lj1_fluid.n_pair

    c_short, _ = load_and_interpolate_cr(lj1_cr_path, n_point, n_pair, z)
    f1 = integral_z_infty_dr_r_c_short(c_short, n_pair, n_point, z)
    f2 = integral_z_infty_dr_r2_c_short(c_short, n_pair, n_point, z)

    wall_zeros = np.zeros((n_point, n_component))
    model = Model(
        z=z,
        z_index=z_index,
        hw=wall_zeros,
        tw=wall_zeros.copy(),
        c_short=c_short,
        f1=f1,
        f2=f2,
    )
    # Touch LJ potential to ensure fixture wiring (not fed to residual directly)
    _ = calc_u_lj(
        lj1_params["epsilon_lj"],
        lj1_params["sigma_lj"],
        n_point,
        n_component,
        n_pair,
        z,
    )

    beta_phiw = lj1_fluid.beta * calc_phiw(z, n_point, n_component)
    beta_psi_charge = np.zeros(n_component)
    tw_initial = np.zeros((n_point, n_component))

    solution = solve_model(
        opt_func,
        tw_initial,
        lj1_fluid,
        model,
        lj1_discrete,
        beta_phiw,
        beta_psi_charge,
        track_history=True,
        disp=False,
    )
    assert solution.success, f"LJ smoke solve failed: {solution.message}"
    hist = getattr(solution, "convergence_history", None)
    assert hist is not None and len(hist.iterations) >= 1
    assert hist.residuals[-1] < hist.residuals[0] or hist.residuals[-1] < 1e-3
    hw = calc_hw(solution.x, n_component, beta_phiw)
    g = hw + 1.0
    assert np.isfinite(g).all()
    # g(z) >= 0; contact value may be ~0 for a strongly repulsive wall
    assert (g >= 0).all()
    # Far from the wall, g should be nearer bulk (order-1)
    assert g[-1, 0] == pytest.approx(1.0, abs=0.5)


def test_convergence_history_update_unit():
    from charged_fluids import ConvergenceHistory

    h = ConvergenceHistory()
    h.update(1, np.array([3.0, 4.0]), step=0.5)
    assert h.residuals[0] == pytest.approx(5.0)
    assert h.max_abs[0] == pytest.approx(4.0)
    assert h.summary()["final_step_size"] == pytest.approx(0.5)
