"""Wall-solve smoke tests across ready and experimental fluids."""

from __future__ import annotations

import numpy as np
import pytest

from charged_fluids import (
    calc_hw,
    clear_fluid_cache,
    opt_func,
    prepare_fluid_run,
    solve_model,
    wall_ready_symbols,
)


@pytest.fixture(scope="module", autouse=True)
def _reload_toml():
    clear_fluid_cache()


def test_wall_ready_symbols_include_new_fluids():
    symbols = wall_ready_symbols()
    for s in ("lj1", "lj2", "kcl", "h2o", "2_2", "1_1"):
        assert s in symbols


@pytest.mark.parametrize("symbol", ["lj1", "1_1", "2_2", "h2o"])
def test_prepare_and_wall_smoke(symbol, data_dir):
    run = prepare_fluid_run(
        symbol,
        z_cutoff=10.0,
        n_point=81,
        tolerance=1e-5,
        max_iteration=25,
        psi_0_mv=0.0,
        data_dir=data_dir,
    )
    assert run.model.c_short.ndim == 2
    assert run.model.c_short.shape[0] == 81
    assert np.isfinite(run.beta_u[1:]).all()

    if run.fluid.status not in {"ready", "experimental"}:
        pytest.skip(f"{symbol} not wall-attemptable")

    sol = solve_model(
        opt_func,
        run.model.tw.copy(),
        run.fluid,
        run.model,
        run.discrete,
        run.beta_phiw,
        run.beta_psi_charge,
    )
    # Accept convergence or a finite residual trajectory (experimental melts/water)
    assert np.isfinite(sol.x).all()
    g = calc_hw(sol.x, run.fluid.n_component, run.beta_phiw) + 1.0
    assert np.isfinite(g).all()
    assert g.shape == (81, run.fluid.n_component)
