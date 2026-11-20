"""TOML parameter loading and shared setup orchestration."""

from __future__ import annotations

import numpy as np
import pytest

from charged_fluids import (
    fluid_specific_parameters,
    list_fluid_symbols,
    load_fluid_table,
    prepare_fluid_run,
    set_fluid_parameters,
)


def test_toml_lists_expected_fluids():
    symbols = list_fluid_symbols()
    assert "lj1" in symbols
    assert "kcl" in symbols
    assert "h2o" in symbols


def test_lj1_from_toml_matches_historical_values():
    fluid = set_fluid_parameters("lj1")
    assert fluid is not None
    assert fluid.name.startswith("Lennard-Jones")
    assert fluid.n_component == 1
    assert fluid.n_pair == 1
    assert fluid.cr_filename == "lj1.cr.dat"
    assert fluid.status == "ready"
    pot = fluid_specific_parameters("lj1")["lj1"]
    np.testing.assert_allclose(pot["epsilon_lj"], [0.1])
    np.testing.assert_allclose(pot["sigma_lj"], [0.5])


def test_kcl_potential_from_toml():
    pot = fluid_specific_parameters("kcl")["kcl"]
    assert pot["type"] == "ionic"
    assert float(pot["b"]) == pytest.approx(0.338e-19)
    assert float(pot["alpha"]) == pytest.approx(1.0 / 0.337)


def test_h2o_cf1_potential_block():
    fluid = set_fluid_parameters("h2o")
    assert fluid is not None
    assert fluid.status == "experimental"
    pot = fluid_specific_parameters("h2o")["h2o"]
    assert pot["type"] == "cf1"


def test_prepare_lj1_run_small_grid(data_dir):
    run = prepare_fluid_run(
        "lj1",
        z_cutoff=8.0,
        n_point=64,
        tolerance=1e-4,
        max_iteration=20,
        data_dir=data_dir,
    )
    assert run.model.tw.shape == (64, 1)
    assert run.model.hw.shape == (64, 1)
    assert run.model.c_short.shape == (64, 1)
    assert run.model.f1.shape == run.model.f2.shape == (64, 1)
    assert np.isfinite(run.beta_u[1:]).all()
    assert run.fluid.beta > 0
    assert run.cr_path.name == "lj1.cr.dat"


def test_load_toml_cache_is_path_keyed(tmp_path):
    # Ensure custom path works and does not collide with default cache key
    src = load_fluid_table()
    custom = tmp_path / "mini.toml"
    custom.write_text(
        '[x]\nname="X"\nsymbol="x"\ncomponent=["A"]\nvalence=[0.0]\n'
        "temperature=300.0\nconcentration=[1.0]\nepsilon_r=1.0\nindex=0\n"
        'cr_filename="none"\nstatus="experimental"\n'
    )
    table = load_fluid_table(custom)
    assert "x" in table
    assert "lj1" in src
