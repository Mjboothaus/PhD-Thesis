"""Shared fixtures for charged_fluids unit tests."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from charged_fluids import (
    calc_beta,
    calc_charge,
    calc_charge_pair,
    calc_epsilon,
    calc_rho,
    fluid_specific_parameters,
    set_fluid_parameters,
    set_num_parameters,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT / "data"


@pytest.fixture(scope="session")
def data_dir() -> Path:
    return DATA_DIR


@pytest.fixture(scope="session")
def lj1_cr_path(data_dir: Path) -> Path:
    path = data_dir / "lj1-cr.data"
    assert path.is_file(), f"Missing bulk c(r) fixture: {path}"
    return path


@pytest.fixture(scope="session")
def kcl_cr_path(data_dir: Path) -> Path:
    """Modern Ng/PY short-range KCl artefact (not the legacy nrcg dump)."""
    path = data_dir / "bulk" / "kcl.cr.dat"
    assert path.is_file(), f"Missing bulk c(r) fixture: {path}"
    return path


@pytest.fixture
def lj1_fluid():
    fluid = set_fluid_parameters("lj1")
    assert fluid is not None
    fluid.beta = calc_beta(fluid.temperature)
    fluid.epsilon = calc_epsilon(fluid.epsilon_r)
    fluid.charge = calc_charge(fluid.valence)
    fluid.charge_pair = calc_charge_pair(
        fluid.beta, fluid.charge, fluid.epsilon, fluid.n_component, fluid.n_pair
    )
    fluid.rho = calc_rho(fluid.concentration)
    return fluid


@pytest.fixture
def kcl_fluid():
    fluid = set_fluid_parameters("kcl")
    assert fluid is not None
    fluid.beta = calc_beta(fluid.temperature)
    fluid.epsilon = calc_epsilon(fluid.epsilon_r)
    fluid.charge = calc_charge(fluid.valence)
    fluid.charge_pair = calc_charge_pair(
        fluid.beta, fluid.charge, fluid.epsilon, fluid.n_component, fluid.n_pair
    )
    fluid.rho = calc_rho(fluid.concentration)
    return fluid


@pytest.fixture
def lj1_params():
    return fluid_specific_parameters("lj1")["lj1"]


@pytest.fixture
def small_grid():
    """Modest grid suitable for fast unit tests (not production resolution)."""
    n_point = 64
    z_cutoff = 8.0
    z = np.linspace(0.0, z_cutoff, n_point)
    z_index = np.arange(n_point, dtype=int)
    return n_point, z_cutoff, z, z_index


@pytest.fixture
def lj1_discrete(lj1_fluid, small_grid):
    n_point, z_cutoff, _, _ = small_grid
    return set_num_parameters(
        n_point=n_point,
        z_cutoff=z_cutoff,
        n_component=lj1_fluid.n_component,
        n_pair=lj1_fluid.n_pair,
        tolerance=1e-4,
        max_iteration=50,
    )
