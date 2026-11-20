"""Unit tests for thermodynamic and pair-index helpers."""

from __future__ import annotations

import numpy as np
import pytest
from scipy.constants import Avogadro, Boltzmann, elementary_charge, epsilon_0

from charged_fluids import (
    calc_beta,
    calc_charge,
    calc_epsilon,
    calc_kappa,
    calc_l_index,
    calc_rho,
    set_fluid_parameters,
)


def test_calc_beta_inverse_boltzmann():
    temperature = 1.0 / Boltzmann
    assert calc_beta(temperature) == pytest.approx(1.0, rel=1e-15)


def test_calc_epsilon_scales_vacuum_permittivity():
    assert calc_epsilon(1.0) == pytest.approx(4.0 * np.pi * epsilon_0)
    assert calc_epsilon(78.3) == pytest.approx(4.0 * np.pi * 78.3 * epsilon_0)


def test_calc_rho_mol_per_dm3_to_per_angstrom3():
    conc = np.array([1.0])
    expected = conc * Avogadro / 1.0e27
    np.testing.assert_allclose(calc_rho(conc), expected)


def test_calc_charge_uses_elementary_charge():
    valence = np.array([1.0, -1.0])
    np.testing.assert_allclose(calc_charge(valence), valence * elementary_charge)


def test_calc_l_index_one_and_two_component():
    assert calc_l_index(0, 0) == 0
    # Historical two-component upper triangle: (0,0)->0, (0,1)->1, (1,1)->2
    pairs = [(0, 0), (0, 1), (1, 1)]
    indices = [calc_l_index(i, j) for i, j in pairs]
    assert indices == [0, 1, 2]


def test_set_fluid_parameters_unknown_returns_none():
    assert set_fluid_parameters("not-a-fluid") is None


def test_set_fluid_parameters_lj1_shape():
    fluid = set_fluid_parameters("lj1")
    assert fluid is not None
    assert fluid.n_component == 1
    assert fluid.n_pair == 1
    assert fluid.cr_filename == "lj1.cr.dat"


def test_set_fluid_parameters_does_not_mutate_global_table():
    first = set_fluid_parameters("lj1")
    second = set_fluid_parameters("lj1")
    assert first is not None and second is not None
    first.temperature = 1.0
    assert second.temperature == pytest.approx(298.15)


def test_kappa_kcl_matches_historical_characterisation(kcl_fluid):
    kappa = calc_kappa(kcl_fluid.beta, kcl_fluid.charge, kcl_fluid.rho, kcl_fluid.epsilon)
    assert kappa == pytest.approx(6.77, abs=0.04)


def test_charge_pair_neutral_for_uncharged_lj(lj1_fluid):
    np.testing.assert_allclose(lj1_fluid.charge_pair, 0.0)


def test_calc_l_index_three_component_no_collision():
    """For n>2, pair indices must be unique across the upper triangle."""
    pairs = [(i, j) for i in range(3) for j in range(i, 3)]
    indices = [calc_l_index(i, j) for i, j in pairs]
    assert len(indices) == len(set(indices))
    assert sorted(indices) == list(range(6))
    # n<=2 historical layout unchanged
    assert [calc_l_index(i, j) for i, j in [(0, 0), (0, 1), (1, 1)]] == [0, 1, 2]
