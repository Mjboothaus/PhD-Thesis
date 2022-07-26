import pytest
from src.modelling import *

def test_calc_beta():
    temperature = 1.0/Boltzmann
    assert calc_beta(temperature) == pytest.approx(1.0, rel=1e-16)


# def test_calc_epsilon(epsilon_r):
    # units same as $\epsilon_0$ (need to allow for distances in angstrom)
    # return 4.0 * np.pi * epsilon_r * epsilon_0


def test_calc_l_index():
    n_component = 2
    l_index = [0, 1, 2]
    for i in range(n_component):
        for j in range(i, n_component):
            l = calc_l_index(i, j)
            assert l == l_index[l]


def test_kappa_kcl():
    kcl = set_fluid_parameters("kcl")
    charge = calc_charge(kcl.valence)
    kcl.beta = calc_beta(kcl.temperature)
    kcl.epsilon = calc_epsilon(kcl.epsilon_r)
    kcl.rho = calc_rho(kcl.concentration)
    assert calc_kappa(kcl.beta, charge, kcl.rho, kcl.epsilon) == pytest.approx(6.77, abs=0.04)