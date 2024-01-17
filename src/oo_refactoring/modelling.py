from dataclasses import dataclass

import numpy as np
import scipy.optimize as optim
from scipy import interpolate
from scipy.constants import Avogadro, Boltzmann, elementary_charge, epsilon_0
from scipy.integrate import trapezoid


@dataclass
class Model:
    z: np.array
    z_index: np.array
    hw: np.array
    c_short: np.array
    f1: np.array
    f2: np.array


class PotentialCalculations:
    @staticmethod
    def calc_beta(temperature):
        return 1.0 / (Boltzmann * temperature)

    # ... other static methods for calculations ...


class Solver:
    @staticmethod
    def solve_model(
        opt_func, tw_initial, fluid, model, discrete, beta_phiw, beta_psi_charge
    ):
        charge_pair = fluid.charge_pair
        n_component = fluid.n_component
        rho = fluid.rho
        f1 = model.f1
        f2 = model.f2
        z = model.z
        z_index = model.z_index
        n_point = discrete.n_point
        tolerance = discrete.tolerance
        max_iteration = discrete.max_iteration
        tw_args = (
            beta_phiw,
            beta_psi_charge,
            charge_pair,
            rho,
            f1,
            f2,
            z,
            n_component,
            n_point,
            z_index,
        )

        return optim.root(
            opt_func,
            tw_initial,
            args=tw_args,
            method="krylov",
            jac=None,
            callback=None,
            options={"disp": True, "maxiter": max_iteration, "fatol": tolerance},
        )
