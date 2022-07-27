
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import interpolate
import scipy.optimize as optim
from scipy.constants import Avogadro, Boltzmann, elementary_charge, epsilon_0
from scipy.integrate import trapezoid
from streamlit import cache

# from helper_functions import read_render_markdown_file


@dataclass
class Model:
    z: np.array
    z_index: np.array
    hw: np.array
    c_short: np.array
    f1: np.array
    f2: np.array


@dataclass
class Fluid:
    name: str
    symbol: str
    component: list[str]
    valence: np.ndarray
    charge: np.ndarray
    temperature: float
    concentration: np.ndarray
    epsilon_r: float
    n_component: int
    n_pair: int
    index: int
    charge_pair: np.array
    rho: np.array
    beta: float
    epsilon: float

# Note: charge & rho need to be calculate and then inserted into the Fluid instance
#       np.array([0]) is a placeholder (which has the right type)


fluid_parameters = dict({"kcl": dict({"name": "Potassium Chloride", "component": ["K", "Cl"], "valence": np.array([
    1.0, -1.0]), "charge": np.array([0]), "temperature": 1075.0, "concentration": np.array([19.5, 19.5]),
    "epsilon_r": 1.0, "index": 1, "charge_pair": np.array([0]), "rho": np.array([0]), "beta": 0.0, "epsilon": 0.0})})

fluid_parameters["h2o"] = dict({"name": "Liquid water", "component": ["H", "2O"], "valence": np.array([
    1.0, -1.0]), "temperature": 298.0, "concentration": np.array([1.0, 1.0]),
    "epsilon_r": 1.0, "index": 2})

fluid_parameters["2_2"] = dict({"name": "2-2 Aqueous electrolyte", "component": ["+2", "-2"], "valence": np.array([
    2.0, -2.0]), "temperature": 298.0, "concentration": np.array([1.0, 1.0]),
    "epsilon_r": 78.3, "index": 3})


def set_fluid_parameters(symbol):
    symbol = symbol.lower()
    if symbol not in fluid_parameters:
        return None
    parameters = fluid_parameters[symbol]
    parameters["symbol"] = "".join(parameters["component"])
    n_component = parameters["n_component"] = len(parameters["component"])
    parameters["n_pair"] = int((n_component+1) * (n_component) / 2)
    return Fluid(name=parameters["name"], symbol=parameters["symbol"], component=parameters["component"], valence=parameters["valence"], charge=parameters["charge"], temperature=parameters["temperature"], concentration=parameters["concentration"], epsilon_r=parameters["epsilon_r"], n_component=parameters["n_component"], n_pair=parameters["n_pair"], index=parameters["index"], charge_pair=parameters["charge_pair"], rho=parameters["rho"], beta=parameters["beta"], epsilon=parameters["epsilon"])


def fluid_specific_parameters(symbol):
    if symbol == "kcl":
        other_params = dict({"kcl": dict({"n_outer_shell": np.array([8., 8.]), "alpha": 1.0 / 0.337,
                                          "b": 0.338e-19, "sigma": [1.463, 1.585],
                                          "cap_c": np.array([24.3, 48.0, 124.5]) * 1e-19,
                                          "cap_d": np.array([24.0, 73.0, 250.0]) * 1e-19})})

    return other_params


def calc_beta(temperature):
    return 1.0 / (Boltzmann * temperature)


def calc_epsilon(epsilon_r):
    # units same as $\epsilon_0$ (need to allow for distances in angstrom)
    return 4.0 * np.pi * epsilon_r * epsilon_0


def calc_l_index(i, j):
    return i + j


def calc_beta_pauling(valence, n_outer_shell, n_component, n_pair):
    beta_pauling = np.zeros(n_pair)
    for i in range(n_component):
        for j in range(i, n_component):
            l = calc_l_index(i, j)
            beta_pauling[l] = 1.0 + valence[i] / \
                n_outer_shell[i] + valence[j] / n_outer_shell[j]
    return beta_pauling


def calc_cap_b(beta_pauling, b, alpha, sigma, n_component, n_pair):
    cap_b = np.zeros(n_pair)
    for i in range(n_component):
        for j in range(i, n_component):
            l = calc_l_index(i, j)
            cap_b[l] = beta_pauling[l] * b * \
                np.exp(alpha * (sigma[i] + sigma[j]))
    return cap_b


def calc_charge(valence):
    return valence * elementary_charge


def calc_charge_pair(beta, charge, epsilon, n_component, n_pair):
    charge_pair = np.zeros(n_pair)
    for i in range(n_component):
        for j in range(i, n_component):
            l = calc_l_index(i, j)
            charge_pair[l] = (2.0e10 * beta * charge[i] * charge[j] / epsilon)
    return charge_pair


def calc_u(charge, cap_b, alpha, cap_c, cap_d, n_point, n_component, n_pair, epsilon, r):
    u = np.zeros((n_point, n_pair))
    for i in range(n_component):
        for j in range(i, n_component):
            l = calc_l_index(i, j)
            u[1:, l] = (charge[i] * charge[j]) / (r[1:] * 1e-10 * epsilon) + cap_c[l] / \
                r[1:]**6 + cap_d[l]/r[1:]**8 + \
                cap_b[l] * np.exp(-alpha * r[1:])
            u[0, l] = u[1, l]
    return u


def calc_rho(concentration):
    return np.array(concentration) / 1.0e27 * Avogadro


def calc_kappa(beta, charge, rho, epsilon):
    return np.sqrt(4.0 * np.pi * beta / epsilon * 1e10 *
                   sum(np.multiply(charge**2, rho)))


def calc_phiw(z, n_point, n_component):
    phiw = np.zeros((n_point, n_component))
    capital_a = 16.274e-19  # joules
    wall_d = 2.97  # inverse Angstrom
    for i in range(n_component):
        phiw[:, i] = np.exp(-wall_d * z) * capital_a * (wall_d * z + 2)
    return phiw

# TODO: Fix up choice of c(r) - in progress
# Read in some c(r) -- currently LJ fluid (assuming c_short ~= c_LJ(r))


def interpolate_cr(r_in, cr_in, n_point, n_pair, z):
    cr = np.zeros((n_point, n_pair))
    for l in range(n_pair):
        f = interpolate.interp1d(r_in, cr_in[:, l])
        r = z
        cr[:, l] = f(r)
        # TODO: Make general - were getting some kinks in c(r) near r=0
        # cr[:10, l] = cr[10, l]
    return cr, r


@cache
def load_and_interpolate_cr(cr_path, n_point, n_pair, z):
    # CR_PATH = "../pyOZ_bulk_fluid/tests/lj/nrcg-cr.dat.orig"
    cr_df = pd.read_csv(cr_path, header=None, delim_whitespace=True)
    cr_df.set_index(0, inplace=True)
    r = cr_df.index.to_numpy()
    cr = cr_df.to_numpy()
    return interpolate_cr(r, cr, n_point, n_pair, z)


def integral_z_infty_dr_r_c_short(c_short, n_pair, z, f1):
    for ij in range(n_pair):
        for k, _ in enumerate(z):
            f1[k, ij] = trapezoid(y=z[:k] * c_short[k, ij], x=z[:k])
    return f1


def integral_z_infty_dr_r2_c_short(c_short, n_pair, z, f2):
    for ij in range(n_pair):
        for k, _ in enumerate(z):
            f2[k, ij] = trapezoid(y=z[:k] * z[:k] * c_short[k, ij], x=z[:k])
    return f2


def calc_hw(tw, n_component, beta_phiw):
    hw = np.zeros((len(tw), n_component))
    for i in range(n_component):
        hw[:, i] = np.exp(tw[:, i] - beta_phiw[:, i]) - 1.0
    return hw


# TODO: Continue from here - check units
# TODO: make initialisation of arrays consistent - not some in Class and others in "calc" functions

def calc_tw(tw_in, beta_phiw, beta_psi_charge, charge_pair, rho, f1, f2, z,
            n_component, n_point, z_index):

    tw = np.zeros((n_point, n_component))
    integral_z_infty = np.zeros((n_point, n_component))
    integral_0_z = np.zeros((n_point, n_component))

    TWO_PI = 2.0 * np.pi
    hw = calc_hw(tw_in, n_component, beta_phiw)

    for i in range(n_component):
        for k in range(n_point):
            integral_0_z[k, i] = trapezoid(y=hw[:k, i], x=z[:k])
            integral_z_infty[k, i] = trapezoid(y=z[k:] * hw[k:, i], x=z[k:])

    for i in range(n_component):
        for k in range(n_point):
            z_minus_t = np.flip(z_index[:k])
            t_minus_z = z_index[k:] - k
            for j in range(i, n_component):
                l = calc_l_index(i, j)
                tw[k, i] = beta_psi_charge[i]
                tw[k, i] += TWO_PI * rho[j] * (z[k] * f1[k, l] - f2[k, l]
                                                + charge_pair[l] * (integral_z_infty[k, j] + z[k] * integral_0_z[k, j])
                                                + trapezoid(y=hw[:k, j] * f1[z_minus_t, l])
                                                + trapezoid(y=hw[k:, j] * f1[t_minus_z, l]))
    return tw


# Documentation: https://scipy.github.io/devdocs/reference/optimize.root-krylov.html

def opt_func(tw_in, beta_phiw, beta_psi_charge, charge_pair, rho, f1, f2, z,
             n_component, n_point, z_index):

    tw = calc_tw(tw_in, beta_phiw, beta_psi_charge, charge_pair, rho, f1, f2, z,
                 n_component, n_point, z_index)

    #dtw_dz = np.zeros((n_point, n_component))
    #for i in range(n_component):
    #    dtw_dz += np.diff(tw[:, i], n=2) / np.diff(z, n=2)
    #print(dtw_dz.shape)

    return tw_in - tw # - 1e-11 * np.sum(np.abs(dtw_dz) * np.repeat(z[1:], 2, axis=0).reshape(n_point-1, n_component))

    #- 0.05 * sum((tw[1:] - tw[:-1])**2 * np.repeat(z[1:], 2, axis=0).reshape(n_point-1, n_component))

#TODO: Work out if regularisation works to keep solution smooth
#TODO: Look at NITSOL and NKSOL parameters to see if any clues?

# https://www.osti.gov/servlets/purl/314885: KINSOL - nonlinear solver based on NKSOL

def solve_model(opt_func, tw_initial, fluid, model, discrete, beta_phiw, beta_psi_charge):
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

    tw_args = (beta_phiw, beta_psi_charge, charge_pair, rho, f1, f2, z,
                n_component, n_point, z_index)
    
    return optim.root(opt_func, tw_initial, args=tw_args,
                      method="krylov", jac=None,
                      tol=tolerance, callback=None, 
                      options={"disp": True, "maxiter": max_iteration})




# callback_func(tw, tw_args) ?



def test_solve_model(opt_func, tw_initial, fluid, model, discrete):
    beta_phiw = model.beta_phiw
    beta_psi_charge = model.beta_psi_charge
    charge_pair = fluid.charge_pair
    rho = fluid.rho
    f1 = model.f1
    f2 = model.f2
    z = model.z
    n_component = fluid.n_component
    n_point = discrete.n_point
    z_index = model.z_index

    tw_args = (beta_phiw, beta_psi_charge, charge_pair, rho, f1, f2, z,
                n_component, n_point, z_index)
    
    tolerance = discrete.tolerance
    max_iteration = discrete.max_iteration

    return tw_args, tolerance, max_iteration