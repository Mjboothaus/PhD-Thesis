from dataclasses import dataclass

import numpy as np
import pandas as pd
import scipy.optimize as optim
from scipy import interpolate
from scipy.constants import Avogadro, Boltzmann, elementary_charge, epsilon_0
from scipy.integrate import trapezoid
from streamlit import cache_data
from utils import monitor_memory


@dataclass
class Model:
    z: np.array
    z_index: np.array
    hw: np.array
    c_short: np.array
    f1: np.array
    f2: np.array


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


@monitor_memory
def calc_u(charge, cap_b, alpha, cap_c, cap_d, n_point, n_component, n_pair, epsilon, r):
    # Pre-allocate output array
    u = np.zeros((n_point, n_pair))
    
    # Pre-compute r-dependent terms to avoid repeated calculations
    r_inv = np.zeros_like(r)
    r_inv[1:] = 1.0 / r[1:]  # Avoid division by zero at r[0]
    r_inv_e10 = r_inv * 1e-10
    r_inv_6 = r_inv ** 6
    r_inv_8 = r_inv ** 8
    exp_term = np.exp(-alpha * r)
    
    # Main calculation loop with minimized memory allocations
    for i in range(n_component):
        for j in range(i, n_component):
            l = calc_l_index(i, j)
            
            # Calculate all components without creating temporary arrays
            charge_term = (charge[i] * charge[j]) * r_inv_e10[1:] / epsilon
            repulsive_term = cap_c[l] * r_inv_6[1:] + cap_d[l] * r_inv_8[1:]
            exp_contribution = cap_b[l] * exp_term[1:]
            
            # Combine terms directly into output array
            u[1:, l] = charge_term + repulsive_term + exp_contribution
            u[0, l] = u[1, l]  # Handle r=0 case
    
    return u


# Bulk LJ potential

def calc_u_lj(epsilon_lj, sigma_lj, n_point, n_component, n_pair, r):
    u = np.zeros((n_point, n_pair))
    
    # Handle r[0] separately
    u[0, :] = np.inf  # Set to infinity or another appropriate value
    
    # Process the rest of the points
    for ir in range(1, n_point):
        for i in range(n_component):
            for j in range(i, n_component):
                l = calc_l_index(i, j)
                u[ir, l] = 4.0 * epsilon_lj[l] * ((sigma_lj[l]/r[ir])**12 - (sigma_lj[l]/r[ir])**6)
    
    return u


# Convert mol / dm3 to number / A^3

def calc_rho(concentration):
    return np.array(concentration) * Avogadro / 1.0e27


def calc_kappa(beta, charge, rho, epsilon):
    return np.sqrt(4.0 * np.pi * beta / epsilon * 1e10 *
                   sum(np.multiply(charge**2, rho)))


@cache_data
def calc_phiw(z, n_point, n_component):
    phiw = np.zeros((n_point, n_component))
    capital_a = 16.274e-19  # joules
    wall_d = 2.97  # inverse Angstrom
    for i in range(n_component):
        phiw[:, i] = np.exp(-wall_d * z) * capital_a * (wall_d * z + 2)
    return phiw


@cache_data
def interpolate_cr(r_in, cr_in, n_point, n_pair, z):
    cr = np.zeros((n_point, n_pair))
    for l in range(n_pair):
        f = interpolate.interp1d(r_in, cr_in[:, l])
        r = z
        cr[:, l] = f(r)
        # TODO: Make general - were getting some kinks in c(r) near r=0
        # cr[:10, l] = cr[10, l]
    return cr, r


@cache_data
def load_and_interpolate_cr(cr_path, n_point, n_pair, z):
    cr_df = pd.read_csv(cr_path, header=None, sep='\s+')
    cr_df.set_index(0, inplace=True)
    r = cr_df.index.to_numpy()
    cr = cr_df.to_numpy()
    return interpolate_cr(r, cr, n_point, n_pair, z)


@cache_data
def calc_f1_integrand(c_short, n_pair, z, n_point):
    f1_integrand = np.zeros((n_point, n_pair))
    for ij in range(n_pair):
        f1_integrand[:, ij] = z * c_short[:, ij]
    return f1_integrand


@cache_data
def calc_f2_integrand(c_short, n_pair, z, n_point):
    f2_integrand = np.zeros((n_point, n_pair))
    for ij in range(n_pair):
        f2_integrand[:, ij] = z*z * c_short[:, ij]
    return f2_integrand


@cache_data
@monitor_memory
def integral_z_infty_dr_r_c_short(c_short, n_pair, n_point, z):
    # Pre-allocate output array only
    f1 = np.zeros((n_point, n_pair))
    
    # Process each pair, computing integrals directly without temporary arrays
    for ij in range(n_pair):
        # Use broadcasting instead of temporary array
        zc = z * c_short[:, ij]  # This creates a view, not a copy
        
        # Compute integrals using cumulative trapz for better memory efficiency
        for k in range(n_point):
            if k < n_point - 1:  # Only compute if there are points remaining
                f1[k:, ij] = np.trapz(zc[k:], z[k:])
    
    return f1


@cache_data
@monitor_memory
def integral_z_infty_dr_r2_c_short(c_short, n_pair, n_point, z):
    # Pre-allocate output array only
    f2 = np.zeros((n_point, n_pair))
    
    # Use broadcasting for z*z calculation once
    z2 = z * z  # This creates a view, not a copy
    
    # Process each pair, computing integrals directly
    for ij in range(n_pair):
        # Use broadcasting instead of temporary array
        z2c = z2 * c_short[:, ij]  # This creates a view, not a copy
        
        # Compute integrals using cumulative trapz for better memory efficiency
        for k in range(n_point):
            if k < n_point - 1:  # Only compute if there are points remaining
                f2[k:, ij] = np.trapz(z2c[k:], z[k:])
    
    return f2


def calc_hw(tw, n_component, beta_phiw):
    hw = np.zeros((len(tw), n_component))
    for i in range(n_component):
        hw[:, i] = np.exp(tw[:, i] - beta_phiw[:, i]) - 1.0
    return hw


# TODO: make initialisation of arrays consistent - not some in Class and others in "calc" functions

@monitor_memory
def calc_tw(tw_in, beta_phiw, beta_psi_charge, charge_pair, rho, f1, f2, z,
            n_component, n_point, z_index):
    
    TWO_PI = 2.0 * np.pi
    
    # Pre-allocate output array only
    tw = np.zeros((n_point, n_component))
    
    # Calculate hw using in-place operations
    hw = np.exp(tw_in - beta_phiw) - 1.0
    
    # Calculate integrals using a generator to avoid storing all intermediate results
    def calc_integrals(hw, z, k, i):
        # Calculate integral_0_z using cumulative trapezoid
        if k > 0:
            integral_0_z = np.trapz(hw[:k, i], z[:k])
        else:
            integral_0_z = 0.0
            
        # Calculate integral_z_infty
        if k < n_point:
            z_slice = z[k:]
            integral_z_infty = np.trapz(z_slice * hw[k:, i], z_slice)
        else:
            integral_z_infty = 0.0
            
        return integral_0_z, integral_z_infty
    
    # Main calculation loop with minimized memory allocations
    for i in range(n_component):
        for k in range(n_point):
            integral_0_z_k, integral_z_infty_k = calc_integrals(hw, z, k, i)
            
            for j in range(i, n_component):
                l = calc_l_index(i, j)
                
                # Avoid temporary array allocations by calculating components directly
                tw[k, i] = beta_psi_charge[i]
                
                # Main terms
                main_terms = z[k] * f1[k, l] - f2[k, l]
                
                # Charge pair terms
                charge_terms = charge_pair[l] * (integral_z_infty_k + z[k] * integral_0_z_k)
                
                # Calculate trapezoid terms efficiently
                if k > 0:
                    trap_term1 = np.trapz(hw[:k, j] * f1[np.flip(z_index[:k]), l])
                else:
                    trap_term1 = 0.0
                    
                if k < n_point:
                    trap_term2 = np.trapz(hw[k:, j] * f1[z_index[k:] - k, l])
                else:
                    trap_term2 = 0.0
                
                tw[k, i] += TWO_PI * rho[j] * (main_terms + charge_terms + trap_term1 + trap_term2)
    return tw


# Documentation: https://scipy.github.io/devdocs/reference/optimize.root-krylov.html

def opt_func(tw_in, beta_phiw, beta_psi_charge, charge_pair, rho, f1, f2, z,
             n_component, n_point, z_index):

    tw = calc_tw(tw_in, beta_phiw, beta_psi_charge, charge_pair, rho, f1, f2, z,
                 n_component, n_point, z_index)
    return tw_in - tw


# TODO: Look at NITSOL and NKSOL parameters to see if any clues?

# C.f. https://www.osti.gov/servlets/purl/314885: KINSOL - nonlinear solver based on NKSOL


@monitor_memory
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
    tw_args = beta_phiw, beta_psi_charge, charge_pair, rho, f1, f2, z, n_component, n_point, z_index

    return optim.root(opt_func, tw_initial, args=tw_args, method="krylov", jac=None, 
            callback=None, options={"disp": True, "maxiter": max_iteration, "fatol": tolerance})
