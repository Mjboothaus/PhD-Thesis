"""Core numerical library for charged fluids near interfaces.

UI frameworks (Streamlit, Marimo) should import from this package rather than
duplicating solver logic.
"""

from charged_fluids.modelling import (
    Model,
    calc_beta,
    calc_beta_pauling,
    calc_cap_b,
    calc_charge,
    calc_charge_pair,
    calc_epsilon,
    calc_f1_integrand,
    calc_f2_integrand,
    calc_hw,
    calc_kappa,
    calc_l_index,
    calc_phiw,
    calc_rho,
    calc_tw,
    calc_u,
    calc_u_lj,
    integral_z_infty_dr_r2_c_short,
    integral_z_infty_dr_r_c_short,
    interpolate_cr,
    load_and_interpolate_cr,
    opt_func,
    solve_model,
)
from charged_fluids.numerics import Discretisation, set_num_parameters
from charged_fluids.parameters import (
    Fluid,
    fluid_specific_parameters,
    set_fluid_parameters,
)

__all__ = [
    "Discretisation",
    "Fluid",
    "Model",
    "calc_beta",
    "calc_beta_pauling",
    "calc_cap_b",
    "calc_charge",
    "calc_charge_pair",
    "calc_epsilon",
    "calc_f1_integrand",
    "calc_f2_integrand",
    "calc_hw",
    "calc_kappa",
    "calc_l_index",
    "calc_phiw",
    "calc_rho",
    "calc_tw",
    "calc_u",
    "calc_u_lj",
    "fluid_specific_parameters",
    "integral_z_infty_dr_r2_c_short",
    "integral_z_infty_dr_r_c_short",
    "interpolate_cr",
    "load_and_interpolate_cr",
    "opt_func",
    "set_fluid_parameters",
    "set_num_parameters",
    "solve_model",
]

__version__ = "0.2.0"
