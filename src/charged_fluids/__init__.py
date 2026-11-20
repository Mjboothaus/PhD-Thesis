"""Core numerical library for charged fluids near interfaces.

UI frameworks (Streamlit, Marimo) should import from this package rather than
duplicating solver logic.
"""

from charged_fluids.modelling import (
    ConvergenceHistory,
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
    n_pairs,
    opt_func,
    solve_model,
)
from charged_fluids.numerics import Discretisation, set_num_parameters
from charged_fluids.parameters import (
    Fluid,
    clear_fluid_cache,
    fluid_specific_parameters,
    fluids_with_status,
    list_fluid_symbols,
    load_fluid_table,
    set_fluid_parameters,
    wall_ready_symbols,
)
from charged_fluids.setup import FluidRun, enrich_fluid, prepare_fluid_run
from charged_fluids.wall_metrics import contact_and_bulk_g, wall_solution_report

__all__ = [
    "ConvergenceHistory",
    "Discretisation",
    "Fluid",
    "FluidRun",
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
    "clear_fluid_cache",
    "contact_and_bulk_g",
    "enrich_fluid",
    "fluid_specific_parameters",
    "fluids_with_status",
    "integral_z_infty_dr_r2_c_short",
    "integral_z_infty_dr_r_c_short",
    "interpolate_cr",
    "list_fluid_symbols",
    "load_and_interpolate_cr",
    "load_fluid_table",
    "n_pairs",
    "opt_func",
    "prepare_fluid_run",
    "set_fluid_parameters",
    "set_num_parameters",
    "solve_model",
    "wall_ready_symbols",
    "wall_solution_report",
]

__version__ = "0.2.0"
