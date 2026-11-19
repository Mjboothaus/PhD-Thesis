"""Shared preparation of fluid + model state for apps and tests."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

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
    calc_phiw,
    calc_rho,
    calc_u,
    calc_u_lj,
    integral_z_infty_dr_r2_c_short,
    integral_z_infty_dr_r_c_short,
    load_and_interpolate_cr,
)
from charged_fluids.numerics import Discretisation, set_num_parameters
from charged_fluids.parameters import Fluid, set_fluid_parameters


def default_data_dir() -> Path:
    return Path(__file__).resolve().parents[2] / "data"


@dataclass
class FluidRun:
    """Fully prepared state for a wall-fluid calculation."""

    fluid: Fluid
    discrete: Discretisation
    model: Model
    potential: dict[str, Any]
    beta_phiw: np.ndarray
    beta_psi_charge: np.ndarray
    beta_u: np.ndarray
    u_r: np.ndarray
    f1_integrand: np.ndarray
    f2_integrand: np.ndarray
    kappa: float
    cr_path: Path
    psi_0_mv: float


def enrich_fluid(fluid: Fluid) -> Fluid:
    """Fill derived thermodynamic fields on a Fluid instance."""
    fluid.beta = calc_beta(fluid.temperature)
    fluid.epsilon = calc_epsilon(fluid.epsilon_r)
    fluid.charge = calc_charge(fluid.valence)
    fluid.charge_pair = calc_charge_pair(
        fluid.beta,
        fluid.charge,
        fluid.epsilon,
        fluid.n_component,
        fluid.n_pair,
    )
    fluid.rho = calc_rho(fluid.concentration)
    return fluid


def _pair_potential(
    fluid: Fluid,
    potential: dict[str, Any],
    n_point: int,
    r: np.ndarray,
) -> np.ndarray:
    pot_type = str(potential.get("type", "")).lower()
    if pot_type == "lj":
        return calc_u_lj(
            np.asarray(potential["epsilon_lj"], dtype=float),
            np.asarray(potential["sigma_lj"], dtype=float),
            n_point,
            fluid.n_component,
            fluid.n_pair,
            r,
        )
    if pot_type == "ionic":
        beta_pauling = calc_beta_pauling(
            fluid.valence,
            np.asarray(potential["n_outer_shell"], dtype=float),
            fluid.n_component,
            fluid.n_pair,
        )
        cap_b = calc_cap_b(
            beta_pauling,
            float(potential["b"]),
            float(potential["alpha"]),
            potential["sigma"],
            fluid.n_component,
            fluid.n_pair,
        )
        return calc_u(
            fluid.charge,
            cap_b,
            float(potential["alpha"]),
            np.asarray(potential["cap_c"], dtype=float),
            np.asarray(potential["cap_d"], dtype=float),
            n_point,
            fluid.n_component,
            fluid.n_pair,
            fluid.epsilon,
            r,
        )
    raise ValueError(f"Unsupported potential type: {pot_type!r} for {fluid.symbol}")


def prepare_fluid_run(
    symbol: str,
    *,
    z_cutoff: float = 25.0,
    n_point: int = 201,
    tolerance: float = 1e-9,
    max_iteration: int = 10,
    psi_0_mv: float = 0.0,
    data_dir: Path | None = None,
    toml_path: Path | None = None,
) -> FluidRun:
    """Build fluid, grid, bulk c(r) integrals, and Model in one place."""
    fluid = set_fluid_parameters(symbol, toml_path=toml_path)
    if fluid is None:
        raise KeyError(f"Unknown fluid symbol: {symbol!r}")
    fluid = enrich_fluid(fluid)
    potential = dict(fluid.potential)

    discrete = set_num_parameters(
        n_point,
        z_cutoff,
        fluid.n_component,
        fluid.n_pair,
        tolerance,
        max_iteration,
    )
    z = np.linspace(0.0, z_cutoff, n_point)
    z_index = np.arange(n_point, dtype=int)
    n_component = fluid.n_component
    n_pair = fluid.n_pair

    wall_zeros = np.zeros((n_point, n_component))

    data_root = Path(data_dir) if data_dir is not None else default_data_dir()
    cr_path = data_root / fluid.cr_filename
    c_short, _ = load_and_interpolate_cr(cr_path, n_point, n_pair, z)
    f1 = integral_z_infty_dr_r_c_short(c_short, n_pair, n_point, z)
    f2 = integral_z_infty_dr_r2_c_short(c_short, n_pair, n_point, z)
    f1_integrand = calc_f1_integrand(c_short, n_pair, z, n_point)
    f2_integrand = calc_f2_integrand(c_short, n_pair, z, n_point)

    u_r = _pair_potential(fluid, potential, n_point, z)
    beta_u = fluid.beta * u_r
    beta_phiw = fluid.beta * calc_phiw(z, n_point, n_component)
    # psi_0 sidebar is in mV; convert to V then multiply by beta and charge
    beta_psi = fluid.beta * float(psi_0_mv) * 1.0e-3
    beta_psi_charge = -beta_psi * fluid.charge

    tw = wall_zeros.copy()
    hw = calc_hw(tw, n_component, beta_phiw)
    model = Model(
        z=z,
        z_index=z_index,
        hw=hw,
        tw=tw,
        c_short=c_short,
        f1=f1,
        f2=f2,
    )

    kappa = calc_kappa(fluid.beta, fluid.charge, fluid.rho, fluid.epsilon)

    return FluidRun(
        fluid=fluid,
        discrete=discrete,
        model=model,
        potential=potential,
        beta_phiw=beta_phiw,
        beta_psi_charge=beta_psi_charge,
        beta_u=beta_u,
        u_r=u_r,
        f1_integrand=f1_integrand,
        f2_integrand=f2_integrand,
        kappa=float(kappa),
        cr_path=cr_path,
        psi_0_mv=float(psi_0_mv),
    )
