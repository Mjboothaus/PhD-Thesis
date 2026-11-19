"""High-level bulk c(r) generation from fluid symbols / TOML."""

from __future__ import annotations

from pathlib import Path
from typing import Literal, cast

import numpy as np
from scipy.constants import Avogadro, Boltzmann

from charged_fluids.bulk.artefacts import default_bulk_dir, write_cr_artefact
from charged_fluids.bulk.oz_solver import (
    BulkResult,
    solve_bulk_binary_coulomb_lj,
    solve_bulk_lj,
)
from charged_fluids.parameters import set_fluid_parameters
from charged_fluids.setup import enrich_fluid

_K_J_PER_KJ_MOL = 1000.0 / Avogadro


def generate_bulk_cr(
    symbol: str,
    *,
    bulk_dir: Path | None = None,
    toml_path: Path | None = None,
    n_points: int | None = None,
    delta_r: float | None = None,
    closure: str | None = None,
    write: bool = True,
) -> tuple[BulkResult, Path | None, Path | None]:
    """Solve bulk OZ for ``symbol`` and optionally write ``data/bulk`` artefacts."""
    fluid = set_fluid_parameters(symbol, toml_path=toml_path)
    if fluid is None:
        raise KeyError(f"Unknown fluid: {symbol!r}")
    fluid = enrich_fluid(fluid)
    pot = fluid.potential
    pot_type = str(pot.get("type", "")).lower()
    from charged_fluids.parameters import load_fluid_table

    table = load_fluid_table(toml_path)
    bulk = dict(table.get(symbol, {}).get("bulk", {}) or {})

    n_pts = int(n_points if n_points is not None else bulk.get("n_points", 1024))
    dr = float(delta_r if delta_r is not None else bulk.get("delta_r", 0.05))
    clos = str(closure if closure is not None else bulk.get("closure", "hnc"))

    if pot_type == "lj" and fluid.n_component == 1:
        sigma = float(bulk.get("sigma_A", float(np.asarray(pot["sigma_lj"]).ravel()[0])))
        if "epsilon_kJ_mol" in bulk:
            epsilon = float(bulk["epsilon_kJ_mol"]) * _K_J_PER_KJ_MOL
        else:
            # Treat epsilon_lj as reduced eps* => epsilon = eps* kT
            eps_star = float(np.asarray(pot["epsilon_lj"]).ravel()[0])
            epsilon = eps_star * Boltzmann * fluid.temperature
        density = float(fluid.rho[0])
        result = solve_bulk_lj(
            temperature=fluid.temperature,
            density=density,
            epsilon=epsilon,
            sigma=sigma,
            n_points=n_pts,
            delta_r=dr,
            closure=cast(Literal["hnc", "py"], clos),
            mix=0.75,
            max_iter=2000,
            tol=1e-9,
        )
    elif pot_type == "ionic" and fluid.n_component == 2:
        sig = np.asarray(pot["sigma"], dtype=float)
        eps_kj = float(bulk.get("epsilon_kJ_mol", 0.5))
        epsilon = eps_kj * _K_J_PER_KJ_MOL
        eps_lj = (epsilon, epsilon, epsilon)
        dens = (float(fluid.rho[0]), float(fluid.rho[1]))
        charges = (float(fluid.charge[0]), float(fluid.charge[1]))
        result = solve_bulk_binary_coulomb_lj(
            temperature=fluid.temperature,
            densities=dens,
            charges=charges,
            epsilon=float(fluid.epsilon),
            sigmas=(float(sig[0]), float(sig[1])),
            epsilons_lj=eps_lj,
            n_points=n_pts,
            delta_r=dr,
            closure=cast(Literal["hnc", "py"], clos),
            mix=0.1,
            max_iter=2000,
            tol=1e-6,
        )
        result.meta["note"] = "c_short = c_full + beta * u_coulomb (thesis long-range split)"
    else:
        raise ValueError(
            f"Bulk generation not configured for symbol={symbol!r} type={pot_type!r} "
            f"n_component={fluid.n_component}"
        )

    result.meta["fluid_symbol"] = fluid.symbol
    result.meta["cr_filename_legacy"] = fluid.cr_filename
    result.meta["bulk_params"] = bulk

    cr_path = meta_path = None
    if write:
        cr_path, meta_path = write_cr_artefact(
            result,
            stem=fluid.symbol,
            bulk_dir=bulk_dir or default_bulk_dir(),
            use_short=True,
        )
    return result, cr_path, meta_path
