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
from charged_fluids.parameters import load_fluid_table, set_fluid_parameters
from charged_fluids.setup import enrich_fluid

_K_J_PER_KJ_MOL = 1000.0 / Avogadro


def c_short_tail_max(result: BulkResult, *, fraction: float = 0.05) -> float:
    """Max |c_short| on the outer ``fraction`` of the radial grid."""
    n = result.r.size
    n_tail = max(5, int(n * fraction))
    return float(np.max(np.abs(result.c_short[:, -n_tail:])))


def wall_usable_c_short(result: BulkResult, *, tail_tol: float = 1e-2) -> bool:
    """Whether short-range c is small enough at large r for wall-fluid input."""
    return c_short_tail_max(result) < tail_tol


def bulk_potential_type(symbol: str, toml_path: Path | None = None) -> str | None:
    """Return potential type for ``symbol``, or None if unknown / no potential."""
    fluid = set_fluid_parameters(symbol, toml_path=toml_path)
    if fluid is None:
        return None
    pot = fluid.potential or {}
    t = str(pot.get("type", "")).lower().strip()
    return t or None


def list_bulk_symbols(toml_path: Path | None = None) -> list[str]:
    """Symbols that ``generate_bulk_cr`` can handle (lj, ionic binary, cf1)."""
    table = load_fluid_table(toml_path)
    out: list[str] = []
    for key, raw in table.items():
        if not isinstance(raw, dict) or "component" not in raw:
            continue
        pot = dict(raw.get("potential", {}) or {})
        pot_type = str(pot.get("type", "")).lower()
        n_comp = len(raw.get("component", []))
        if pot_type == "lj" or (pot_type == "ionic" and n_comp == 2) or pot_type == "cf1":
            out.append(str(raw.get("symbol", key)))
    return sorted(set(out))


def generate_all_bulk_cr(
    *,
    bulk_dir: Path | None = None,
    toml_path: Path | None = None,
    n_points: int | None = None,
    delta_r: float | None = None,
    closure: str | None = None,
    write: bool = True,
) -> list[tuple[str, BulkResult, Path | None, Path | None]]:
    """Generate bulk artefacts for every supported fluid symbol."""
    rows: list[tuple[str, BulkResult, Path | None, Path | None]] = []
    for symbol in list_bulk_symbols(toml_path=toml_path):
        result, cr_path, meta_path = generate_bulk_cr(
            symbol,
            bulk_dir=bulk_dir,
            toml_path=toml_path,
            n_points=n_points,
            delta_r=delta_r,
            closure=closure,
            write=write,
        )
        rows.append((symbol, result, cr_path, meta_path))
    return rows


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

    table = load_fluid_table(toml_path)
    bulk = dict(table.get(symbol, {}).get("bulk", {}) or {})

    n_pts = int(n_points if n_points is not None else bulk.get("n_points", 1024))
    dr = float(delta_r if delta_r is not None else bulk.get("delta_r", 0.05))
    clos = cast(
        Literal["hnc", "py"],
        str(closure if closure is not None else bulk.get("closure", "hnc")),
    )

    if pot_type == "lj":
        sigma = float(bulk.get("sigma_A", float(np.asarray(pot["sigma_lj"]).ravel()[0])))
        if "epsilon_kJ_mol" in bulk:
            epsilon = float(bulk["epsilon_kJ_mol"]) * _K_J_PER_KJ_MOL
        else:
            eps_star = float(np.asarray(pot["epsilon_lj"]).ravel()[0])
            epsilon = eps_star * Boltzmann * fluid.temperature
        density = float(np.sum(fluid.rho))
        result = solve_bulk_lj(
            temperature=fluid.temperature,
            density=density,
            epsilon=epsilon,
            sigma=sigma,
            n_points=n_pts,
            delta_r=dr,
            closure=clos,
            mix=0.75,
            max_iter=2000,
            tol=1e-9,
        )
        if fluid.n_component > 1:
            n_pair = fluid.n_pair
            c1 = result.c_short[0]
            result.c_full = np.vstack([c1] * n_pair)
            result.c_short = result.c_full.copy()
            result.g = np.vstack([result.g[0]] * n_pair)
            result.gamma = np.vstack([result.gamma[0]] * n_pair)
            result.rho = fluid.rho.copy()
            result.meta["n_component"] = fluid.n_component
            result.meta["note"] = "Equal multi-comp LJ: pair c(r) from 1-comp OZ at total density"
    elif pot_type == "ionic" and fluid.n_component == 2:
        sig = np.asarray(pot["sigma"], dtype=float)
        eps_kj = float(bulk.get("epsilon_kJ_mol", 0.5))
        epsilon = eps_kj * _K_J_PER_KJ_MOL
        eps_lj = (epsilon, epsilon, epsilon)
        dens = (float(fluid.rho[0]), float(fluid.rho[1]))
        result = solve_bulk_binary_coulomb_lj(
            temperature=fluid.temperature,
            densities=dens,
            charge_pair=np.asarray(fluid.charge_pair, dtype=float),
            sigmas=(float(sig[0]), float(sig[1])),
            epsilons_lj=eps_lj,
            n_points=n_pts,
            delta_r=dr,
            closure=clos,
            mix=float(bulk.get("mix", 0.12)),
            max_iter=int(bulk.get("max_iter", 2000)),
            tol=float(bulk.get("tol", 1e-6)),
            ng_alpha=float(bulk.get("ng_alpha", 1.2)),
            anneal_steps=int(bulk.get("anneal_steps", 1)),
            use_anderson=bool(bulk.get("use_anderson", True)),
            t_hot=float(bulk["t_hot"]) if "t_hot" in bulk else None,
        )
        result.meta["note"] = (
            "Ng bulk OZ (+ optional T anneal); c_short uses wall-fluid charge_pair Coulomb scale"
        )
    elif pot_type == "cf1":
        from charged_fluids.bulk.oz_solver import solve_bulk_cf1

        source = str(bulk.get("source", "oz")).lower()
        if source in {"digitised_fig3_12", "digitised", "fig3_12"}:
            # Explicit opt-in to thesis digitised bulk (legacy path)
            from charged_fluids.bulk.artefacts import read_meta
            from charged_fluids.bulk.build_digitised_h2o import build_h2o_cr
            from charged_fluids.bulk.oz_solver import BulkResult

            out_dir = bulk_dir or default_bulk_dir()
            cr_path, meta_path = build_h2o_cr(
                bulk_dir=out_dir,
                n_points=n_pts,
                r_max=float(bulk.get("r_max", 9.0)),
            )
            data = np.loadtxt(cr_path)
            r = data[:, 0]
            c_short = data[:, 1:].T
            meta = read_meta(meta_path) if meta_path.is_file() else {}
            result = BulkResult(
                r=r,
                c_full=c_short.copy(),
                c_short=c_short,
                g=np.ones_like(c_short),
                gamma=np.zeros_like(c_short),
                n_iter=0,
                converged=True,
                residual=0.0,
                beta=fluid.beta,
                rho=fluid.rho.copy(),
                meta={
                    "fluid": "cf1_digitised",
                    "note": "Digitised thesis Fig 3.12 bulk c(r); not a live OZ solve",
                    **{k: v for k, v in meta.items() if k not in {"created_at"}},
                },
            )
            result.meta["fluid_symbol"] = fluid.symbol
            result.meta["cr_filename_legacy"] = fluid.cr_filename
            result.meta["bulk_params"] = bulk
            if write:
                return result, cr_path, meta_path
            return result, None, None

        dens = (float(fluid.rho[0]), float(fluid.rho[1]))
        result = solve_bulk_cf1(
            temperature=fluid.temperature,
            densities=dens,
            charge_pair=np.asarray(fluid.charge_pair, dtype=float),
            n_points=n_pts,
            delta_r=dr,
            closure=clos,
            mix=float(bulk.get("mix", 0.03)),
            max_iter=int(bulk.get("max_iter", 4000)),
            tol=float(bulk.get("tol", 5.0e-5)),
            ng_alpha=float(bulk.get("ng_alpha", 1.5)),
            anneal_steps=int(bulk.get("anneal_steps", 24)),
            use_anderson=bool(bulk.get("use_anderson", True)),
            t_hot=float(bulk["t_hot"]) if "t_hot" in bulk else None,
        )
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
