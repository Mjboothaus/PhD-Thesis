"""Fluid parameter definitions loaded from TOML (single source of truth)."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover
    import tomli as tomllib  # type: ignore


def _default_toml_path() -> Path:
    # src/charged_fluids/parameters.py -> repo root data/
    return Path(__file__).resolve().parents[2] / "data" / "fluid_parameters.toml"


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
    charge_pair: np.ndarray
    rho: np.ndarray
    beta: float
    epsilon: float
    cr_filename: str = ""
    status: str = "experimental"
    potential: dict[str, Any] = field(default_factory=dict)


@lru_cache(maxsize=4)
def _load_toml(path_str: str) -> dict[str, Any]:
    path = Path(path_str)
    with path.open("rb") as fh:
        return tomllib.load(fh)


def load_fluid_table(toml_path: Path | None = None) -> dict[str, Any]:
    """Return the raw TOML table (section key -> mapping)."""
    path = Path(toml_path) if toml_path is not None else _default_toml_path()
    if not path.is_file():
        raise FileNotFoundError(f"Fluid parameter TOML not found: {path}")
    return _load_toml(str(path.resolve()))


def list_fluid_symbols(toml_path: Path | None = None) -> list[str]:
    return sorted(load_fluid_table(toml_path).keys())


def clear_fluid_cache() -> None:
    """Drop cached TOML tables (call after editing the parameter file in-process)."""
    _load_toml.cache_clear()


def fluids_with_status(*statuses: str, toml_path: Path | None = None) -> list[str]:
    """Return fluid symbols whose TOML ``status`` is in ``statuses``."""
    wanted = {s.lower() for s in statuses}
    table = load_fluid_table(toml_path)
    out: list[str] = []
    for key, raw in table.items():
        if not isinstance(raw, dict):
            continue
        if "component" not in raw:
            continue  # skip nested bulk-only mistakes
        status = str(raw.get("status", "")).lower()
        if status in wanted:
            out.append(str(raw.get("symbol", key)))
    return sorted(set(out))


def wall_ready_symbols(toml_path: Path | None = None) -> list[str]:
    """Fluids that can attempt a wall solve (ready or experimental with bulk data)."""
    return fluids_with_status("ready", "experimental", toml_path=toml_path)


def set_fluid_parameters(symbol: str, toml_path: Path | None = None) -> Fluid | None:
    """Return a :class:`Fluid` for ``symbol``, or ``None`` if unknown."""
    symbol = symbol.lower()
    table = load_fluid_table(toml_path)
    if symbol not in table:
        return None

    raw = deepcopy(table[symbol])
    potential = dict(raw.pop("potential", {}) or {})
    component = list(raw["component"])
    n_component = len(component)
    n_pair = int((n_component + 1) * n_component / 2)
    valence = np.asarray(raw["valence"], dtype=float)
    concentration = np.asarray(raw["concentration"], dtype=float)

    return Fluid(
        name=str(raw["name"]),
        symbol=str(raw.get("symbol", symbol)),
        component=component,
        valence=valence,
        charge=np.zeros(n_component, dtype=float),
        temperature=float(raw["temperature"]),
        concentration=concentration,
        epsilon_r=float(raw["epsilon_r"]),
        n_component=n_component,
        n_pair=n_pair,
        index=int(raw.get("index", 0)),
        charge_pair=np.zeros(n_pair, dtype=float),
        rho=np.zeros(n_component, dtype=float),
        beta=0.0,
        epsilon=0.0,
        cr_filename=str(raw.get("cr_filename", "")),
        status=str(raw.get("status", "experimental")),
        potential=_arrays_in_potential(potential),
    )


def _arrays_in_potential(potential: dict[str, Any]) -> dict[str, Any]:
    out = dict(potential)
    for key, value in list(out.items()):
        if isinstance(value, list):
            out[key] = np.asarray(value, dtype=float)
    return out


def fluid_specific_parameters(symbol: str, toml_path: Path | None = None) -> dict[str, Any]:
    """Potential-model extras keyed by fluid symbol (legacy shape).

    Returns ``{symbol: potential_dict}`` for call-site compatibility with the
    Streamlit pages. Arrays are materialised as :class:`numpy.ndarray`.
    """
    fluid = set_fluid_parameters(symbol, toml_path=toml_path)
    if fluid is None:
        raise KeyError(f"No fluid-specific parameters defined for symbol={symbol!r}")
    if not fluid.potential:
        raise KeyError(f"No fluid-specific parameters defined for symbol={symbol!r}")

    pot = deepcopy(fluid.potential)
    # Promote list fields to arrays where numerical code expects them
    for key, value in list(pot.items()):
        if isinstance(value, list):
            pot[key] = np.asarray(value, dtype=float)
    return {fluid.symbol: pot}


# Back-compat alias used by older call sites / docs
def fluid_parameters(toml_path: Path | None = None) -> dict[str, Any]:
    return load_fluid_table(toml_path)
