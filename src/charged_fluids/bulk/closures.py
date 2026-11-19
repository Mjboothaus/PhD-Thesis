"""Named bulk OZ closures (HNC, PY) and a small registry.

Additional closures (MSA, Rogers-Young, ...) can register here without touching
the Picard loop. Extend the registry as new closures are validated.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Literal

import numpy as np

ClosureName = Literal["hnc", "py"]

ClosureFn = Callable[[np.ndarray, np.ndarray], tuple[np.ndarray, np.ndarray]]


def hnc_closure(beta_u: np.ndarray, gamma: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Hypernetted-chain: c = exp(-betau + gamma) - 1 - gamma, g = exp(-betau + gamma)."""
    bu = np.where(np.isfinite(beta_u), beta_u, 1.0e6)
    t = np.clip(-bu + gamma, -100.0, 50.0)
    e = np.exp(t)
    return e - 1.0 - gamma, e


def py_closure(beta_u: np.ndarray, gamma: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Percus-Yevick: c = (e^{-betau} - 1)(1 + gamma), g = e^{-betau}(1 + gamma)."""
    bu = np.where(np.isfinite(beta_u), beta_u, 1.0e6)
    t = np.clip(-bu, -100.0, 50.0)
    f = np.exp(t) - 1.0
    c = f * (1.0 + gamma)
    g = (1.0 + f) * (1.0 + gamma)
    return c, g


CLOSURE_REGISTRY: dict[str, ClosureFn] = {
    "hnc": hnc_closure,
    "py": py_closure,
}


def get_closure(name: str) -> ClosureFn:
    """Return a closure callable by lower-case name."""
    key = name.lower().strip()
    try:
        return CLOSURE_REGISTRY[key]
    except KeyError as exc:
        known = ", ".join(sorted(CLOSURE_REGISTRY))
        raise ValueError(f"Unknown closure {name!r}; choose one of: {known}") from exc


def apply_closure(
    closure: str | ClosureFn,
    beta_u: np.ndarray,
    gamma: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Apply a named or callable closure; return ``(c, g)``."""
    fn: ClosureFn = get_closure(closure) if isinstance(closure, str) else closure
    return fn(beta_u, gamma)


def list_closures() -> list[str]:
    return sorted(CLOSURE_REGISTRY)
