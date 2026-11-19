"""Persist bulk c(r) artefacts with JSON metadata sidecars."""

from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np

from charged_fluids.bulk.oz_solver import BulkResult


def default_bulk_dir() -> Path:
    return Path(__file__).resolve().parents[3] / "data" / "bulk"


def write_cr_artefact(
    result: BulkResult,
    *,
    stem: str,
    bulk_dir: Path | None = None,
    use_short: bool = True,
) -> tuple[Path, Path]:
    """Write ``stem.cr.dat`` and ``stem.meta.json`` under bulk_dir.

    Format matches wall-fluid loader: whitespace columns ``r c0 [c1 c2 ...]``.
    """
    out_dir = Path(bulk_dir) if bulk_dir is not None else default_bulk_dir()
    out_dir.mkdir(parents=True, exist_ok=True)
    cr_path = out_dir / f"{stem}.cr.dat"
    meta_path = out_dir / f"{stem}.meta.json"

    c = result.c_short if use_short else result.c_full
    # shape (n_pair, n_r) → columns
    data = np.column_stack([result.r, c.T])
    np.savetxt(cr_path, data, fmt="%12.6f")

    meta: dict[str, Any] = {
        "created_at": datetime.now(UTC).isoformat(),
        "stem": stem,
        "use_short": use_short,
        "converged": result.converged,
        "n_iter": result.n_iter,
        "residual": result.residual,
        "beta": result.beta,
        "rho": result.rho.tolist(),
        "n_pair": int(c.shape[0]),
        "n_points": int(result.r.size),
        "r_max": float(result.r[-1]),
        "solver_meta": result.meta,
        "tail_abs_max_c_short": float(
            np.max(np.abs(result.c_short[:, -max(5, result.r.size // 20) :]))
        ),
    }
    meta_path.write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    return cr_path, meta_path


def read_meta(meta_path: Path) -> dict[str, Any]:
    return json.loads(Path(meta_path).read_text(encoding="utf-8"))


def resolve_cr_path(
    fluid_symbol: str,
    *,
    data_dir: Path,
    bulk_dir: Path | None = None,
    legacy_name: str | None = None,
) -> Path:
    """Prefer ``data/bulk/{symbol}.cr.dat`` then legacy ``data/{legacy_name}``."""
    bdir = Path(bulk_dir) if bulk_dir is not None else default_bulk_dir()
    preferred = bdir / f"{fluid_symbol}.cr.dat"
    if preferred.is_file():
        return preferred
    if legacy_name:
        legacy = Path(data_dir) / legacy_name
        if legacy.is_file():
            return legacy
    raise FileNotFoundError(
        f"No bulk c(r) for {fluid_symbol!r}: tried {preferred}"
        + (f" and {Path(data_dir) / legacy_name}" if legacy_name else "")
    )


def validate_cr_artefact(
    cr_path: Path,
    *,
    expected_n_pair: int | None = None,
    tail_tol: float = 1e-2,
) -> dict[str, Any]:
    """Load a ``.cr.dat`` file and return basic validation metrics.

    Raises ``ValueError`` if the file is empty, malformed, or the tail is too large
    when ``tail_tol`` is finite.
    """
    path = Path(cr_path)
    data = np.loadtxt(path)
    if data.ndim != 2 or data.shape[0] < 2 or data.shape[1] < 2:
        raise ValueError(
            f"Malformed bulk c(r) artefact: {path} shape={getattr(data, 'shape', None)}"
        )
    r = data[:, 0]
    c = data[:, 1:]
    n_pair = int(c.shape[1])
    if expected_n_pair is not None and n_pair != expected_n_pair:
        raise ValueError(f"{path}: expected {expected_n_pair} pair columns, got {n_pair}")
    if not np.all(np.diff(r) > 0):
        raise ValueError(f"{path}: radial grid is not strictly increasing")
    n_tail = max(5, c.shape[0] // 20)
    tail = float(np.max(np.abs(c[-n_tail:, :])))
    if tail_tol is not None and tail >= tail_tol:
        raise ValueError(f"{path}: |c_short| tail max={tail:.3e} exceeds tol={tail_tol}")
    return {
        "path": str(path),
        "n_points": int(r.size),
        "n_pair": n_pair,
        "r_min": float(r[0]),
        "r_max": float(r[-1]),
        "tail_abs_max": tail,
    }
