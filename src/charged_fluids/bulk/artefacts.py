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
