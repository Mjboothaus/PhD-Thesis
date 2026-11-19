"""Build data/bulk/h2o.cr.dat from thesis digitised Fig 3.12 pair curves.

Pair order: OO, OH, HH (components O, H).
"""

from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
from scipy import interpolate

from charged_fluids.bulk.artefacts import default_bulk_dir


def _load_pair(csv_path: Path) -> tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(csv_path, delimiter=",")
    r = data[:, 0]
    y = data[:, 1]
    order = np.argsort(r)
    r, y = r[order], y[order]
    # Drop non-positive r (digitising noise)
    mask = r > 1.0e-4
    return r[mask], y[mask]


def build_h2o_cr(
    *,
    data_dir: Path | None = None,
    bulk_dir: Path | None = None,
    n_points: int = 512,
    r_max: float = 9.0,
) -> tuple[Path, Path]:
    root = Path(data_dir) if data_dir is not None else Path(__file__).resolve().parents[3] / "data"
    bulk = Path(bulk_dir) if bulk_dir is not None else default_bulk_dir()
    bulk.mkdir(parents=True, exist_ok=True)

    pairs = []
    for name in ("OO", "OH", "HH"):
        r, y = _load_pair(root / f"Fig3.12_{name}_digitised.csv")
        pairs.append((r, y))

    r_grid = np.linspace(r_max / n_points, r_max, n_points)
    cols = [r_grid]
    for r, y in pairs:
        f = interpolate.interp1d(r, y, bounds_error=False, fill_value=(y[0], 0.0))
        cols.append(f(r_grid))
    table = np.column_stack(cols)

    cr_path = bulk / "h2o.cr.dat"
    meta_path = bulk / "h2o.meta.json"
    np.savetxt(cr_path, table, fmt="%12.6f")
    meta = {
        "created_at": datetime.now(UTC).isoformat(),
        "stem": "h2o",
        "use_short": True,
        "source": "digitised_fig3_12",
        "pair_order": ["OO", "OH", "HH"],
        "n_points": n_points,
        "r_max": r_max,
        "note": "Thesis Fig 3.12 digitised bulk c(r) for CF1 water pairs",
    }
    meta_path.write_text(json.dumps(meta, indent=2) + "\n")
    return cr_path, meta_path


if __name__ == "__main__":
    cr, meta = build_h2o_cr()
    print(f"wrote {cr}")
    print(f"meta  {meta}")
