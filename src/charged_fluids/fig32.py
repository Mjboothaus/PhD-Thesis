"""Thesis Fig 3.2 wall–O / wall–H digitised profiles (ψ₀ = 0 mV).

Digitised by hand from ``data/Fig3.2.jpg`` (Case 1, zero surface potential).
Axes: z = 0–15 Å; g_wO ≈ 0–8 (peak ~7); g_wH ≈ 0–9 (peak ~8.5).
Bulk c(r) is held fixed for all CF1 wall calculations (thesis practice).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.interpolate import interp1d

DEFAULT_DATA_DIR = Path(__file__).resolve().parents[2] / "data"


def fig32_csv_paths(data_dir: Path | None = None) -> tuple[Path, Path]:
    root = Path(data_dir) if data_dir is not None else DEFAULT_DATA_DIR
    return root / "Fig3.2_wO_digitised.csv", root / "Fig3.2_wH_digitised.csv"


def load_fig32_profiles(
    data_dir: Path | None = None,
) -> dict[str, np.ndarray]:
    """Return ``z``, ``g_wO``, ``g_wH`` on the digitised grid."""
    p_o, p_h = fig32_csv_paths(data_dir)
    o = np.loadtxt(p_o, delimiter=",", skiprows=1)
    h = np.loadtxt(p_h, delimiter=",", skiprows=1)
    # Prefer oxygen z grid; hydrogen should match
    z = o[:, 0]
    g_o = o[:, 1]
    g_h = interp1d(h[:, 0], h[:, 1], bounds_error=False, fill_value=(h[0, 1], h[-1, 1]))(z)
    return {"z": z, "g_wO": g_o, "g_wH": g_h}


def compare_to_fig32(
    z_live: np.ndarray,
    g_live: np.ndarray,
    *,
    data_dir: Path | None = None,
    z_max: float = 15.0,
) -> dict[str, float | dict]:
    """RMSE / peak metrics of live ``g`` (n_z, 2) vs digitised Fig 3.2.

    Component order must be O, H.
    """
    thesis = load_fig32_profiles(data_dir)
    z_t = thesis["z"]
    mask = z_t <= z_max
    z_t = z_t[mask]
    out: dict[str, float | dict] = {}
    names = ("g_wO", "g_wH")
    for i, name in enumerate(names):
        g_th = thesis[name][mask]
        g_lv = interp1d(
            np.asarray(z_live, float),
            np.asarray(g_live[:, i], float),
            bounds_error=False,
            fill_value=np.nan,
        )(z_t)
        m = np.isfinite(g_th) & np.isfinite(g_lv)
        rmse = float(np.sqrt(np.mean((g_lv[m] - g_th[m]) ** 2))) if m.any() else float("nan")
        ith = int(np.nanargmax(g_th))
        ilv = int(np.nanargmax(g_lv)) if np.isfinite(g_lv).any() else 0
        out[name] = {
            "rmse": rmse,
            "thesis_peak": float(g_th[ith]),
            "thesis_peak_z": float(z_t[ith]),
            "live_peak": float(g_lv[ilv]) if np.isfinite(g_lv[ilv]) else float("nan"),
            "live_peak_z": float(z_t[ilv]) if np.isfinite(g_lv[ilv]) else float("nan"),
            "thesis_g_inf": float(g_th[-1]),
            "live_g_inf": float(g_lv[-1]) if np.isfinite(g_lv[-1]) else float("nan"),
        }
    return out
