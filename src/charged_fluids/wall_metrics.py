"""Wall-solution diagnostics (UI-free).

Summarise Newton–Krylov residual history and wall–fluid ``g_wi(z)`` contact /
far-field values for tests, Streamlit, and Marimo.
"""

from __future__ import annotations

from typing import Any

import numpy as np


def contact_and_bulk_g(
    z: np.ndarray,
    g: np.ndarray,
    *,
    near_cut: float = 5.0,
) -> list[dict[str, float]]:
    """Per-component contact (first g>0), near-wall peak, and g(∞)."""
    z = np.asarray(z, dtype=float)
    g = np.asarray(g, dtype=float)
    if g.ndim == 1:
        g = g.reshape(-1, 1)
    rows: list[dict[str, float]] = []
    for i in range(g.shape[1]):
        gi = g[:, i]
        finite = np.isfinite(gi) & (gi > 1.0e-12)
        if finite.any():
            j = int(np.argmax(finite))
            g_c, z_c = float(gi[j]), float(z[j])
        else:
            g_c, z_c = float("nan"), float("nan")
        near = z <= min(near_cut, float(z[-1]) if z.size else near_cut)
        if near.any() and np.isfinite(gi[near]).any():
            jpk = int(np.nanargmax(np.where(near, gi, np.nan)))
            g_pk, z_pk = float(gi[jpk]), float(z[jpk])
        else:
            g_pk, z_pk = float("nan"), float("nan")
        rows.append(
            {
                "g_contact": g_c,
                "z_contact": z_c,
                "g_peak": g_pk,
                "z_peak": z_pk,
                "g_inf": float(gi[-1]) if gi.size else float("nan"),
            }
        )
    return rows


def wall_solution_report(
    solution: Any,
    *,
    z: np.ndarray | None = None,
    g: np.ndarray | None = None,
    components: list[str] | None = None,
) -> dict[str, Any]:
    """Compact report from ``solve_model`` result (+ optional g profiles)."""
    summary = dict(getattr(solution, "convergence_summary", {}) or {})
    hist = getattr(solution, "convergence_history", None)
    report: dict[str, Any] = {
        "success": bool(getattr(solution, "success", False)),
        "nit": int(getattr(solution, "nit", 0) or 0),
        "nfev": int(getattr(solution, "nfev", 0) or 0),
        "message": str(getattr(solution, "message", "") or ""),
        "convergence_summary": summary,
        "n_history": len(hist.iterations) if hist is not None else 0,
    }
    if z is not None and g is not None:
        rows = contact_and_bulk_g(z, g)
        if components is not None:
            for i, row in enumerate(rows):
                if i < len(components):
                    row["species"] = components[i]
        report["g_profiles"] = rows
        g_arr = np.asarray(g, dtype=float)
        if g_arr.ndim == 1:
            g_arr = g_arr.reshape(-1, 1)
        report["g_inf_max_abs_err"] = float(np.max(np.abs(g_arr[-1, :] - 1.0)))
    return report
