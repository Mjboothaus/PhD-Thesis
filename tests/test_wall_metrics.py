"""Wall diagnostics helpers."""

from __future__ import annotations

import numpy as np
import pytest

from charged_fluids import contact_and_bulk_g, wall_solution_report


def test_contact_and_bulk_g_peak_and_inf():
    z = np.linspace(0.0, 10.0, 21)
    g = np.zeros((21, 2))
    g[2:, 0] = 1.1
    g[3, 0] = 4.0
    g[2:, 1] = 0.9
    g[5, 1] = 2.5
    g[-1, :] = 1.0
    rows = contact_and_bulk_g(z, g)
    assert rows[0]["g_peak"] == pytest.approx(4.0)
    assert rows[1]["g_peak"] == pytest.approx(2.5)
    assert rows[0]["g_inf"] == pytest.approx(1.0)


def test_wall_solution_report_with_fake_solution():
    class _Sol:
        success = True
        nit = 3
        nfev = 10
        message = "ok"
        convergence_summary = {"final_residual": 1e-6}
        convergence_history = type("H", (), {"iterations": [1, 2, 3]})()

    z = np.linspace(0.0, 5.0, 11)
    g = np.ones((11, 1))
    g[0, 0] = 0.0
    g[2, 0] = 3.0
    rep = wall_solution_report(_Sol(), z=z, g=g, components=["LJ"])
    assert rep["success"] is True
    assert rep["n_history"] == 3
    assert rep["g_profiles"][0]["species"] == "LJ"
    assert rep["g_inf_max_abs_err"] == pytest.approx(0.0)
