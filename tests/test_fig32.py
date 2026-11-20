"""Digitised thesis Fig 3.2 wall profiles."""

from __future__ import annotations

import numpy as np
import pytest

from charged_fluids import compare_to_fig32, load_fig32_profiles


def test_load_fig32_peaks_and_range():
    p = load_fig32_profiles()
    z = p["z"]
    assert z[0] == pytest.approx(0.0)
    assert z[-1] == pytest.approx(15.0)
    # Thesis contact peaks (hand-digitised; allow modest tolerance)
    assert float(np.max(p["g_wO"])) == pytest.approx(7.0, abs=0.35)
    assert float(np.max(p["g_wH"])) == pytest.approx(8.5, abs=0.5)
    assert float(z[np.argmax(p["g_wO"])]) == pytest.approx(1.05, abs=0.25)
    assert float(z[np.argmax(p["g_wH"])]) == pytest.approx(1.0, abs=0.25)
    # Far field order-1 (not bulk 1.0 exactly — figure asymptote ~1.6–1.8)
    assert 1.2 <= float(p["g_wO"][-1]) <= 2.2
    assert 1.2 <= float(p["g_wH"][-1]) <= 2.2


def test_compare_to_fig32_self_rmse_near_zero():
    p = load_fig32_profiles()
    g = np.column_stack([p["g_wO"], p["g_wH"]])
    m = compare_to_fig32(p["z"], g)
    assert m["g_wO"]["rmse"] == pytest.approx(0.0, abs=1e-6)
    assert m["g_wH"]["rmse"] == pytest.approx(0.0, abs=1e-6)
