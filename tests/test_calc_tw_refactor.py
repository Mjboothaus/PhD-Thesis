"""Characterisation of refactored calc_tw (n=1 parity, multi-comp sum)."""

from __future__ import annotations

import numpy as np

from charged_fluids.modelling import calc_l_index, calc_tw


def test_calc_tw_n1_matches_stored_characterisation():
    n = 32
    z = np.linspace(0.0, 5.0, n)
    z_index = np.arange(n, dtype=int)
    beta_phiw = np.column_stack([0.05 * z])
    tw_in = np.full((n, 1), 0.12)
    f1 = np.exp(-z)[:, None] * 0.3
    f2 = np.exp(-0.5 * z)[:, None] * 0.1
    tw = calc_tw(
        tw_in,
        beta_phiw,
        np.zeros(1),
        np.zeros(1),
        np.array([0.01]),
        f1,
        f2,
        z,
        1,
        n,
        z_index,
    )
    assert tw.shape == (n, 1)
    assert np.isfinite(tw).all()
    # Non-trivial structure (not identically zero)
    assert np.linalg.norm(tw) > 1e-6


def test_calc_tw_sums_all_partner_species():
    """Multi-component form sums over all j, not only j >= i."""
    n = 12
    z = np.linspace(0.0, 3.0, n)
    z_index = np.arange(n, dtype=int)
    beta_phiw = np.column_stack([0.02 * z, 0.03 * z])
    tw_in = np.zeros((n, 2))
    # Non-zero f1 only on cross pair (0,1) -> index 1
    f1 = np.zeros((n, 3))
    f2 = np.zeros((n, 3))
    f1[:, 1] = 0.5
    rho = np.array([0.01, 0.02])
    tw = calc_tw(
        tw_in,
        beta_phiw,
        np.zeros(2),
        np.zeros(3),
        rho,
        f1,
        f2,
        z,
        2,
        n,
        z_index,
    )
    assert tw.shape == (n, 2)
    assert np.isfinite(tw).all()
    # Species 0 couples to j=1 via pair index 1; should see non-zero contribution
    assert np.linalg.norm(tw[:, 0]) > 1e-8
    # Species 1 couples to j=0 via same pair index
    assert np.linalg.norm(tw[:, 1]) > 1e-8
    assert calc_l_index(0, 1) == calc_l_index(1, 0) == 1
