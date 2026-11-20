"""CF1 central-force model of water (Haymet / Stillinger-Rahman family).

Pair potentials follow Duh, Perera & Haymet, J. Chem. Phys. 102, 3736 (1995)
and the CF1 forms used in later MC studies (e.g. CMP 24, 33503 (2021)):

    U_OO, U_HH, U_OH in kcal/mol with r in angstrom.

Partial charges: q_H = 0.32983 e, q_O = -2 q_H.

Component order in this package is always **O then H** (pair index OO, OH, HH).
"""

from __future__ import annotations

import numpy as np
from scipy.constants import Avogadro, calorie

# Elementary-charge fractions (dimensionless)
CF1_Q_H = 0.32983
CF1_Q_O = -2.0 * CF1_Q_H

# kcal/mol -> joule (per pair interaction, not per mole of molecules)
_KCAL_MOL_TO_J = (calorie * 1.0e3) / Avogadro


def _u_oo_kcal(r: np.ndarray) -> np.ndarray:
    r = np.asarray(r, dtype=float)
    out = np.full_like(r, np.inf, dtype=float)
    m = r > 1.0e-8
    rm = r[m]
    out[m] = (
        144.538 / rm
        + 24082.38 / rm**8.8591
        - 0.25 * np.exp(-4.0 * (rm - 3.4) ** 2)
        - 0.25 * np.exp(-1.5 * (rm - 4.5) ** 2)
    )
    return out


def _u_hh_kcal(r: np.ndarray) -> np.ndarray:
    r = np.asarray(r, dtype=float)
    out = np.full_like(r, np.inf, dtype=float)
    m = r > 1.0e-8
    rm = r[m]
    exp_hh = np.exp(np.clip(40.0 * (rm - 2.0), -60.0, 60.0))
    out[m] = (
        36.1345 / rm
        + 18.0 / (1.0 + exp_hh)
        - 17.0 * np.exp(np.clip(-7.62177 * (rm - 1.45251) ** 2, -60.0, 0.0))
    )
    return out


def _u_oh_kcal(r: np.ndarray) -> np.ndarray:
    r = np.asarray(r, dtype=float)
    out = np.full_like(r, np.inf, dtype=float)
    m = r > 1.0e-8
    rm = r[m]
    exp_oh1 = np.exp(np.clip(40.0 * (rm - 1.05), -60.0, 60.0))
    exp_oh2 = np.exp(np.clip(5.49305 * (rm - 2.2), -60.0, 60.0))
    out[m] = (
        -72.269 / rm
        + 6.23403 / rm**9.19912
        - 10.0 / (1.0 + exp_oh1)
        - 4.0 / (1.0 + exp_oh2)
    )
    return out


def cf1_pair_potential_kcal(r: np.ndarray) -> np.ndarray:
    """Return U_ij(r) in kcal/mol with shape (n_r, 3) for pairs OO, OH, HH."""
    r = np.asarray(r, dtype=float)
    return np.column_stack([_u_oo_kcal(r), _u_oh_kcal(r), _u_hh_kcal(r)])


def calc_u_cf1(n_point: int, r: np.ndarray) -> np.ndarray:
    """CF1 pair potential on grid ``r`` in joules; shape (n_point, 3)."""
    del n_point
    u_kcal = cf1_pair_potential_kcal(r)
    u = u_kcal * _KCAL_MOL_TO_J
    # Match other potentials: copy first finite interior point to r=0 slot if needed
    if r.size and r[0] == 0.0:
        u[0, :] = u[1, :] if r.size > 1 else np.inf
    return u
