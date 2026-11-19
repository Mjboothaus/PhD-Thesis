"""Unit conventions for wall-fluid and bulk integral-equation paths.

Australian English in prose; SI base units unless noted.

## Length
- **Wall grid ``z`` / bulk radial ``r``**: angstrom (A).
- Bulk LJ ``sigma_A`` in TOML is also A (historical pyOZ used nm in some input
  files; generation code converts explicitly when needed).

## Energy and temperature
- Pair potentials ``u`` and wall potentials are in **joules** (J).
- ``beta = 1 / (k_B T)`` with ``T`` in kelvin, so ``beta * u`` is dimensionless.
- TOML LJ ``epsilon_kJ_mol`` is converted to J via ``1e3 / N_A`` in bulk generation.

## Charges and Coulomb
- ``charge = valence * e`` (coulombs).
- ``epsilon = 4 pi eps_r eps_0`` (SI), with distances still carried in A, so Coulomb
  pieces include an explicit ``1e-10`` (m/A) factor - see ``calc_u`` and
  ``coulomb_potential``.
- Wall/bulk dimensionless Coulomb scale::

      charge_pair[l] = 2e10 * beta * q_i * q_j / epsilon
      beta * u_c(r)  = charge_pair[l] / (2 r)     # r in A

  The factor 2 is historical bookkeeping so that wall integrals written with
  ``charge_pair`` match thesis formulae that absorbed a factor of two.

## Density and concentration
- TOML ``concentration`` is mol/dm3 (M).
- ``rho = concentration * N_A / 1e27`` -> number density per A^3.

## Bulk short-range direct correlation
- Neutral LJ: ``c_short = c_full``.
- Ionic (Ng split): ``c_short = c_s + beta u_c_short`` (= ``c_full + beta u_c``),
  which vanishes at large ``r`` when ``c_full -> -beta u_c``.

## Pair indexing
- Pair arrays use upper-triangle order including the diagonal via
  :func:`charged_fluids.modelling.calc_l_index` (symmetric in ``i, j``).
  For one- and two-component fluids this matches the historical ``i + j`` layout
  used by thesis data files.

## Wall convolution kernels
- Charge-path integrals (int hw dz) use the physical ``z`` grid (angstrom).
- The hw*f1 trapezoid pieces in ``calc_tw`` follow the historical thesis code and
  use **unit sample spacing** (``trapezoid(y)`` without ``z``). Scaling with delta z
  is absorbed into the discrete formulation as written in 1998.
"""

from __future__ import annotations

# Factors used when converting between laboratory and model units.
ANGSTROM_TO_METRE = 1.0e-10
LITRE_TO_ANGSTROM3 = 1.0e27  # 1 dm3 = 1e27 A^3
KJ_MOL_TO_J = 1.0e3  # divide by N_A at call site
