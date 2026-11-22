# Physics status and methods

Australian English. This note tracks **what the code solves**, **what is trusted**,
and **known gaps** versus the 1998 thesis and literature.

Prioritised plans and parked options: [`ROADMAP.md`](ROADMAP.md).

## Pipeline

```text
TOML parameters
  → bulk OZ (Picard / Ng / optional T anneal) → data/bulk/{symbol}.cr.dat
  → prepare_fluid_run (interpolate c_short, f1, f2, wall φ)
  → Newton–Krylov singlet residual F[t_w]=0
  → g_wi(z) = h_wi(z)+1
```

Wall Coulomb scale matches bulk via `charge_pair` (`β u_c = charge_pair/(2r)`).
Ng short-range bulk input is `c_short = c_s + β u_c,short` so the large-r tail
vanishes when `c → −β u_c`.

## Fluid matrix (2026-11)

| Symbol | Bulk method | Wall | Status | Notes |
|--------|-------------|------|--------|-------|
| **lj1** | 1-comp HNC LJ | NK | **ready** | Demo quality |
| **lj2** | equal multi-comp from 1-comp OZ | NK | **ready** | Pair c copied across pairs |
| **kcl** | Ng/PY + **T anneal** (soft LJ core) | NK | experimental | Melt density; bulk residual now small with anneal |
| **h2o** | CF1 Ng/HNC + **T anneal** | NK | experimental | Bare HNC under-associates (coord_OH ≲ 2); bridge IE open |
| **1_1** | Ng/HNC + T anneal, soft RPM core | NK | experimental | Aqueous 1:1 continuum solvent (εᵣ=78.3) |
| **2_2** | Ng/PY + T anneal, soft RPM core | NK | experimental | Aqueous 2:2 continuum solvent |

Bulk artefacts: `data/bulk/{symbol}.cr.dat` (+ `.meta.json`). Regenerate with
`just bulk-cr <symbol>` or `just bulk-cr-all`.

## Bulk OZ methods

### Closures

- **HNC:** `c = exp(−βu+γ) − 1 − γ`
- **PY:** `c = (e^{−βu}−1)(1+γ)`
- Exp arguments are clipped for numerical stability (association / deep wells).

### Ng Coulomb split

Long-range piece `β u_long ∝ charge_pair · erf(α r)/(2 r)`; short-range potential
is iterated; wall input uses thesis-style short-range `c`.

### Temperature annealing

Dense ionic melts and CF1 association are stiff for plain Picard at the target
state point. Geometric cool-down (hot → target) with continuous `γ` (rescaled by
`β`) plus optional Anderson mixing is used for:

- CF1 water (`solve_bulk_cf1`)
- Binary ionic + LJ (`solve_bulk_binary_coulomb_lj` when `anneal_steps > 1`)

### Ionic bulk model limitation

Ionic bulk currently uses **soft LJ + Coulomb**, not the full Pauling
exponential + dispersion core from the wall `calc_u` path. That is intentional
for Ng/OZ stability; melt structure is therefore semi-quantitative.

## Wall singlet

- Residual: `F[t_w] = t_w − calc_tw(...)` (multi-component form).
- Solver: SciPy Newton–Krylov (`solve_model`), optional `ConvergenceHistory`.
- Diagnostics: `wall_solution_report` / `contact_and_bulk_g` (contact, peak, g(∞)).

Hard-wall-like contact (`g≈0` at z=0) is expected for the historical wall
potential form used here.

## Thesis / literature gaps

| Topic | Thesis / lit | This stack |
|-------|--------------|------------|
| CF1 bulk | Duh–Haymet HNC **+ bridge** (JCP 1995) | HNC only + anneal; coord_OH ~1.2 vs ~2 |
| CF1 wall Fig 3.2 | Digitised wall–O/H at ψ₀=0 (peaks ~7 / ~8.5 near 1 Å) | **Digitised in-repo** (`data/Fig3.2_*`); live NK wall does **not** yet recover contact peaks (residual blow-up / g≈0–1) |
| CF1 bulk Fig 3.12 | Digitised c(r) | Live OZ default; digitised opt-in |
| KCl melt g(r) | Full ionic potential / MC | Soft LJ+Coulomb Ng/PY; improved residual with anneal |
| 1-1 / 2-2 | Soft-core PM in continuum water | Runnable Ng bulk + wall smoke; not full free-energy suite |
| LJ benchmarks | Hansen–McDonald tables | Qualitative; no automated table regression |

Bridge functions for CF1 (and better ionic cores) remain the main physics
extension. Fig 3.2 is digitised for regression; matching it requires a working
CF1 wall residual (and likely bulk bridge), not only bulk c_short.

### Fig 3.2 assets

- Source scan: `data/Fig3.2.jpg` (Case 1, ψ₀ = 0 mV)
- Curves: `data/Fig3.2_wO_digitised.csv`, `data/Fig3.2_wH_digitised.csv`
- Manual knots: `data/Fig3.2_*_manual_points.csv`
- Overlay: `data/Fig3.2_comparison.png` (thesis vs current live solve)
- API: `load_fig32_profiles`, `compare_to_fig32`

Thesis bulk c(r) is held **constant** across CF1 wall calculations (same bulk
for all ψ₀ / wall cases).

## Validation commands

```bash
just bulk-cr kcl          # annealed Ng/PY melt bulk
just bulk-cr h2o          # CF1 anneal (slow at full grid)
just test                 # unit + smoke
uv run python -m charged_fluids.bulk.cli --list
```

## References (selected)

- Booth, PhD thesis (1998), *Charged fluids near interfaces*
- Duh, Perera & Haymet, J. Chem. Phys. **102**, 3736 (1995) — CF1 bulk IE
- Ichiye & Haymet, J. Chem. Phys. **89**, 4315 (1988) — CF water
- Ng, J. Chem. Phys. **61**, 2680 (1974) — long-range renormalisation
- Hansen & McDonald, *Theory of Simple Liquids*
