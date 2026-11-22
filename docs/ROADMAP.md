# Roadmap — Charged fluids near interfaces

Australian English. Living plan for the modern reimplementation of the 1998 PhD
(*Charged fluids near interfaces: Integral equation theory*).

**Last updated:** 2026-11-22  
**Production:** https://phd-thesis-mjb.onrender.com (Render)  
**Physics detail:** [`PHYSICS.md`](PHYSICS.md) · **Fluids how-to:** [`FLUIDS.md`](FLUIDS.md) · **Legacy:** [`LEGACY.md`](LEGACY.md)

This workstream is **parked** after the physics/docs arc through PR #42, except
for the Marimo Studio note below (evaluation complete; no install).

---

## Where we are (snapshot)

### Done and trusted enough for demos

| Area | Status |
|------|--------|
| **uv + Just** toolchain, pytest, ruff/ty/prek | Stable |
| **`charged_fluids` core** (wall kernels, TOML params, `prepare_fluid_run`) | Stable |
| **Bulk OZ pipeline** → `data/bulk/*.cr.dat` (LJ, ionic Ng, CF1 anneal) | Usable |
| **Streamlit multipage** + Render Docker deploy | Live |
| **Marimo explorers** (wall + bulk, Plotly) | Live (`just marimo` / `marimo-bulk`) |
| **Fluid coverage** lj1, lj2, kcl, h2o, 1_1, 2_2 | All bulk- and wall-wired |
| **KCl bulk** Ng/PY + T anneal | Converges; soft LJ+Coulomb core |
| **Fig 3.2 digitised** wall–O/H (ψ₀=0) | In-repo CSVs + comparison API |
| **Legacy cleanup** oo_refactoring, empty packages, pyOZ mirror, `src/bulk/` | Removed |
| **Solver residual history** `ConvergenceHistory` | In core + Streamlit plot |

### Intentionally experimental

| Area | Why |
|------|-----|
| **h2o wall** | Newton–Krylov does not recover Fig 3.2 contact peaks (residual blow-up / g≈0–1) |
| **CF1 bulk** | Bare HNC under-associates (coord_OH ~1.2 vs ~2); no Duh–Haymet bridge |
| **Ionic bulk structure** | Soft LJ+Coulomb, not full Pauling core; semi-quantitative g(r) |
| **1_1 / 2_2** | Runnable RPM-style smoke; no full thermodynamics suite |
| **lj2 multi-pair c** | Equal multi-comp from 1-comp OZ (pair c copied) |

Architecture (one line):

```text
TOML → bulk Picard OZ (HNC/PY, Ng, optional T anneal) → c_short artefact
     → prepare_fluid_run → Newton–Krylov F[t_w]=0 → g_wi(z)
     → Streamlit / Marimo
```

---

## Open items (priority order)

### P0 — CF1 wall vs thesis Fig 3.2

**Goal:** Live wall solve at ψ₀ = 0 reproduces digitised contact peaks (~7 O, ~8.5 H near 1 Å).

| Step | Notes |
|------|--------|
| Diagnose residual blow-up | `calc_tw` / charge path / wall φ scale vs thesis FORTRAN |
| Fix wall convergence | Possibly softer wall, better initial guess, mixing, or bulk input scale |
| Regression | `compare_to_fig32` already in package; turn into hard assert when green |
| Bulk held fixed | Thesis practice: same bulk c(r) for all CF1 wall cases |

**Assets ready:** `data/Fig3.2.jpg`, `Fig3.2_w{O,H}_digitised.csv`, `load_fig32_profiles`, overlay PNG.

### P1 — CF1 bulk bridge IE

**Goal:** Intramolecular coordination closer to molecules intact (coord_OH ≈ 2, coord_HH ≈ 1).

| Step | Notes |
|------|--------|
| Literature bridge | Duh–Haymet JCP 1995; Ichiye–Haymet 1988 |
| Implement B(r) in bulk OZ | Beyond bare HNC; may need density-dependent form |
| Re-anneal + artefact | Regenerate `h2o.cr.dat`; re-test wall with new bulk |

### P2 — Ionic bulk realism

**Goal:** Melt/electrolyte g(r) closer to MC / full potential.

| Step | Notes |
|------|--------|
| Pauling + dispersion core in bulk OZ | Match wall `calc_u` ionic path, not only soft LJ |
| KCl g(r) vs published melt structure | Residual is already small with T anneal |
| 1_1 / 2_2 free energy / pressure routes | Beyond c_short artefacts |

### P3 — Benchmarks and hygiene

| Item | Notes |
|------|--------|
| Hansen–McDonald LJ tables | Automated regression for lj1 bulk/wall |
| Modern closures | e.g. Rogers–Young, Pihlajamaa–Janssen-style options in registry |
| Collapse Streamlit shims | `src/modelling.py` etc. → import `charged_fluids` only |
| Data hygiene | Archive `nrcg-cr-approx-maybe.dat`, optional drop `fluid_parameters.json` / legacy `lj1-cr.data` |
| Defer or archive `src/pyoz/` | Once bulk trust is high enough |
| Render custom domain | Optional; Cloud Run already deprecated |

---

## Parked / evaluated (no active work)

### Marimo Studio

**Status:** Evaluated; **not adopted**. Plain Marimo explorers remain the research UI.

| Option | Verdict |
|--------|---------|
| Replace Streamlit production with Studio | **No** — Render stays Streamlit |
| Day-to-day science / OZ solves | Keep **plain Marimo** (`just marimo`) + local SciPy |
| Teaching / thesis narrative (HTML, slides, PDF) | **Optional later** — good fit for presentation views only |
| molab / pure WASM export | Weak (needs `data/bulk` + installed package) |
| marimohub multi-user | Overkill unless a shared hub is required |

**If revisited:** pilot only — `uv add --group notebook 'marimo-studio[deno]'`, one explainer view on `notebooks/wall_fluid_explorer.py`, views under e.g. `studio/`. Do not migrate production.

### Other parked cleanups

- Blog long-form already has capstone + Warp cost post (`docs/blog/`)
- `WARP.md` may lag; prefer README + PHYSICS + this roadmap

---

## Suggested sequencing (when un-parking)

```text
1. CF1 wall residual diagnosis  ──►  Fig 3.2 regression green
2. CF1 bulk bridge              ──►  re-bulk h2o, re-check wall
3. Pauling ionic bulk           ──►  KCl / electrolyte structure
4. LJ table benchmarks          ──►  “ready” confidence band
5. Shim collapse + data hygiene ──►  thinner repo
6. Optional Marimo Studio pilot ──►  presentation only
```

Parallel-safe anytime: docs, Render tweaks, shim collapse, data hygiene (no physics dependency).

---

## Potential stretch ideas (not committed)

- **RISM / site–site** water as alternate to CF1 central force  
- **3D wall** or patterned electrodes (far beyond thesis scope)  
- **Machine-learned bridge** functions trained on simulation  
- **Interactive thesis companion** (Marimo Studio slides + theory page)  
- **Published package** on PyPI (`charged-fluids`) with minimal deps  
- **CI artefact cache** for bulk c(r) so PRs need not re-anneal CF1  
- **Digitise more thesis figures** (other ψ₀, KCl wall, bulk g(r)) for regression suite  
- **Fortran reconstruction notes** blog series from characterisation tests  

---

## How to use this doc

- **Starting a sprint:** pick one P0–P2 item; link the PR to this heading  
- **Physics truth:** update [`PHYSICS.md`](PHYSICS.md) when methods or gaps change  
- **Parking again:** leave a one-line “parked as of DATE” under *Where we are*  
- **Do not expand** legacy trees (`LEGACY.md`)

---

## Recent milestone trail (for context)

| When | Milestone |
|------|-----------|
| 2026-11 | uv-only stack, characterisation tests, `charged_fluids` package |
| 2026-11 | Bulk c(r) CLI + artefacts; Streamlit on Render |
| 2026-11 | Marimo wall/bulk explorers; all six fluids wired |
| 2026-11 | CF1 live bulk OZ (T anneal); ionic bulk T anneal (KCl converges) |
| 2026-11 | Fig 3.2 digitised; live wall still fails contact peaks |
| 2026-11 | oo_refactoring / src/bulk removed; PHYSICS.md + this roadmap |

---

*Parked workstream note: resume at **P0 CF1 wall vs Fig 3.2** unless priorities change.*
