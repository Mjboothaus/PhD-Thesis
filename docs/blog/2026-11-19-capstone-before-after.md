# Before and after: modernising a lost-Fortran PhD into a tested Python stack

**Date:** 2026-11-19  
**Status:** Capstone comparison (post Phases 0–5)  
**App:** [https://phd.databooth.com.au](https://phd.databooth.com.au)  
**PRs:** [#22](https://github.com/Mjboothaus/PhD-Thesis/pull/22) · [#23](https://github.com/Mjboothaus/PhD-Thesis/pull/23) · [#24](https://github.com/Mjboothaus/PhD-Thesis/pull/24) · [#25](https://github.com/Mjboothaus/PhD-Thesis/pull/25)

---

## Why this post exists

In 1998 I submitted *Charged fluids near interfaces: Integral equation theory*.
The numerical engine was Fortran 77. Those sources are gone. A Python/Streamlit
reimplementation grew over the years as a useful demo — but with almost no
tests, mixed packaging, and charged-fluid bulk inputs that never quite
behaved.

This post is not a changelog. It is a **before/after** of the codebase after a
deliberate modernisation arc: tooling → tests → package boundaries → bulk c(r)
pipeline → multi-fluid readiness → Marimo.

**Git anchors**

| Label | Ref |
|-------|-----|
| Before | `879767c` (pre-modernisation `main`) |
| After | `852b028` (`main` after PR #25) |

```bash
git show 879767c:pyproject.toml | head -40
git show 852b028:pyproject.toml | head -50
```

---

## 1. Story arc at a glance

| Theme | Before | After |
|-------|--------|-------|
| Confidence without Fortran | “LJ seems to work in the app” | Characterisation tests, physical invariants, explicit xfails |
| Python / deps | `>=3.13`, `requirements*.txt` + uv | **Python 3.14**, **uv-only**, dependency groups |
| Quality | black / isort / flake8 / mypy | **prek + ruff + ty** |
| Science vs UI | `modelling.py` imports Streamlit `@cache_data` | **`charged_fluids`** package is UI-free |
| Tests | ~3 helpers under `src/tests/` | Top-level `tests/`, ~99% coverage on core |
| Bulk c(r) | Hand-dropped files; KCl short-range dubious | `data/bulk/` artefacts + `charged_fluids.bulk` OZ |
| Fluids | LJ1 credible | **LJ1 + LJ2 ready**; KCl experimental; water/2-2 stubs |
| Exploration | Streamlit only | Streamlit **and** Marimo on the same core |

---

## 2. Tooling and packaging

### Before

- Streamlit sat in **core** dependencies (pulling pyarrow into every install).
- Dev extras: black, isort, flake8, mypy.
- No installable science package; pages relied on `PYTHONPATH=src`.

### After

- Core: numpy, scipy, pandas, psutil, loguru.
- **`app`**: streamlit, plotly, matplotlib, duckdb, `pyarrow>=26` (binary wheels on 3.14).
- **`dev`**: pytest, ruff, ty, prek.
- **`notebook`**: marimo.
- hatchling builds `src/charged_fluids`.

**Why it matters:** unit tests and notebooks never need Streamlit. Docker runs
`uv sync --group app`. One lockfile is the source of truth.

```bash
just quality   # ruff + ty + pytest
just app
just marimo
just bulk-cr lj1
```

---

## 3. Architecture: where the equations live

### Before

```
src/modelling.py     # numerics + streamlit.cache_data
src/parameters.py     # hard-coded dict, mutated on read
src/pages/*.py       # star-imports; duplicated LJ/KCl setup
src/oo_refactoring/  # parallel half-finished rewrite
```

```python
# 879767c:src/modelling.py
from streamlit import cache_data
from utils import monitor_memory
```

### After

```
src/charged_fluids/
  modelling.py      # UI-free singlet OZ helpers
  parameters.py     # TOML loader
  setup.py          # prepare_fluid_run
  bulk/             # bulk OZ + artefacts
src/components/fluid_page.py   # shared Streamlit shell
notebooks/wall_fluid_explorer.py
docs/LEGACY.md
```

```python
# charged_fluids.modelling
"""This module is UI-free. Streamlit caching belongs at the app edge."""
from charged_fluids.memory import monitor_memory
```

**`Model` is consistent** (old issue #21): `hw`, `tw`, `c_short`, `f1`, `f2`
all live on one dataclass. **`data/fluid_parameters.toml`** is the sole
parameter source (issue #16).

---

## 4. Tests: from three asserts to a safety net

| Before | After |
|--------|-------|
| `calc_beta`, trivial pair index, KCl κ ≈ 6.77 | Thermo, potentials, integrals, residual identity, LJ smoke solve |
| Nested under `src/`; pyarrow install pain | Top-level `tests/`; core install without Streamlit |
| No bulk generation tests | Bulk sine transforms, LJ generate, KCl c_short **tail** |

Scale: on the order of **hundreds of lines** of tests vs **~25** before;
`charged_fluids` line coverage **~99%** (excluding experimental ionic residual
paths).

```bash
uv run pytest --cov=charged_fluids --cov-report=term-missing
```

We still **do not** claim bit-identical Fortran parity. The safety net is
algebraic and physical: signs, limits, analytic checks, and regression markers.

---

## 5. Bulk c(r): the quiet dependency of the whole thesis

Wall–fluid integral equations need bulk direct correlations c(r) as input.
Historically those arrived as opaque data files.

### After Phase 3–5

| Fluid | Bulk path | Notes |
|-------|-----------|-------|
| **lj1** | `data/bulk/lj1.cr.dat` | Picard HNC; converges; tail ~1e-5 |
| **lj2** | `data/bulk/lj2.cr.dat` | Equal mixture; pair columns from 1-comp OZ at total density |
| **kcl** | `data/bulk/kcl.cr.dat` | **Ng + PY**; **c_short tail ~1e-6**; Picard residual still elevated / g(r) not fully physical |
| h2o / 2_2 | missing | Stubs only |

```bash
just bulk-cr lj1
just bulk-cr kcl
```

`prepare_fluid_run` prefers `data/bulk/{symbol}.cr.dat` over legacy
`data/*.data`.

### KCl stabilisation (honest account)

The original gap: short-range charged c(r) should → 0 at large r. A naive
multi-comp HNC at melt density **diverged** when Coulomb was scaled correctly.

What worked:

1. **Dimensionless Coulomb** must use the wall-fluid `charge_pair` scale  
   (`βu_c = charge_pair / (2r)`), not a bare SI 1/r that underflows to ~0.
2. **Ng split**: iterate on short-ranged potential  
   `βu_s = βu_LJ + βu_c·erfc(αr)/…` with γ initialised from the long-range piece.
3. **PY closure** + small mixing at melt density yields a **finite c_short with
   vanishing tail**, even when the Picard residual has not yet hit a tight tol
   and g(r) may still show unphysical negatives.

So: the **pipeline requirement for the wall code** (short-range c tail) is met
for KCl artefacts; full thermodynamic convergence of bulk g(r) at melt density
remains **experimental** (issues #13 / #17). That is progress you can measure,
not a claim that the salt problem is “solved”.

---

## 6. Two UIs, one engine

| | Streamlit | Marimo |
|--|-----------|--------|
| Entry | `just app` | `just marimo` |
| Role | Deployed demo (Cloud Run) | Research exploration |
| Code | `render_fluid_page("lj1")` etc. | `prepare_fluid_run` + `solve_model` |

Neither UI owns the physics. That was the point of Phase 0–2.

---

## 7. Developer experience

| | Before | After |
|--|--------|-------|
| Run app | streamlit / just | `just app` via uv |
| Test | sparse / fragile | `just test` / `just quality` |
| Format / lint | black + flake8 | ruff |
| Types | mypy | ty |
| Hooks | none | prek |
| Bulk inputs | tribal knowledge | `just bulk-cr` + metadata JSON |
| Docs trail | README TODOs | `docs/blog/` progressive notes + this capstone |

---

## 8. What we deliberately did *not* do

- Rewrite bulk OZ to match bit-exact pyoz/Fortran (sources incomplete; units
  archaeology continues).
- Delete `src/pyoz` / `src/bulk` / `oo_refactoring` yet (see `docs/LEGACY.md`).
- Implement water or 2-2 electrolytes end-to-end.
- Pretend KCl melt HNC is production-ready.

---

## 9. Figures worth adding when you publish externally

1. Architecture: Streamlit/Marimo → `charged_fluids` → `data/bulk`.
2. Side-by-side imports (Streamlit-coupled vs UI-free).
3. Test-count or coverage before/after.
4. LJ `g_wi(z)` from Marimo or the app.
5. KCl `|c_short(r)|` tail: legacy file vs Ng artefact.

---

## Related progressive notes

- [uv + Astral tooling](2026-11-19-uv-astral-tooling.md)
- [Characterisation tests](2026-11-19-characterisation-tests.md)
- [calc_tw coverage](2026-11-19-calc-tw-coverage.md)
- [Phase 2 core/UI](2026-11-19-phase-2-core-ui.md)
- [Phase 3 bulk c(r)](2026-11-19-phase-3-bulk-cr.md)
- [Phase 4–5 LJ2 + Marimo](2026-11-19-phase-4-5-fluids-marimo.md)
- Issue triage: [`docs/ISSUES_TRIAGE.md`](../ISSUES_TRIAGE.md)

---

*Australian English in documentation. Code identifiers keep US spelling.*
