# DRAFT — Before and after: modernising the PhD fluid code

**Status:** working draft — fill in as phases complete; publish when the
modernisation arc (tooling → tests → refactor → bulk → fluids → Marimo) is
in a shape worth comparing end-to-end.

**Audience:** readers of the progressive notes who want a single comparative
view of *what changed and why*, not a changelog dump.

**Git anchors (Phase 0/1 baseline):**

| Label | Ref |
|-------|-----|
| Before | `main` at pre-modernisation tip (parent of PR #22), e.g. `879767c` |
| After (tooling+tests) | `feature/phase-0-1-tooling-and-tests` / PR #22 |
| After (final) | _TBD — tip of main once remaining phases land_ |

Use `git show <ref>:path` when refreshing snippets so the post stays honest.

---

## 1. Story arc

| Theme | Before | After (target) |
|-------|--------|----------------|
| Lost Fortran | No gold master; confidence = “LJ seems to work in Streamlit” | Characterisation tests + physical invariants + digitised figures |
| Python env | `>=3.13`, mixed `requirements*.txt` + uv | Pinned **3.14**, **uv-only**, locked latest wheels |
| Quality | black / isort / flake8 / mypy | **prek + ruff + ty** |
| Science vs UI | `modelling.py` imports Streamlit `@cache_data` | **`charged_fluids`** package is UI-free |
| Tests | 3 helpers under `src/tests/` | Top-level `tests/`, residual/`calc_tw` coverage ~99% |
| Bulk c(r) | Hand-dropped files; charged tail wrong | **Phase 3:** `data/bulk/` artefacts + `charged_fluids.bulk` OZ (LJ ready; KCl experimental) |
| Fluids | LJ credible; KCl/water/2-2 incomplete | _Phase 4: staged readiness flags_ |
| Exploration UI | Streamlit only | _Phase 5: Marimo on same core_ |

---

## 2. Tooling and packaging

### Before (`main:pyproject.toml`)

- Project name `PhD-Thesis`, `requires-python = ">=3.13"`
- Streamlit in **core** dependencies
- Optional extras: black, isort, flake8, mypy
- No build backend / installable science package

### After (PR #22)

- `phd-thesis` 0.2.0, `requires-python = ">=3.14,<3.15"`
- Core deps: numpy, scipy, pandas, psutil, loguru only
- **`app`** group: streamlit, plotly, matplotlib, duckdb, `pyarrow>=26`
- **`dev`**: pytest, ruff, ty, prek
- hatchling package: `src/charged_fluids`
- `.python-version` tracked; `.pre-commit-config.yaml` for prek

**Why it matters:** unit tests and Marimo can import the physics without pulling
Streamlit/pyarrow; Docker uses `uv sync --group app`.

### Commands readers can run

```bash
git show main:pyproject.toml | head -40
git show HEAD:pyproject.toml | head -50
just quality   # ruff + ty + pytest
```

---

## 3. Architecture: where the equations live

### Before

```
src/modelling.py     # numerics + streamlit.cache_data
src/parameters.py     # hard-coded fluid dict (mutated on read)
src/pages/*.py       # star-imports; duplicated LJ/KCl orchestration
src/oo_refactoring/  # parallel partial rewrite
```

Opening lines of the old core:

```python
# main:src/modelling.py (excerpt)
from streamlit import cache_data
from utils import monitor_memory
```

### After (Phase 0/1)

```
src/charged_fluids/          # installable, UI-free science core
  modelling.py
  parameters.py
  numerics.py
  memory.py
src/modelling.py             # thin shim → charged_fluids
src/pages/*.py               # still Streamlit; shims keep imports working
docs/LEGACY.md               # oo_refactoring / bulk / pyoz deferred
```

```python
# HEAD:src/charged_fluids/modelling.py (excerpt)
"""This module is UI-free. Streamlit caching belongs at the app edge, not here."""
from charged_fluids.memory import monitor_memory
```

**Why it matters:** refactor and vectorisation can proceed behind tests;
Streamlit and future Marimo share one engine.

**Phase 2 update:** TOML is the sole parameter source; `prepare_fluid_run` owns
setup; Streamlit LJ/KCl pages are thin; `Model` holds `hw`/`tw`/`c_short`/`f1`/`f2`.

---

## 4. Tests: from three asserts to a safety net

### Before

- `src/tests/test_modelling.py` — ~25 lines
- `calc_beta`, trivial `calc_l_index`, KCl kappa ≈ 6.77
- Nested under `src/`; awkward collection; env often broken on pyarrow builds

### After (Phase 1a+1b)

| Area | Examples |
|------|----------|
| Thermo helpers | beta, epsilon, rho, charge, kappa, pair index, non-mutating fluid factory |
| Potentials | LJ minimum at \(2^{1/6}\sigma\); wall decay; KCl Pauling / \(B\) |
| Bulk I/O | load `lj1-cr.data`; analytic f1 integral; hw Boltzmann factor |
| Residual | `opt_func ≡ tw_in - calc_tw`; 1- and 2-component neutral `calc_tw` |
| Solver smoke | small-grid LJ Newton–Krylov (`@pytest.mark.slow`) |
| Known gap | KCl short-range c(r) tail **xfail** until bulk fix |

Rough scale at PR #22: **~600+ lines** of tests vs **~25** before; package
coverage **~99%** on `charged_fluids`.

```bash
uv run pytest --cov=charged_fluids --cov-report=term-missing
```

---

## 5. Numerical core behaviour (call out carefully)

Document intentional behaviour locks vs cleanups:

| Topic | Before | After | Note |
|-------|--------|-------|------|
| Integrals | `np.trapz`; scalar broadcast onto `f1[k:]` | `scipy.integrate.trapezoid` on `f1[k]` | End-state equivalent; clearer |
| Fluid factory | Mutated global `fluid_parameters` dict | `deepcopy` per call | Tests isolated |
| Streamlit cache on c(r)/phiw | Yes | No in core | Cache at UI edge later |
| Pair index `i+j` | Historical | Unchanged (characterised) | Revisit before LJ2/water |

Do **not** claim bit-identical Fortran parity — sources are lost.

---

## 6. Bulk inputs and multi-fluid status

_Fill after Phases 3–4._

| Fluid | Before | After target |
|-------|--------|--------------|
| LJ1 | Working with `lj1-cr.data` | Reference + tests |
| LJ2 | Params only; missing bulk file | Bulk artefact + tests |
| KCl | Page exists; c(r) tail not → 0 | Fixed bulk short-range |
| 2-2 / water | `TO_BE_DEFINED.dat` | Staged experimental |

---

## 7. Developer experience

| | Before | After |
|--|--------|-------|
| Run app | `just app` / streamlit | same, via uv |
| Test | fragile / sparse | `just test` / `just quality` |
| Format/lint | black + flake8 | `just format` / ruff |
| Types | mypy | ty |
| Hooks | none project-local | `just hooks-install` (prek) |
| Blog trail | none | `docs/blog/` dated notes |

---

## 8. Suggested figures for the final post

1. Tree diagram: Streamlit/Marimo → `charged_fluids` → data/bulk  
2. Side-by-side import blocks (Streamlit coupling vs UI-free)  
3. Coverage bar or test-count before/after  
4. LJ \(g_{wi}(z)\) plot (app or Marimo) once reproducible from tests  
5. KCl c(r) tail plot illustrating the remaining scientific gap  

---

## 9. Checklist before publishing

- [ ] Refresh all `git show` snippets against final `main`
- [ ] Replace “target” rows with actual Phase 2–5 outcomes
- [ ] Link progressive notes (tooling, characterisation, calc_tw, …)
- [ ] Link PR(s) and production app https://phd.databooth.com.au
- [ ] Australian English pass; no overclaim on Fortran equivalence
- [ ] Optional: assemble into Databooth long-form from this draft + dated notes

---

## Related progressive notes

- [uv + Astral tooling](2026-11-19-uv-astral-tooling.md)
- [Characterisation tests](2026-11-19-characterisation-tests.md)
- [calc_tw coverage](2026-11-19-calc-tw-coverage.md)
