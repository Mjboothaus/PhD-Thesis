# Blog notes — PhD thesis modernisation

Progressive notes captured while modernising the 1998 *Charged fluids near
interfaces* work as a Python project. Each phase/PR should add a dated stub;
a long-form post can be assembled later for Databooth or a personal site.

## Index

| Date | Note |
|------|------|
| 2026-11-19 | [From mixed requirements to uv-only + Astral quality gates](2026-11-19-uv-astral-tooling.md) |
| 2026-11-19 | [Rebuilding confidence without the Fortran](2026-11-19-characterisation-tests.md) |
| 2026-11-19 | [Phase 1b — residual coverage without a Fortran oracle](2026-11-19-calc-tw-coverage.md) |

## Themes

- Lost Fortran 77 sources and what “correct” means without a gold master
- Integral equation theory (singlet OZ + bulk c(r) inputs)
- uv-only Python management; prek + ruff + ty
- Characterisation tests before refactor
- Bulk PyOZ pipeline as a first-class artefact
- Streamlit for deployment, Marimo for research exploration
