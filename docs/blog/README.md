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
| 2026-11-19 | [Phase 2 — core package depth and thinner Streamlit pages](2026-11-19-phase-2-core-ui.md) |

## Capstone comparison (in progress)

At the end of the modernisation arc, publish a single comparative post built from
[DRAFT-before-after-comparison.md](DRAFT-before-after-comparison.md). It is kept
updated with git anchors (`main` vs feature tips) so before/after claims stay
verifiable. Progressive notes below feed that draft; they are not a substitute
for it.

## Themes

- Lost Fortran 77 sources and what “correct” means without a gold master
- Integral equation theory (singlet OZ + bulk c(r) inputs)
- uv-only Python management; prek + ruff + ty
- Characterisation tests before refactor
- Bulk PyOZ pipeline as a first-class artefact
- Streamlit for deployment, Marimo for research exploration
- **Before/after comparison of key codebase elements** (capstone post)
