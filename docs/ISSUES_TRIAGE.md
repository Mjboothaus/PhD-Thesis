# GitHub issues triage (2026-11-19)

Open issues on https://github.com/Mjboothaus/PhD-Thesis/issues were last
touched ~4 years ago. Status relative to the 2026 modernisation (PR #22+).

| # | Title | Verdict | Notes |
|---|-------|---------|-------|
| 21 | Treat model fields consistently (f1/f2 vs hw/tw) | **Addressed in Phase 2** | `Model` now holds `hw`, `tw`, `c_short`, `f1`, `f2`; setup + solver write back |
| 20 | Implement pytest | **Done (Phase 1)** | Top-level `tests/`, ~99% `charged_fluids` coverage |
| 19 | Update SciPy 1.8→1.9 | **Superseded** | Lock has SciPy 1.18.x on Python 3.14 |
| 18 | Cloud Run connect error 60000 | **Unknown / retest** | Needs fresh deploy after Docker 3.14 image; not verified this sprint |
| 17 | Bulk KCl at 2 state points | **Still open (Phase 3)** | Needs correct c_sr(r) generation |
| 16 | Parameters in TOML/JSON | **Done (Phase 2)** | `data/fluid_parameters.toml` is sole source |
| 14 | Max iterations vs reported nits | **Improved (Phase 2)** | UI shows `nit` explicitly from SciPy result |
| 13 | Get KCl working (c_sr formula) | **Still open (Phase 3)** | Same bulk long-range separation as #17; xfail test tracks tail |
| 11 | Mobile experience | **Still open / low** | UX polish; not blocking science |
| 10 | Streamlit tabs jumping back | **Still open / low** | Likely Streamlit session state; retest on 1.67 |
| 9 | Chart legend / axis labels | **Partial** | Bulk plots use clearer names; more polish welcome |
| 8 | README + in-app instructions | **Partial** | README updated in Phase 0; in-app copy still thin |
| 6 | S3 or on-the-fly pyoz | **Still open (Phase 3)** | Prefer on-the-fly/cached bulk artefacts over S3 |
| 5 | Convergence y-axis labels | **wontfix / low** | Already labelled; Plotly quirks |
| 4 | More thesis markdown pages | **Still open / docs** | Content backlog |
| 3 | beta*u bulk curves not visible | **Addressed in Phase 2** | LJ bulk tab plots `beta*u` with sensible ylim |

## Suggested GitHub actions

- Close: #19, #20, #16 (with pointers to PRs)
- Close or note fixed: #21, #3, #14
- Keep open with milestone "bulk-c(r)": #13, #17, #6
- Keep open backlog: #4, #8, #9, #10, #11, #18
- Leave #5 as wontfix
