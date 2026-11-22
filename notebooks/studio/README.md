# Marimo Studio views (presentation pilot)

Australian English. Views for `notebooks/wall_fluid_presentation.py` only.

| View | Format | URL path (when running) |
|------|--------|-------------------------|
| **explainer** (default) | Vanilla HTML | `/` or `/explainer/` |
| **lecture** | React + Reveal.js slides | `/lecture/` |

## Run

```bash
# From repo root (PYTHONPATH=src via justfile)
just marimo-presentation
```

In the marimo editor, open the **Studio** pane and select **Explainer** or **Lecture**.
Preview rebuilds when view source is saved.

Validate bindings without the UI:

```bash
just marimo-studio-validate
```

## Design rules

- **Notebook** holds narrative markdown + Fig 3.2 Plotly (and any shared figures).
- **Views** hold layout only — no Python.
- Heavy interactive wall/bulk **solves** stay in `wall_fluid_explorer.py` / Streamlit.
- Production remains Streamlit on Render.

## Layout on disk

```text
notebooks/
  wall_fluid_presentation.py
  studio/wall_fluid_presentation/
    explainer/     # index.html, view.toml
    lecture/       # React Reveal deck
```

Build artefacts under `notebooks/__marimo__/studio/` are gitignored.
