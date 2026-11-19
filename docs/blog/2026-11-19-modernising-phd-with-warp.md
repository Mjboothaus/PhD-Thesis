# Modernising a lost-Fortran PhD with Warp: equations first, agents second

**Status:** Ready for Databooth publication  
**Date:** 2026-11-19  
**Suggested slug:** `modernising-phd-with-warp`  
**Suggested URL:** https://www.databooth.com.au/posts/modernising-phd-with-warp/  
**Tags:** AI-assisted development, Python, scientific computing, Warp, uv, Marimo, engineering practices  
**App:** [https://phd-thesis-mjb.onrender.com](https://phd-thesis-mjb.onrender.com) (Render; Cloud Run deprecated)  
**Repo:** [Mjboothaus/PhD-Thesis](https://github.com/Mjboothaus/PhD-Thesis)  
**PRs:** [#22](https://github.com/Mjboothaus/PhD-Thesis/pull/22)–[#32](https://github.com/Mjboothaus/PhD-Thesis/pull/32) (merged); [#33](https://github.com/Mjboothaus/PhD-Thesis/pull/33) bulk CLI; [#34](https://github.com/Mjboothaus/PhD-Thesis/pull/34) notebooks + this post

---

## Executive summary

- In 1998 I submitted *Charged fluids near interfaces: Integral equation theory*. The numerical engine was Fortran 77; those sources are gone. From at least **July 2022** (see this repository’s git history) I started a **new** Python reimplementation—not a salvage of leftover thesis code—building wall–fluid numerics, pyOZ experiments, and a Streamlit app from scratch.
- By the **2026 modernisation tip** that stack was a useful demo, but still thin on tests, mixed packaging, and bulk correlation inputs that did not quite behave for charged fluids. A concentrated modernisation arc then rebuilt it around a UI-free science core (`charged_fluids`), uv-only packaging on Python 3.14, Astral quality gates (prek + ruff + ty), a bulk Ornstein–Zernike pipeline, multi-fluid readiness, and research notebooks in Marimo with Plotly.
- Warp’s agentic environment accelerated the work: phased pull requests, characterisation tests before refactors, bulk CLI design, KCl short-range stabilisation, and notebook consolidation. Humans still owned physics judgement, merge discipline, and “do not overclaim Fortran parity”.
- **Cost of Warp assistance is real and metered.** On current Warp pricing, agent usage is charged at provider API rates (plus optional cloud compute and platform charges). For a day of dense multi-PR agent work on a mid-size repository, a practical planning band is **on the order of tens of US dollars of included or add-on agent usage**—not “free coding”, and not a substitute for domain skill. Exact invoices vary with model choice, conversation length, and tool-call depth; check *Settings → Billing and usage* for your own run.

---

## Why this matters

Most AI coding stories stop at “we shipped faster”. Scientific and long-lived research code has a harder test:

1. **There may be no gold master.** Lost Fortran means a greenfield Python rewrite cannot diff against the original engine; you characterise behaviour and lock physical invariants.
2. **UI is not the product.** Streamlit is a demo surface on top of a multi-year reimplementation. The equations must live somewhere importable by tests, notebooks, and CLIs.
3. **Speed without guardrails is just faster entropy.** Agents will happily duplicate fluid selectors, leave xfail markers forever, and regenerate bulk metadata timestamps if you let them.

This project is a small case study in that combination: **intent-driven development with deterministic checks**, applied to integral equation theory rather than a CRUD app.

---

## Timeline (so the story stays honest)

| When | What |
|------|------|
| **1998** | PhD submitted; numerical work in **Fortran 77** |
| **Later** | F77 sources no longer available (not recoverable from the old online-resources path) |
| **2022-07 onward** | **New** Python implementation in this repository (`Initial setup`, pyOZ bring-up, Streamlit app, digitised thesis figures)—a rewrite, not residual Fortran glue |
| **Pre-2026 modernisation tip** | Working demo; sparse tests; packaging and UI still tangled (`879767c`) |
| **2026-11** | Phased modernisation PRs #22–#34 |

---

## The before picture (honest)

“Before” here means **immediately prior to the 2026 modernisation arc**—after years of the 2022+ Python rewrite—not “whatever was left of the thesis tarball”.

| Theme | Before (≈ `879767c`) |
|-------|--------|
| Confidence | “LJ seems to work in the app” |
| Packaging | `requirements*.txt` + uv; Streamlit in core dependencies |
| Quality | black / isort / flake8 / mypy |
| Science vs UI | `modelling.py` imported Streamlit `@cache_data` |
| Tests | ~3 helpers under `src/tests/` |
| Bulk c(r) | Hand-dropped files; KCl short-range tail dubious |
| Exploration | Streamlit only |

Git anchor for that snapshot: `879767c` on `main`.

---

## What we built (the after picture)

Architecture in one line:

> TOML parameters → bulk Picard OZ (HNC/PY, optional Ng) → wall singlet residual (Newton–Krylov) → Streamlit **or** Marimo

| Theme | After |
|-------|--------|
| Python | **3.14**, **uv-only**, dependency groups (`dev`, `app`, `notebook`) |
| Quality | **prek + ruff + ty**; `just quality` |
| Science core | installable **`src/charged_fluids/`** (UI-free) |
| Parameters | single source: `data/fluid_parameters.toml` |
| Bulk | `charged_fluids.bulk` + `data/bulk/{symbol}.cr.dat` (+ `.meta.json`) |
| Fluids | lj1/lj2 **ready**; kcl, h2o (CF1), 2_2, 1_1 **experimental** |
| UIs | production Streamlit + Marimo explorers (Plotly) |
| Tests | top-level `tests/`; full suite green (no permanent KCl xfail) |

Production demo: [phd-thesis-mjb.onrender.com](https://phd-thesis-mjb.onrender.com).

```bash
just init
just quality
just app                 # Streamlit
just marimo              # wall–fluid explorer (Plotly)
just marimo-bulk         # bulk artefacts / parameters
just bulk-cr lj1         # regenerate bulk c(r)
```

---

## Story arc in phases (not a changelog dump)

The work was deliberately **phased**, each landing as a pull request with a short progressive note under `docs/blog/`:

| Phase | Outcome | PR |
|-------|---------|-----|
| 0–1 | uv/Astral tooling, `charged_fluids`, characterisation tests | [#22](https://github.com/Mjboothaus/PhD-Thesis/pull/22) |
| 2 | TOML params, `prepare_fluid_run`, thinner Streamlit pages | [#23](https://github.com/Mjboothaus/PhD-Thesis/pull/23) |
| 3 | Bulk c(r) OZ pipeline and artefacts | [#24](https://github.com/Mjboothaus/PhD-Thesis/pull/24) |
| 4–5 | LJ2 ready + first Marimo explorer | [#25](https://github.com/Mjboothaus/PhD-Thesis/pull/25) |
| Capstone + KCl Ng/PY | honest bulk short-range tail; before/after write-up | [#26](https://github.com/Mjboothaus/PhD-Thesis/pull/26) |
| Hardening | bulk CLI TOML defaults, CF1/2-2/1-1, core `calc_tw`, Marimo display fixes | [#27](https://github.com/Mjboothaus/PhD-Thesis/pull/27)–[#32](https://github.com/Mjboothaus/PhD-Thesis/pull/32) |
| Bulk CLI polish | list/all/validate, wall-usable tails | [#33](https://github.com/Mjboothaus/PhD-Thesis/pull/33) (companion) |
| Notebooks + narrative | Plotly Marimo, single fluid selector, legacy Jupyter archived, modern KCl tests, this post | [#34](https://github.com/Mjboothaus/PhD-Thesis/pull/34) |

From the pre-modernisation tip through this arc: on the order of **~90 files touched**, **~9k insertions / ~3k deletions**—mostly package extraction, tests, bulk solver, and documentation, not vanity rewrites.

---

## Three technical decisions that actually mattered

### 1. UI-free science core

Streamlit caching does not belong next to Newton–Krylov residuals. Once `charged_fluids` stopped importing Streamlit, the same `prepare_fluid_run` / `solve_model` path could serve:

- unit tests without pyarrow pain,
- bulk CLI,
- Marimo research notebooks,
- thin Streamlit pages.

That is the hyper-team pattern applied to research code: **one engine, multiple surfaces**.

### 2. Characterisation before “fix the bulk file”

For KCl, the wall code needs short-range \(c(r)\) that **vanishes at large r**. A naive multi-component HNC at melt density diverged when Coulomb was scaled correctly. What worked:

1. dimensionless Coulomb on the thesis `charge_pair / (2r)` scale;
2. **Ng split** (iterate short-ranged potential; long-range in γ);
3. **PY + small mixing** → finite `c_short` with tail ~**1e-6** on shipped artefacts.

Honest caveat: Picard residual / \(g(r)\) at melt density can still be imperfect. Meeting the **pipeline requirement for the wall solver** is not the same as claiming the salt problem is “solved”.

We also **removed a permanent xfail** that still pointed at a legacy dump (`nrcg-cr-approx-maybe.dat`) and replaced it with a hard assert on `data/bulk/kcl.cr.dat`. Xfail markers that outlive the fix are how technical debt hides in plain sight.

### 3. Marimo as research UI, Streamlit as demo

Streamlit stays the deployed face of the project (**Render**; Cloud Run is deprecated). Marimo is where you poke grids, bulk OZ dry-runs, and wall solves with **one fluid selector** (no separate “bulk fluid” dropdown fighting the wall choice) and **Plotly** figures for interactive inspection—consistent with the app’s Plotly style.

Legacy Jupyter notebooks moved under `notebooks/legacy/` with a mapping table. Do not extend them; do not pretend they are the product.

---

## Guardrails: agents generate, policy enforces

This matches the stance in earlier Databooth writing on deterministic guardrails and “more than just the vibe”:

- **prek + ruff + ty + pytest** as non-optional gates (`just quality`);
- feature branches and small pull requests (standard git workflow);
- progressive `docs/blog/` notes so the narrative is reconstructable;
- Australian English in docs; US spelling in identifiers;
- explicit LEGACY trees rather than silent dual implementations.

Agents are excellent at drafting the bulk CLI, wiring Marimo cells, and proposing Ng mixing parameters. They are poor at knowing when a physical residual is “good enough” or when an xfail should die. That judgement stayed human.

---

## Indicative cost of using Warp on this work

### How Warp bills agent work (as of late 2025 / 2026)

Warp’s paid **Build** plan is commonly cited around **US$20/month** with **included agent usage metered at public API rates** (higher tiers include more allowance). Usage is not a flat “one prompt = one dollar”:

| Component | What it is |
|-----------|------------|
| **Inference** | Tokens in/out (and cache) at the chosen model’s API rates when Warp supplies the model |
| **Compute** | Warp-hosted sandboxes for **cloud** agents (local runs use your machine) |
| **Platform** | Billable agent time on some cloud / team configurations |

See Warp’s current [pricing](https://www.warp.dev/pricing) and [usage and billing](https://docs.warp.dev/support-and-community/plans-and-billing/credits/) docs—the product surface moves; **your invoice is the source of truth**.

### What this modernisation looked like as agent load

Qualitative drivers that **raise** cost:

- multi-hour conversations with large repository context;
- frontier models for architecture and numerical debugging;
- many tool rounds (tests, bulk OZ solves, pull request creation);
- re-reading long transcripts instead of starting a fresh thread per phase.

Drivers that **lower** cost:

- phased prompts (“Phase 3 only”) and new conversations per phase;
- cheaper / auto-efficient models for mechanical edits;
- strong local gates so the agent does not thrash on format and lint;
- project rules (`WARP.md`, uv-only, Australian English) so preferences are not re-explained every turn;
- bring-your-own-key if you already pay Anthropic or OpenAI directly (inference leaves Warp’s meter; harness and platform rules still apply).

### Order-of-magnitude planning numbers (not a receipt)

I am **not** pasting a private billing export into this post. For planning peer work of similar density, these bands are useful:

| Scenario | Rough Warp-related spend (USD) | Notes |
|----------|--------------------------------|--------|
| Light assist (a few hours, efficient models, local agent) | **~$5–20** | Fits inside a light Build allowance day |
| This project’s dense multi-PR arc (tooling → tests → bulk → fluids → Marimo), mostly local agents | **~$20–80** agent usage band is a sensible planning range | Often one Build month + modest reload, **or** a chunk of a higher allowance—**highly** model-dependent |
| Same arc with heavy frontier models + long single chat + cloud agents | **$100+** easily | Cloud compute and platform stack on top of inference |

**Relative value frame:** even at the upper end of the dense band, agent spend is typically **small compared with a day of senior scientific-software time**, provided you keep tests and physics review in the loop. The failure mode is not the invoice—it is **accepting plausible numerics without characterisation tests**.

Practical habit: after a heavy session, open **Settings → Billing and usage**, expand conversation totals, and write one line in the project blog note (“~X of agent usage on Phase 3”). Future-you will thank present-you.

---

## What I would repeat (and what I would not)

**Repeat**

- Phase the work; one pull-request narrative each.
- Characterise before rewriting bulk inputs.
- Keep the science package UI-free from day one of the refactor.
- Treat Marimo and Streamlit as skins.
- Kill stale xfails when artefacts catch up.

**Avoid**

- One mega-conversation for the whole thesis modernisation.
- Letting the agent “also tidy” bulk `.meta.json` timestamps into commits.
- Dual fluid selectors and dual notebook stacks without a consolidation plan.
- Claiming Fortran bit-parity you cannot prove.

---

## Closing position

Lost sources force humility. Agents force discipline.

The durable outcome is not that Warp “rewrote my PhD”. The 2022 Python rewrite already did the hard restart without Fortran. This arc made that rewrite **operable**: equations in a tested package, bulk \(c(r)\) a first-class artefact, charged short-range tails measurable, exploration UIs sharing one engine—with a **known, metered cost** for the agent assistance that accelerated the last mile.

If you maintain research or “demo that became production” codebases, the pattern is portable:

> **LLMs for generation and exploration · deterministic gates for enforcement · humans for physical and product judgement.**

---

## Related reading

- Progressive notes in-repo: [`docs/blog/`](https://github.com/Mjboothaus/PhD-Thesis/tree/main/docs/blog)
- Capstone comparison: [Before and after](https://github.com/Mjboothaus/PhD-Thesis/blob/main/docs/blog/2026-11-19-capstone-before-after.md)
- Production app: [phd-thesis-mjb.onrender.com](https://phd-thesis-mjb.onrender.com)
- Earlier Databooth themes: [deterministic guardrails with prek](https://www.databooth.com.au/posts/deterministic-llm-guardrails-with-prek/), [hyper-team software](https://www.databooth.com.au/posts/hyper-team-software/), [coding beyond the vibe](https://www.databooth.com.au/posts/coding-vibe/)
- Warp pricing: [warp.dev/pricing](https://www.warp.dev/pricing)

---

*Australian English in documentation. Code identifiers keep US spelling.*
