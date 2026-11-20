"""Liquid water (CF1 central-force model)."""

from components.fluid_page import render_fluid_page

render_fluid_page(
    "h2o",
    title_prefix="Charged fluids near an interface",
    experimental_message=(
        "CF1 water: live Ng/HNC bulk OZ with temperature anneal "
        "(`data/bulk/h2o.cr.dat`). Bare HNC under-associates vs Duh–Haymet bridge; "
        "wall residual vs thesis Fig 3.2 remains experimental. "
        "Digitised Fig 3.12 bulk is opt-in via `[h2o.bulk] source`."
    ),
    beta_u_ylim=[-50, 200],
)
