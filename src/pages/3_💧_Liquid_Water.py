"""Liquid water (CF1 central-force model)."""

from components.fluid_page import render_fluid_page

render_fluid_page(
    "h2o",
    title_prefix="Charged fluids near an interface",
    experimental_message=(
        "CF1 water: bulk c(r) from digitised thesis Fig 3.12 (OO, OH, HH). "
        "Full bulk OZ regeneration of CF1 remains experimental."
    ),
    beta_u_ylim=[-50, 200],
)
