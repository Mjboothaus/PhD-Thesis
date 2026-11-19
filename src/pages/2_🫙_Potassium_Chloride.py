"""Potassium chloride melt near an interface."""

from components.fluid_page import render_fluid_page

render_fluid_page(
    "kcl",
    experimental_message=(
        "KCl remains **experimental**: bulk short-range c(r) generation at melt "
        "density is not yet fully stable (issues #13 / #17)."
    ),
    beta_u_ylim=[0, 10000],
)
