"""Potassium chloride melt near an interface."""

from components.fluid_page import render_fluid_page

render_fluid_page(
    "kcl",
    experimental_message=(
        "KCl remains **experimental**: Ng/PY bulk c_short is wall-usable "
        "(near-zero large-r tail in `data/bulk/kcl.cr.dat`), but bulk g(r) / "
        "Picard residual at melt density is still imperfect."
    ),
    beta_u_ylim=[0, 10000],
)
