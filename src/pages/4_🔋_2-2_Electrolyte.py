"""2-2 aqueous electrolyte (soft-core primitive model)."""

from components.fluid_page import render_fluid_page

render_fluid_page(
    "2_2",
    title_prefix="Charged fluids near an interface",
    experimental_message=(
        "2-2 electrolyte: Ng/PY bulk c_short artefact (experimental residual). "
        "Soft-core ionic parameters are a starting point for wall smoke tests."
    ),
    beta_u_ylim=[0, 500],
)
