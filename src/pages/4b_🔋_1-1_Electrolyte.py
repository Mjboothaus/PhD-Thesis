"""1-1 aqueous electrolyte (soft-core primitive model)."""

from components.fluid_page import render_fluid_page

render_fluid_page(
    "1_1",
    title_prefix="Charged fluids near an interface",
    experimental_message=(
        "1-1 electrolyte: Ng/HNC bulk c_short (experimental). "
        "Companion to the 2-2 page; parameters are a soft-core RPM starting point."
    ),
    beta_u_ylim=[0, 200],
)
