"""Potassium chloride melt near an interface."""

from components.fluid_page import render_fluid_page

render_fluid_page(
    "kcl",
    experimental_message=(
        "KCl remains **experimental**: Ng/PY bulk with temperature anneal "
        "(`data/bulk/kcl.cr.dat`) now reaches a small Picard residual and usable "
        "c_short tail. Bulk still uses a soft LJ+Coulomb core (not full Pauling); "
        "wall profiles vs melt data are not yet regression-tested."
    ),
    beta_u_ylim=[0, 10000],
)
