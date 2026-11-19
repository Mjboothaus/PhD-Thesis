"""Streamlit sidebar controls for fluid pages."""

from __future__ import annotations

import streamlit as st


def create_sidebar(fluid):
    st.sidebar.markdown("## Inputs")
    st.sidebar.markdown(f"**Fluid:** {fluid.name} (`{fluid.symbol}`)")
    st.sidebar.markdown(f"Temperature (K): {fluid.temperature}")
    st.sidebar.markdown(f"Concentration[s] (mol/L): {fluid.concentration}")
    st.sidebar.caption(f"Status: `{getattr(fluid, 'status', 'unknown')}`")

    if fluid.symbol not in {"lj1", "lj2"}:
        psi_0 = st.sidebar.number_input(
            label="ψ₀ (mV)",
            value=0.0,
            min_value=-1000.0,
            max_value=1000.0,
            step=10.0,
        )
    else:
        psi_0 = 0.0

    st.sidebar.markdown("### Numerical parameters")
    z_cutoff = st.sidebar.number_input(
        label="z cutoff (Å)", value=20.0, min_value=10.0, max_value=50.0, step=1.0
    )
    n_point = st.sidebar.number_input(
        label="n_point", value=201, min_value=51, max_value=4001, step=50
    )
    grid_size = float(z_cutoff) / (int(n_point) - 1)
    st.sidebar.markdown(f"**Grid size (Å):** {grid_size:.4g}")

    st.sidebar.markdown("### Optimisation parameters")
    max_iteration = st.sidebar.number_input(
        "Maximum iterations", min_value=1, max_value=1000, value=30
    )
    tolerance = st.sidebar.number_input(
        "Convergence tolerance (× 1e-9)",
        min_value=0.001,
        max_value=1000.0,
        value=1.0,
    )
    tolerance = float(tolerance) * 1.0e-9

    return float(z_cutoff), int(n_point), float(psi_0), tolerance, int(max_iteration)
