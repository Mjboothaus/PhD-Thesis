"""Discretisation parameters for the wall-fluid solver grid."""

from dataclasses import dataclass


@dataclass
class Discretisation:
    n_point: int
    n_component: int
    n_pair: int
    z_cutoff: float
    grid_size: float
    tolerance: float
    max_iteration: int


def set_num_parameters(
    n_point: int,
    z_cutoff: float,
    n_component: int,
    n_pair: int,
    tolerance: float,
    max_iteration: int,
) -> Discretisation:
    """Build a :class:`Discretisation` from UI or CLI numeric controls."""
    grid_size = z_cutoff / (n_point - 1)
    return Discretisation(
        n_point, n_component, n_pair, z_cutoff, grid_size, tolerance, max_iteration
    )
