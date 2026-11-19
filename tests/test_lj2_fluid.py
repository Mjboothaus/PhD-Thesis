"""LJ2 equal two-component fluid readiness."""

from __future__ import annotations

import numpy as np

from charged_fluids import prepare_fluid_run, set_fluid_parameters
from charged_fluids.bulk.generate import generate_bulk_cr


def test_lj2_parameters_ready():
    fluid = set_fluid_parameters("lj2")
    assert fluid is not None
    assert fluid.n_component == 2
    assert fluid.n_pair == 3
    assert fluid.status == "ready"


def test_lj2_bulk_artefact_three_pairs(tmp_path):
    result, cr_path, _ = generate_bulk_cr(
        "lj2", bulk_dir=tmp_path, n_points=256, delta_r=0.05, write=True
    )
    assert result.converged
    assert result.c_short.shape[0] == 3
    assert cr_path is not None
    data = np.loadtxt(cr_path)
    assert data.shape[1] == 4  # r + 3 pairs


def test_prepare_lj2_run(data_dir, tmp_path):
    # Ensure bulk artefact exists in tmp then point bulk via writing there
    generate_bulk_cr("lj2", bulk_dir=tmp_path, n_points=256, delta_r=0.05, write=True)
    # prepare uses default bulk dir; copy into data/bulk if missing is ok — use real bulk
    import shutil

    from charged_fluids.bulk.artefacts import default_bulk_dir

    dest = default_bulk_dir()
    dest.mkdir(parents=True, exist_ok=True)
    for name in ("lj2.cr.dat", "lj2.meta.json"):
        src = tmp_path / name
        if src.is_file():
            shutil.copy(src, dest / name)

    run = prepare_fluid_run("lj2", z_cutoff=8.0, n_point=64, data_dir=data_dir)
    assert run.model.c_short.shape == (64, 3)
    assert run.fluid.n_pair == 3
