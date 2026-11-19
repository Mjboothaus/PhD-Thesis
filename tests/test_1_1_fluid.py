"""1-1 aqueous electrolyte bulk + setup."""

from __future__ import annotations

import numpy as np

from charged_fluids import clear_fluid_cache, prepare_fluid_run, set_fluid_parameters
from charged_fluids.bulk.generate import generate_bulk_cr


def setup_module():
    clear_fluid_cache()


def test_1_1_toml():
    fluid = set_fluid_parameters("1_1")
    assert fluid is not None
    assert fluid.status == "experimental"
    assert fluid.n_component == 2
    np.testing.assert_allclose(fluid.valence, [1.0, -1.0])
    assert fluid.potential["type"] == "ionic"


def test_1_1_bulk_generate(tmp_path):
    result, cr, _ = generate_bulk_cr("1_1", bulk_dir=tmp_path, n_points=256, write=True)
    assert cr is not None
    assert result.c_short.shape[0] == 3
    tail = float(abs(result.c_short[:, -max(5, result.r.size // 20) :]).max())
    assert tail < 1e-2


def test_prepare_1_1(data_dir, tmp_path):
    import shutil

    from charged_fluids.bulk.artefacts import default_bulk_dir

    generate_bulk_cr("1_1", bulk_dir=tmp_path, n_points=256, write=True)
    dest = default_bulk_dir()
    dest.mkdir(parents=True, exist_ok=True)
    for name in ("1_1.cr.dat", "1_1.meta.json"):
        src = tmp_path / name
        if src.is_file():
            shutil.copy(src, dest / name)
    run = prepare_fluid_run("1_1", z_cutoff=8.0, n_point=64, data_dir=data_dir)
    assert run.model.c_short.shape == (64, 3)
    assert run.kappa > 0.0
