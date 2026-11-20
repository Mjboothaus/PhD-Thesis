"""CF1 water and 2-2 electrolyte readiness tests."""

from __future__ import annotations

import numpy as np
import pytest

from charged_fluids import prepare_fluid_run, set_fluid_parameters
from charged_fluids.bulk.build_digitised_h2o import build_h2o_cr
from charged_fluids.bulk.generate import generate_bulk_cr
from charged_fluids.potentials.cf1 import CF1_Q_H, CF1_Q_O, calc_u_cf1, cf1_pair_potential_kcal


def test_cf1_charges_match_literature():
    assert pytest.approx(0.32983) == CF1_Q_H
    assert pytest.approx(-2.0 * CF1_Q_H) == CF1_Q_O


def test_cf1_oh_minimum_near_bond_length():
    r = np.linspace(0.5, 1.5, 400)
    u = cf1_pair_potential_kcal(r)[:, 1]  # OH
    imin = int(np.argmin(u))
    assert r[imin] == pytest.approx(0.9584, rel=0.05)


def test_h2o_toml_is_cf1():
    fluid = set_fluid_parameters("h2o")
    assert fluid is not None
    assert fluid.status == "experimental"
    assert fluid.component == ["O", "H"]
    assert fluid.n_pair == 3
    assert fluid.potential.get("type") == "cf1"
    np.testing.assert_allclose(fluid.valence, [CF1_Q_O, CF1_Q_H], rtol=1e-5)


def test_calc_u_cf1_shape_and_finite():
    r = np.linspace(0.0, 6.0, 61)
    u = calc_u_cf1(len(r), r)
    assert u.shape == (61, 3)
    assert np.isfinite(u[1:]).all()


def test_build_digitised_h2o_and_prepare_from_tmp(tmp_path, data_dir, monkeypatch):
    """Digitised Fig 3.12 path stays available; do not clobber data/bulk/h2o.cr.dat."""
    from charged_fluids.bulk import artefacts as bulk_artefacts

    cr, meta = build_h2o_cr(bulk_dir=tmp_path, n_points=256, r_max=8.0)
    assert cr.is_file() and meta.is_file()
    data = np.loadtxt(cr)
    assert data.shape[1] == 4
    # Point bulk resolution at tmp so prepare_fluid_run uses digitised data only here
    monkeypatch.setattr(bulk_artefacts, "default_bulk_dir", lambda: tmp_path)
    run = prepare_fluid_run("h2o", z_cutoff=8.0, n_point=64, data_dir=data_dir)
    assert run.model.c_short.shape == (64, 3)
    assert run.u_r.shape == (64, 3)
    assert np.isfinite(run.beta_u[1:]).all()


@pytest.mark.slow
def test_cf1_bulk_oz_anneal_converges(tmp_path):
    """Live CF1 Ng/HNC temperature anneal → finite c_short with small tail."""
    result, cr_path, _ = generate_bulk_cr(
        "h2o",
        bulk_dir=tmp_path,
        n_points=512,
        delta_r=0.05,
        write=True,
    )
    assert cr_path is not None and cr_path.is_file()
    assert result.converged
    assert np.isfinite(result.c_short).all()
    assert np.isfinite(result.g).all()
    tail = float(np.max(np.abs(result.c_short[:, -max(10, result.r.size // 20) :])))
    assert tail < 1e-2, f"CF1 |c_short| tail max={tail}"
    # Intramolecular OH peak near 1 Å; bare HNC under-associates (coord < 2)
    g_oh = result.g[1]
    r_peak = float(result.r[int(np.argmax(g_oh))])
    assert 0.85 <= r_peak <= 1.15
    assert float(result.meta.get("coord_OH_1.3A", 0.0)) > 0.5


def test_2_2_parameters_and_bulk(tmp_path):
    fluid = set_fluid_parameters("2_2")
    assert fluid is not None
    assert fluid.status == "experimental"
    assert fluid.potential.get("type") == "ionic"
    np.testing.assert_allclose(fluid.valence, [2.0, -2.0])
    result, cr_path, _ = generate_bulk_cr(
        "2_2", bulk_dir=tmp_path, n_points=256, delta_r=0.05, write=True
    )
    assert cr_path is not None
    assert result.c_short.shape[0] == 3
    tail = float(abs(result.c_short[:, -max(5, result.r.size // 20) :]).max())
    assert tail < 1e-2


def test_prepare_2_2_with_artefact(data_dir, tmp_path):
    import shutil

    from charged_fluids.bulk.artefacts import default_bulk_dir

    generate_bulk_cr("2_2", bulk_dir=tmp_path, n_points=256, delta_r=0.05, write=True)
    dest = default_bulk_dir()
    dest.mkdir(parents=True, exist_ok=True)
    for name in ("2_2.cr.dat", "2_2.meta.json"):
        src = tmp_path / name
        if src.is_file():
            shutil.copy(src, dest / name)
    run = prepare_fluid_run("2_2", z_cutoff=8.0, n_point=64, data_dir=data_dir)
    assert run.model.c_short.shape == (64, 3)
    assert run.kappa > 0.0
