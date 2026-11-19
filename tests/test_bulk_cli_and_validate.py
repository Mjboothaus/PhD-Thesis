"""Bulk generation CLI helpers and artefact validation."""

from __future__ import annotations

import numpy as np
import pytest

from charged_fluids.bulk.artefacts import default_bulk_dir, validate_cr_artefact
from charged_fluids.bulk.cli import main as bulk_cli_main
from charged_fluids.bulk.generate import (
    c_short_tail_max,
    generate_bulk_cr,
    list_bulk_symbols,
    wall_usable_c_short,
)
from charged_fluids.parameters import clear_fluid_cache


@pytest.fixture(autouse=True)
def _clear_cache():
    clear_fluid_cache()


def test_list_bulk_symbols_includes_known_fluids():
    symbols = list_bulk_symbols()
    for s in ("lj1", "lj2", "kcl", "h2o", "2_2", "1_1"):
        assert s in symbols


def test_cli_list_exits_zero(capsys):
    code = bulk_cli_main(["--list"])
    assert code == 0
    out = capsys.readouterr().out
    assert "lj1" in out
    assert "h2o" in out


def test_cli_single_dry_run_lj1(capsys):
    code = bulk_cli_main(["lj1", "--n-points", "256", "--delta-r", "0.05", "--dry-run"])
    assert code == 0
    out = capsys.readouterr().out
    assert "lj1:" in out
    assert "tail max" in out


def test_validate_existing_lj1_artefact():
    path = default_bulk_dir() / "lj1.cr.dat"
    if not path.is_file():
        pytest.skip("lj1 bulk artefact missing")
    metrics = validate_cr_artefact(path, expected_n_pair=1, tail_tol=1e-2)
    assert metrics["n_pair"] == 1
    assert metrics["n_points"] > 10
    assert metrics["tail_abs_max"] < 1e-2


def test_validate_rejects_bad_tail(tmp_path):
    r = np.linspace(0.1, 5.0, 50)
    c = np.ones((50, 1))  # does not decay
    data = np.column_stack([r, c])
    path = tmp_path / "bad.cr.dat"
    np.savetxt(path, data)
    with pytest.raises(ValueError, match="tail"):
        validate_cr_artefact(path, tail_tol=1e-3)


def test_wall_usable_helper_on_lj_generate(tmp_path):
    result, cr, _ = generate_bulk_cr(
        "lj1", bulk_dir=tmp_path, n_points=256, delta_r=0.05, write=True
    )
    assert wall_usable_c_short(result)
    assert c_short_tail_max(result) < 1e-2
    assert cr is not None
    validate_cr_artefact(cr, expected_n_pair=1)


def test_cli_h2o_digitised(tmp_path, capsys):
    code = bulk_cli_main(["h2o", "--output-dir", str(tmp_path)])
    assert code == 0
    assert (tmp_path / "h2o.cr.dat").is_file()
    out = capsys.readouterr().out
    assert "h2o:" in out
