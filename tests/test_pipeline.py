"""Portable smoke tests for the current parsing and detection pipeline."""

from __future__ import annotations

from docklens import batch_runner as br
from docklens import export
from docklens.entity_resolver import resolve, set_sides
from docklens.interaction_core import compute_interactions
from docklens.parser_mol2 import parse_mol2
from docklens.parser_pdb import parse_pdb
from docklens.parser_pdbqt import parse_pdbqt


def test_mol2_ccdc_fixture_is_portable(fixture_path):
    path = fixture_path("minimal_ccdc_sol16.mol2")
    poses = parse_mol2(path)

    assert len(poses) == 1
    assert poses[0].sol == 16
    resolutions = resolve(poses[0])
    assert len(resolutions) == 1
    assert resolutions[0].method == "ccdc"
    assert resolutions[0].needs_confirmation is False


def test_pdb_fixture_resolves_and_detects_hbond(fixture_path):
    poses = parse_pdb(fixture_path("minimal_complex.pdb"))

    assert len(poses) == 1
    resolution = resolve(poses[0])[0]
    assert resolution.method == "hetatm"
    set_sides(resolution)
    interactions = compute_interactions(
        resolution.receptor_atoms,
        resolution.ligand_atoms,
        resolution.waters,
    )
    assert "hbond" in {item["type"] for item in interactions}


def test_pdbqt_fixture_preserves_models_scores_and_solution(fixture_path):
    poses = parse_pdbqt(fixture_path("two_poses_sol3.pdbqt"))

    assert [pose.pose_index for pose in poses] == [0, 1]
    assert [pose.score for pose in poses] == [-7.5, -6.25]
    assert [pose.sol for pose in poses] == [3, 3]


def test_current_csv_and_xlsx_exports_still_work(fixture_path, tmp_path):
    result = br.run([fixture_path("minimal_complex.pdb")])

    csv_paths = export.export_csv(result, tmp_path / "interactions")
    assert all(path.exists() and path.stat().st_size for path in map(__import__("pathlib").Path, csv_paths))

    xlsx_path = export.export_xlsx(result, str(tmp_path / "interactions.xlsx"))
    assert __import__("pathlib").Path(xlsx_path).stat().st_size > 0
