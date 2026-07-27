from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from docklens import batch_runner as br
from docklens.integration_manifest import ManifestError, load_manifest
from docklens.integration_result import write_integration_result


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _paired_files(tmp_path: Path, fixture_path):
    root = tmp_path / "project"
    inputs = root / "inputs"
    runs = root / "runs"
    reports = root / "reports"
    inputs.mkdir(parents=True)
    runs.mkdir()
    reports.mkdir()
    receptor = inputs / "receptor.pdb"
    receptor.write_text(
        "".join(
            line + "\n"
            for line in Path(fixture_path("minimal_complex.pdb")).read_text().splitlines()
            if line.startswith("ATOM")
        ),
        encoding="utf-8",
    )
    poses = runs / "poses.pdbqt"
    poses.write_bytes(Path(fixture_path("two_poses_sol3.pdbqt")).read_bytes())
    manifest = reports / "docklens.json"
    manifest.write_text(json.dumps({
        "schema": "docklens-launch-v1",
        "project_root": str(root),
        "run_id": "run-1",
        "receptor": {"path": "inputs/receptor.pdb", "sha256": _sha256(receptor)},
        "poses": {"path": "runs/poses.pdbqt", "sha256": _sha256(poses)},
        "result": {"path": "reports/docklens_result_run-1.json", "schema": "docklens-result-v1"},
        "options": {"hbond_preset": "plip", "key_residues": ["SER1A"]},
    }), encoding="utf-8")
    return manifest, receptor, poses


def test_vinalab_manifest_loads_confined_hashed_pair(tmp_path, fixture_path):
    manifest_path, receptor, poses = _paired_files(tmp_path, fixture_path)

    manifest = load_manifest(manifest_path)

    assert manifest.receptor_path == receptor.resolve()
    assert manifest.poses_path == poses.resolve()
    assert manifest.hbond_preset == "plip"
    assert manifest.key_residues == ("SER1A",)
    assert manifest.result_path == manifest_path.parent / "docklens_result_run-1.json"


def test_vinalab_manifest_rejects_path_escape_and_hash_mismatch(tmp_path, fixture_path):
    manifest_path, _receptor, _poses = _paired_files(tmp_path, fixture_path)
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    payload["poses"]["path"] = "../../outside.pdbqt"
    manifest_path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ManifestError) as escaped:
        load_manifest(manifest_path)
    assert escaped.value.code == "path_outside_project"

    manifest_path, _receptor, _poses = _paired_files(tmp_path / "second", fixture_path)
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    payload["receptor"]["sha256"] = "0" * 64
    manifest_path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ManifestError) as mismatched:
        load_manifest(manifest_path)
    assert mismatched.value.code == "hash_mismatch"


def test_vinalab_manifest_rejects_blank_run_id(tmp_path, fixture_path):
    manifest_path, _receptor, _poses = _paired_files(tmp_path, fixture_path)
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    payload["run_id"] = "   "
    manifest_path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ManifestError) as invalid:
        load_manifest(manifest_path)

    assert invalid.value.code == "run_id_invalid"


def test_vinalab_manifest_confines_result_to_reports(tmp_path, fixture_path):
    manifest_path, _receptor, _poses = _paired_files(tmp_path, fixture_path)
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    payload["result"]["path"] = "other.json"
    manifest_path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ManifestError) as outside_reports:
        load_manifest(manifest_path)

    assert outside_reports.value.code == "result_outside_reports"


def test_vinalab_manifest_rejects_result_collision(tmp_path, fixture_path):
    manifest_path, _receptor, _poses = _paired_files(tmp_path, fixture_path)
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    payload["result"]["path"] = "reports/docklens.json"
    manifest_path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ManifestError) as collision:
        load_manifest(manifest_path)

    assert collision.value.code == "result_path_collision"


def test_paired_runner_uses_explicit_receptor_for_each_pose(tmp_path, fixture_path):
    manifest_path, _receptor, _poses = _paired_files(tmp_path, fixture_path)
    manifest = load_manifest(manifest_path)

    result = br.run_paired(
        manifest.receptor_path,
        manifest.poses_path,
        key_residues=manifest.key_residues,
        hbond_preset=manifest.hbond_preset,
    )

    assert len(result.summaries) == 2
    assert not result.pending
    assert all(summary.resolution_method == "paired-manifest" for summary in result.summaries)
    assert all(record.receptor_atoms == 4 for record in result.input_qc)
    assert all(record.ligand_atoms == 2 for record in result.input_qc)


def test_paired_result_is_written_as_versioned_atomic_roundtrip(tmp_path, fixture_path):
    manifest_path, _receptor, _poses = _paired_files(tmp_path, fixture_path)
    manifest = load_manifest(manifest_path)
    result = br.run_paired(manifest.receptor_path, manifest.poses_path)

    output = write_integration_result(manifest, result)
    payload = json.loads(output.read_text(encoding="utf-8"))

    assert output == manifest.result_path
    assert payload["schema"] == "docklens-result-v1"
    assert payload["run_id"] == "run-1"
    assert payload["launch_manifest_sha256"] == _sha256(manifest_path)
    assert payload["totals"]["poses"] == 2
    assert payload["totals"]["interactions"] == len(result.details)
    assert len(payload["poses"]) == 2
    assert payload["inputs"]["receptor_sha256"] == manifest.receptor_sha256
