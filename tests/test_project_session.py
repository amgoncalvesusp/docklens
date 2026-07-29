"""Reproducible DockLens project save/load and provenance contracts."""

from __future__ import annotations

from dataclasses import replace
import json
import zipfile

import pytest

from docklens.project_session import (
    ProjectDataset,
    ProjectInput,
    ProjectState,
    build_project_input,
    load_project,
    methods_summary,
    save_project,
    validate_project_inputs,
)
from docklens.observation_series import ObservationPoint, ObservationSeries
from docklens.results import Detail, Endpoint, Summary, make_result


def _project(source):
    ligand = Endpoint("ligand", "atom", "C1", (1,), "LIG", "1")
    receptor = Endpoint("receptor", "atom", "OE1", (2,), "GLU", "166")
    result = make_result(
        summaries=(
            Summary(
                ligand_id="LIG",
                source_file=source.name,
                sol=None,
                pose=1,
                docking_score=-8.5,
                n_total_interactions=1,
                n_key_residue_interactions=1,
                counts={"hbond": 1},
                source_id="source-1",
                pose_id="pose-1",
                source_path=str(source),
            ),
        ),
        details=(
            Detail(
                ligand_id="LIG",
                source_file=source.name,
                interaction_type="hbond",
                subtype="Conventional Hydrogen Bond",
                ligand=ligand,
                receptor=receptor,
                distance_A=2.8,
                source_id="source-1",
                pose_id="pose-1",
                interaction_id="interaction-1",
                pose=1,
                docking_score=-8.5,
                source_path=str(source),
                is_key_residue=True,
            ),
        ),
        key_residues=("GLU166",),
        receptor_residues=("GLU166",),
    )
    return ProjectState(
        app_version="1.0.0",
        analysis_profile="ds_like",
        hbond_preset="dsv",
        key_residues=("GLU166", "SER70"),
        selected_types=("hbond", "saltbridge"),
        active_workspace="states",
        selected_residue="GLU166",
        state_threshold=0.7,
        bootstrap_iterations=500,
        primary=ProjectDataset(
            label="System A",
            mode="md",
            time_step_ns=0.25,
            inputs=(build_project_input(source),),
            result=result,
            observation_series=ObservationSeries(
                mode="md",
                points=(
                    ObservationPoint(
                        "pose-1",
                        0,
                        frame_index=5,
                        time_ns=1.25,
                        replica_id="replica-A",
                    ),
                ),
                time_step_ns=0.25,
            ),
        ),
        bootstrap_block_size=4,
        bootstrap_seed=77,
        confidence_level=0.95,
    )


def test_project_round_trip_preserves_settings_hashes_and_methods(tmp_path):
    source = tmp_path / "frames.pdb"
    source.write_text("MODEL 1\nENDMDL\n", encoding="utf-8")
    project = _project(source)

    outputs = save_project(project, tmp_path / "analysis.docklens")
    loaded = load_project(tmp_path / "analysis.docklens")

    assert loaded == project
    assert len(outputs) == 2
    assert outputs[1].name == "analysis_methods.txt"
    assert "Discovery Studio-like" in outputs[1].read_text(encoding="utf-8")
    assert validate_project_inputs(loaded) == ()
    assert loaded.primary.result == project.primary.result
    assert loaded.primary.observation_series == project.primary.observation_series


def test_project_detects_missing_or_changed_sources(tmp_path):
    source = tmp_path / "frames.pdb"
    source.write_text("original", encoding="utf-8")
    project = _project(source)

    source.write_text("changed", encoding="utf-8")
    messages = validate_project_inputs(project)

    assert len(messages) == 1
    assert "changed" in messages[0].lower()


def test_project_loader_rejects_oversized_unknown_or_invalid_documents(tmp_path):
    oversized = tmp_path / "oversized.docklens"
    oversized.write_bytes(b" " * (5 * 1024 * 1024 + 1))
    with pytest.raises(ValueError, match="too large"):
        load_project(oversized)

    unknown = tmp_path / "unknown.docklens"
    with zipfile.ZipFile(unknown, "w") as archive:
        archive.writestr(
            "manifest.json",
            json.dumps({"schema_version": "99", "primary": {}}),
        )
    with pytest.raises(ValueError, match="schema"):
        load_project(unknown)

    malformed = tmp_path / "malformed.docklens"
    with zipfile.ZipFile(malformed, "w") as archive:
        archive.writestr("manifest.json", "[]")
    with pytest.raises(ValueError, match="object"):
        load_project(malformed)


def test_methods_summary_discloses_counting_state_and_bootstrap_assumptions(tmp_path):
    source = tmp_path / "poses.mol2"
    source.write_text("@<TRIPOS>MOLECULE\n", encoding="utf-8")

    text = methods_summary(_project(source))

    assert "one presence per observation" in text
    assert "Jaccard/Tanimoto" in text
    assert "threshold 0.700" in text
    assert "circular moving-block bootstrap" in text
    assert "4 saved frames" in text
    assert "seed 77" in text
    assert "95.0% confidence" in text
    assert "1 replica" in text


def test_comparison_series_round_trip_and_methods_are_disclosed(tmp_path):
    source = tmp_path / "frames.pdb"
    source.write_text("MODEL 1\nENDMDL\n", encoding="utf-8")
    base = _project(source)
    comparison_series = ObservationSeries(
        mode="md",
        points=(
            ObservationPoint(
                "pose-1",
                0,
                frame_index=8,
                time_ns=2.0,
                replica_id="run-B",
            ),
        ),
        time_step_ns=0.25,
    )
    project = replace(
        base,
        comparison=ProjectDataset(
            label="System B",
            mode="md",
            time_step_ns=0.25,
            inputs=base.primary.inputs,
            result=base.primary.result,
            observation_series=comparison_series,
        ),
    )

    save_project(project, tmp_path / "comparison.docklens")
    loaded = load_project(tmp_path / "comparison.docklens")

    assert loaded.comparison is not None
    assert loaded.comparison.observation_series == comparison_series
    assert "System B trajectory mapping" in methods_summary(loaded)


def test_project_rejects_internal_hash_tampering_and_path_traversal(tmp_path):
    source = tmp_path / "frames.pdb"
    source.write_text("MODEL 1\nENDMDL\n", encoding="utf-8")
    original = tmp_path / "original.docklens"
    save_project(_project(source), original)

    with zipfile.ZipFile(original, "r") as archive:
        members = {
            info.filename: archive.read(info.filename)
            for info in archive.infolist()
        }
    payload_name = next(
        name for name in members if name.endswith("run_result.json")
    )
    members[payload_name] += b" "
    tampered = tmp_path / "tampered.docklens"
    with zipfile.ZipFile(tampered, "w") as archive:
        for name, payload in members.items():
            archive.writestr(name, payload)

    with pytest.raises(ValueError, match="integrity|hash"):
        load_project(tampered)

    traversal = tmp_path / "traversal.docklens"
    with zipfile.ZipFile(traversal, "w") as archive:
        archive.writestr("../outside.txt", "unsafe")
        archive.writestr("manifest.json", "{}")
    with pytest.raises(ValueError, match="unsafe"):
        load_project(traversal)


def test_project_input_rejects_network_and_uri_paths_before_validation():
    for unsafe in (
        r"\\server\share\frames.pdb",
        "//server/share/frames.pdb",
        "file://server/share/frames.pdb",
        "https://example.test/frames.pdb",
    ):
        with pytest.raises(ValueError, match="local"):
            ProjectInput(
                path=unsafe,
                sha256="0" * 64,
                size_bytes=1,
                modified_ns=1,
            )
