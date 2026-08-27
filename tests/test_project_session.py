"""Reproducible DockLens project save/load and provenance contracts."""

from __future__ import annotations

from dataclasses import replace
import json
import zipfile

import pytest

from docklens import project_result_codec
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
from docklens.input_plan import InputJob, InputPlan
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
        primary_ligand_group="source-1",
        comparison_ligand_group=None,
        observation_label_mode="file",
        heatmap_group_by="observation",
        heatmap_feature_level="residue",
        heatmap_top_n=80,
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
    assert loaded.primary_ligand_group == "source-1"
    assert loaded.observation_label_mode == "file"
    assert loaded.heatmap_group_by == "observation"
    assert loaded.heatmap_feature_level == "residue"
    assert loaded.heatmap_top_n == 80


def test_project_round_trip_preserves_conservative_hybrid_profile(tmp_path):
    source = tmp_path / "frames.pdb"
    source.write_text("MODEL 1\nENDMDL\n", encoding="utf-8")
    project = replace(_project(source), hbond_preset="luna_dsv")

    save_project(project, tmp_path / "hybrid.docklens")
    loaded = load_project(tmp_path / "hybrid.docklens")

    assert loaded.hbond_preset == "luna_dsv"


def test_project_round_trip_preserves_explicit_input_plan(tmp_path):
    receptor = tmp_path / "protein.pdb"
    ligand = tmp_path / "poses.pdbqt"
    receptor.write_text("ATOM\n", encoding="utf-8")
    ligand.write_text("MODEL 1\nENDMDL\n", encoding="utf-8")
    base = _project(ligand)
    plan = InputPlan(
        (
            InputJob(
                "paired",
                source_path=str(ligand),
                receptor_path=str(receptor),
                group_id="receptor_1",
            ),
        )
    )
    project = replace(
        base,
        primary=replace(
            base.primary,
            inputs=(build_project_input(receptor), build_project_input(ligand)),
            input_plan=plan,
        ),
    )

    save_project(project, tmp_path / "planned.docklens")
    loaded = load_project(tmp_path / "planned.docklens")

    assert loaded.primary.input_plan == plan
    assert "1 paired job(s)" in methods_summary(loaded)
    assert "1 external receptor(s)" in methods_summary(loaded)


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


@pytest.mark.parametrize("legacy_schema", ("1", "2", "3"))
def test_legacy_project_migrates_to_current_chart_defaults(
    tmp_path, legacy_schema
):
    source = tmp_path / "frames.pdb"
    source.write_text("MODEL 1\nENDMDL\n", encoding="utf-8")
    current = tmp_path / "current.docklens"
    save_project(_project(source), current)
    with zipfile.ZipFile(current, "r") as archive:
        members = {
            item.filename: archive.read(item.filename)
            for item in archive.infolist()
        }
    manifest = json.loads(members["manifest.json"])
    manifest["schema_version"] = legacy_schema
    manifest["project"].pop("primary_ligand_group", None)
    manifest["project"].pop("comparison_ligand_group", None)
    manifest["project"].pop("observation_label_mode", None)
    manifest["project"].pop("heatmap_group_by", None)
    manifest["project"].pop("heatmap_feature_level", None)
    manifest["project"].pop("heatmap_top_n", None)
    members["manifest.json"] = json.dumps(manifest).encode("utf-8")
    legacy = tmp_path / "legacy.docklens"
    with zipfile.ZipFile(legacy, "w") as archive:
        for name, payload in members.items():
            archive.writestr(name, payload)

    loaded = load_project(legacy)

    assert loaded.primary_ligand_group is None
    assert loaded.comparison_ligand_group is None
    assert loaded.observation_label_mode == "ligand"
    assert loaded.heatmap_group_by == "source"
    assert loaded.heatmap_feature_level == "residue_type"
    assert loaded.heatmap_top_n == 40


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
    assert "Primary chart scope" in text
    assert "poses.mol2" in text
    assert "Observation labels: Uploaded file name" in text
    assert "Heatmap: individual observations" in text


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
    payload_name = next(name for name in members if name.endswith(".ndjson"))
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


def test_schema5_writes_incremental_ndjson_chunks_and_roundtrips(
    tmp_path, monkeypatch
):
    source = tmp_path / "frames.pdb"
    source.write_text("MODEL 1\nENDMDL\n", encoding="utf-8")
    base = _project(source)
    details = tuple(
        replace(
            base.primary.result.details[0],
            interaction_id=f"interaction-{index}",
        )
        for index in range(80)
    )
    result = replace(base.primary.result, details=details)
    project = replace(base, primary=replace(base.primary, result=result))
    monkeypatch.setattr(project_result_codec, "MAX_RESULT_CHUNK_BYTES", 4096)

    destination = tmp_path / "chunked.docklens"
    save_project(project, destination)

    with zipfile.ZipFile(destination) as archive:
        names = set(archive.namelist())
        manifest = json.loads(archive.read("manifest.json"))
    assert manifest["schema_version"] == "5"
    assert not any(name.endswith("run_result.json") for name in names)
    assert len(
        [name for name in names if name.startswith("results/system-a/details/")]
    ) > 1
    assert load_project(destination) == project


def test_schema5_chunk_tampering_raises_project_integrity_error(tmp_path):
    source = tmp_path / "frames.pdb"
    source.write_text("MODEL 1\nENDMDL\n", encoding="utf-8")
    destination = tmp_path / "original.docklens"
    save_project(_project(source), destination)

    with zipfile.ZipFile(destination) as archive:
        members = {
            info.filename: archive.read(info.filename)
            for info in archive.infolist()
        }
    chunk_name = next(name for name in members if name.endswith(".ndjson"))
    members[chunk_name] = members[chunk_name] + b"\n"
    tampered = tmp_path / "tampered.docklens"
    with zipfile.ZipFile(tampered, "w") as archive:
        for name, payload in members.items():
            archive.writestr(name, payload)

    with pytest.raises(ValueError, match="integrity|hash"):
        load_project(tampered)
