"""Publication-figure and sidecar export contracts."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
import pandas as pd
import pytest

matplotlib.use("Agg")

from docklens.figure_export import export_figure_bundle
from docklens.plotting import (
    build_comparison_chart,
    build_fingerprint_chart,
    build_residue_chart,
)
from docklens.results import Detail, Endpoint, Summary, make_result


def _result(prefix, observations, events):
    ligand = Endpoint("ligand", "atom", "C1", (1,), "LIG", "1")
    summaries = tuple(
        Summary(
            ligand_id=prefix,
            source_file=f"{pose_id}.pdb",
            sol=None,
            pose=index,
            docking_score=-8.0 + index,
            n_total_interactions=0,
            n_key_residue_interactions=0,
            counts={},
            pose_id=pose_id,
            source_id=f"{prefix}-{index}",
        )
        for index, pose_id in enumerate(observations, 1)
    )
    details = []
    for index, (pose_id, residue, kind) in enumerate(events, 1):
        letters = "".join(char for char in residue if char.isalpha())
        digits = "".join(char for char in residue if char.isdigit())
        receptor = Endpoint(
            "receptor", "atom", "X", (index + 1,), letters, digits
        )
        details.append(
            Detail(
                ligand_id=prefix,
                source_file=f"{pose_id}.pdb",
                interaction_type=kind,
                subtype="",
                ligand=ligand,
                receptor=receptor,
                distance_A=3.0,
                source_id=f"{prefix}-{observations.index(pose_id) + 1}",
                pose_id=pose_id,
                interaction_id=f"{prefix}-event-{index}",
                pose=observations.index(pose_id) + 1,
            )
        )
    return make_result(details=details, summaries=summaries)


@pytest.fixture
def systems():
    left = _result(
        "A",
        ("a1", "a2"),
        (
            ("a1", "SER70", "hbond"),
            ("a2", "GLU166", "saltbridge"),
        ),
    )
    right = _result(
        "B",
        ("b1", "b2", "b3"),
        (
            ("b1", "SER70", "hbond"),
            ("b2", "SER70", "hbond"),
            ("b3", "TRP105", "pi_stacking"),
        ),
    )
    return left, right


def test_residue_chart_exposes_plot_data_denominator_and_axis_label(systems):
    artifact = build_residue_chart(systems[0], mode="docking")

    assert artifact.kind == "residue-prevalence"
    assert set(artifact.data.columns) >= {
        "receptor_residue",
        "interaction_type",
        "observation_count",
        "total_observations",
    }
    assert artifact.metadata["counting_unit"] == (
        "observation × receptor residue × interaction type"
    )
    assert artifact.metadata["total_observations"] == 2
    assert artifact.figure.axes[0].get_xlabel() == "Docking frequency (% of poses)"


def test_fingerprint_and_comparison_charts_have_auditable_source_rows(systems):
    fingerprint = build_fingerprint_chart(systems[0])
    comparison = build_comparison_chart(*systems)

    assert len(fingerprint.data) == 4
    assert set(fingerprint.data["observation_id"]) == {"a1", "a2"}
    assert fingerprint.figure.axes[0].get_ylabel() == "Pose"
    assert set(comparison.data["receptor_residue"]) == {
        "GLU166",
        "SER70",
        "TRP105",
    }
    assert comparison.figure.axes[0].get_xlabel() == (
        "Δ prevalence (B − A, percentage points)"
    )


def test_figure_bundle_exports_png_svg_data_and_manifest(systems, tmp_path):
    artifact = build_residue_chart(systems[0], mode="docking")

    paths = export_figure_bundle(
        artifact,
        tmp_path / "residue-profile",
        formats=("png", "svg"),
        dpi=300,
        extra_metadata={"analysis_profile": "ds_like"},
    )

    assert {Path(path).suffix for path in paths} == {
        ".png",
        ".svg",
        ".csv",
        ".json",
    }
    assert all(Path(path).stat().st_size > 0 for path in paths)
    manifest_path = next(Path(path) for path in paths if path.endswith(".json"))
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["kind"] == "residue-prevalence"
    assert manifest["dpi"] == 300
    assert manifest["analysis_profile"] == "ds_like"
    assert manifest["total_observations"] == 2


def test_figure_bundle_validates_format_and_dpi(systems, tmp_path):
    artifact = build_residue_chart(systems[0])

    with pytest.raises(ValueError, match="formats"):
        export_figure_bundle(artifact, tmp_path / "bad", formats=("jpg",))
    with pytest.raises(ValueError, match="dpi"):
        export_figure_bundle(artifact, tmp_path / "bad", dpi=0)


def test_figure_bundle_neutralizes_formula_like_chart_labels(systems, tmp_path):
    artifact = build_residue_chart(systems[0])
    unsafe = artifact.data.copy()
    unsafe.loc[unsafe.index[0], "receptor_residue"] = "=CMD|' /C calc'!A0"
    guarded = type(artifact)(
        kind=artifact.kind,
        figure=artifact.figure,
        data=unsafe,
        metadata=artifact.metadata,
    )

    paths = export_figure_bundle(guarded, tmp_path / "safe")
    data_path = next(path for path in paths if path.endswith("_data.csv"))
    written = pd.read_csv(data_path)

    assert written.loc[0, "receptor_residue"].startswith("'=")
