"""Interaction heatmaps compare normalized sources or individual observations."""

from __future__ import annotations

from dataclasses import replace

import pytest

import docklens.interaction_heatmap as heatmap_module
from docklens.interaction_heatmap import build_interaction_heatmap_data
from docklens.observation_series import ObservationPoint, ObservationSeries
from docklens.plotting import build_interaction_heatmap_chart


def test_source_heatmap_uses_independent_denominators_and_zero_contact_poses(
    multi_source_result,
):
    heatmap = build_interaction_heatmap_data(
        multi_source_result,
        group_by="source",
        feature_level="residue_type",
        label_mode="ligand",
        mode="docking",
        top_n=None,
    )

    assert tuple(heatmap.matrix.index) == ("LIG-A", "LIG-B")
    assert tuple(heatmap.matrix.columns) == (("GLU166", "hbond"),)
    assert heatmap.matrix.loc["LIG-A", ("GLU166", "hbond")] == 50.0
    assert heatmap.matrix.loc["LIG-B", ("GLU166", "hbond")] == 100.0
    assert heatmap.metadata["row_denominators"] == {
        "LIG-A": 2,
        "LIG-B": 1,
    }
    assert heatmap.metadata["zero_contact_observations_included"] is True


def test_observation_heatmap_uses_human_labels_and_binary_presence(
    multi_source_result,
):
    heatmap = build_interaction_heatmap_data(
        multi_source_result,
        group_by="observation",
        feature_level="residue_type",
        label_mode="file",
        mode="docking",
        top_n=None,
    )

    assert tuple(heatmap.matrix.index) == (
        "ligand_a.mol2 · Pose 1",
        "ligand_a.mol2 · Pose 2",
        "ligand_b.mol2",
    )
    assert tuple(
        heatmap.matrix[("GLU166", "hbond")].astype(float)
    ) == (100.0, 0.0, 100.0)
    assert set(heatmap.cell_data["observation_id"]) == {
        "S000001:P0001:R001",
        "S000001:P0002:R001",
        "S000002:P0001:R001",
    }


def test_residue_only_heatmap_uses_any_interaction_not_channel_sum(
    multi_source_result,
):
    second_channel = replace(
        multi_source_result.details[0],
        interaction_id="I-A-hydrophobic",
        interaction_type="hydrophobic",
    )
    result = replace(
        multi_source_result,
        details=multi_source_result.details + (second_channel,),
    )
    heatmap = build_interaction_heatmap_data(
        result,
        group_by="source",
        feature_level="residue",
        label_mode="file",
        mode="docking",
        top_n=None,
    )

    assert tuple(heatmap.matrix.columns) == ("GLU166",)
    assert heatmap.matrix.loc["ligand_a.mol2", "GLU166"] == 50.0
    assert heatmap.metadata["feature_counting"] == (
        "binary any-interaction presence per observation and residue"
    )


def test_top_features_rank_normalized_ligand_prevalence_not_raw_pose_count(
    multi_source_result,
):
    b_hydrophobic = replace(
        multi_source_result.details[1],
        interaction_type="hydrophobic",
    )
    result = replace(
        multi_source_result,
        details=(multi_source_result.details[0], b_hydrophobic),
    )
    heatmap = build_interaction_heatmap_data(
        result,
        group_by="source",
        feature_level="residue_type",
        label_mode="ligand",
        mode="docking",
        top_n=1,
    )

    assert heatmap.metadata["top_n"] == 1
    assert heatmap.metadata["features_before_limit"] == 2
    assert heatmap.metadata["features_displayed"] == 1
    assert tuple(heatmap.matrix.columns) == (("GLU166", "hydrophobic"),)


def test_multiligand_file_is_split_without_splitting_multipose_ligand(
    multi_source_result,
):
    first_summary = replace(
        multi_source_result.summaries[0],
        ligand_id="PEP_pose_0001",
    )
    second_summary = replace(
        multi_source_result.summaries[1],
        ligand_id="PEP_pose_0002",
    )
    third_summary = replace(
        multi_source_result.summaries[2],
        source_id="S000001",
        source_file="ligand_a.mol2",
        source_path="C:/data/ligand_a.mol2",
        ligand_id="LIG-C",
    )
    third_detail = replace(
        multi_source_result.details[1],
        source_id="S000001",
        source_file="ligand_a.mol2",
        source_path="C:/data/ligand_a.mol2",
        ligand_id="LIG-C",
    )
    result = replace(
        multi_source_result,
        summaries=(first_summary, second_summary, third_summary),
        details=(multi_source_result.details[0], third_detail),
    )

    heatmap = build_interaction_heatmap_data(
        result,
        group_by="source",
        label_mode="ligand",
        top_n=None,
    )

    assert tuple(heatmap.matrix.index) == ("PEP", "LIG-C")
    assert heatmap.metadata["row_denominators"] == {
        "PEP": 2,
        "LIG-C": 1,
    }


def test_md_observation_heatmap_uses_explicit_series_order(
    multi_source_result,
):
    series = ObservationSeries(
        mode="md",
        points=(
            ObservationPoint(
                "S000002:P0001:R001", 0, frame_index=10, replica_id="B"
            ),
            ObservationPoint(
                "S000001:P0001:R001", 1, frame_index=20, replica_id="A"
            ),
            ObservationPoint(
                "S000001:P0002:R001", 2, frame_index=21, replica_id="A"
            ),
        ),
    )

    heatmap = build_interaction_heatmap_data(
        multi_source_result,
        group_by="observation",
        label_mode="index",
        mode="md",
        series=series,
        top_n=None,
    )

    assert tuple(heatmap.cell_data["observation_id"].unique()) == (
        "S000002:P0001:R001",
        "S000001:P0001:R001",
        "S000001:P0002:R001",
    )
    assert tuple(heatmap.matrix.index) == (
        "Frame 10 · B",
        "Frame 20 · A",
        "Frame 21 · A",
    )


def test_md_heatmap_rejects_incomplete_series(multi_source_result):
    series = ObservationSeries(
        mode="md",
        points=(
            ObservationPoint(
                "S000001:P0001:R001", 0, frame_index=1
            ),
        ),
    )

    with pytest.raises(ValueError, match="series"):
        build_interaction_heatmap_data(
            multi_source_result,
            group_by="observation",
            mode="md",
            series=series,
        )


def test_heatmap_chart_is_exportable_and_auditable(multi_source_result):
    artifact = build_interaction_heatmap_chart(
        multi_source_result,
        group_by="source",
        feature_level="residue_type",
        label_mode="ligand",
        mode="docking",
        top_n=40,
    )

    assert artifact.kind == "interaction-comparison-heatmap"
    assert artifact.metadata["label_mode"] == "ligand"
    assert artifact.metadata["group_by"] == "source"
    assert set(artifact.data.columns) >= {
        "row_id",
        "row_label",
        "receptor_residue",
        "interaction_type",
        "value_pct",
        "total_observations",
    }
    assert artifact.figure.axes[0].get_ylabel() == "Ligand / uploaded file"


def test_heatmap_rejects_invalid_scientific_options(multi_source_result):
    with pytest.raises(ValueError, match="group_by"):
        build_interaction_heatmap_data(
            multi_source_result,
            group_by="atom_pair",
        )
    with pytest.raises(ValueError, match="label_mode"):
        build_interaction_heatmap_data(
            multi_source_result,
            group_by="source",
            label_mode="unsafe",
        )


def test_heatmap_applies_disclosed_hard_cell_limit(
    multi_source_result, monkeypatch
):
    monkeypatch.setattr(heatmap_module, "MAX_HEATMAP_ROWS", 10)
    monkeypatch.setattr(heatmap_module, "MAX_HEATMAP_CELLS", 2)

    heatmap = build_interaction_heatmap_data(
        multi_source_result,
        group_by="observation",
        feature_level="residue_type",
        label_mode="file",
        top_n=None,
    )

    assert heatmap.matrix.size <= 2
    assert heatmap.metadata["rows_before_limit"] == 3
    assert heatmap.metadata["rows_sampled"] is True
    assert heatmap.metadata["hard_cell_limit"] == 2
    assert tuple(heatmap.cell_data["observation_id"].unique()) == (
        "S000001:P0001:R001",
        "S000002:P0001:R001",
    )


def test_adversarial_ligand_labels_remain_unique(
    multi_source_result,
):
    template = multi_source_result.summaries[0]
    summaries = tuple(
        replace(
            template,
            ligand_id=ligand_id,
            source_id=source_id,
            source_file=f"{source_id}.mol2",
            source_path=f"C:/data/{source_id}.mol2",
            pose_id=f"{source_id}:P0001:R001",
        )
        for source_id, ligand_id in (
            ("S1", "A"),
            ("S2", "A"),
            ("S3", "A · S1 [2]"),
            ("S4", "A · S1"),
        )
    )
    result = replace(
        multi_source_result,
        summaries=summaries,
        details=(),
    )

    heatmap = build_interaction_heatmap_data(
        result,
        group_by="source",
        label_mode="ligand",
        top_n=None,
    )

    assert len(heatmap.matrix.index) == len(set(heatmap.matrix.index))
