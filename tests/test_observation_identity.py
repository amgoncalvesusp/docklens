"""Human-readable labels must never replace stable observation identifiers."""

from __future__ import annotations

from docklens.observation_identity import observation_labels
from docklens.observation_series import ObservationPoint, ObservationSeries


def test_ligand_and_file_labels_disambiguate_multipose_observations(
    multi_source_result,
):
    ligand = observation_labels(
        multi_source_result,
        mode="docking",
        label_mode="ligand",
    )
    source_file = observation_labels(
        multi_source_result,
        mode="docking",
        label_mode="file",
    )

    assert tuple(ligand) == tuple(
        summary.pose_id for summary in multi_source_result.summaries
    )
    assert tuple(ligand.values()) == (
        "LIG-A · Pose 1",
        "LIG-A · Pose 2",
        "LIG-B",
    )
    assert tuple(source_file.values()) == (
        "ligand_a.mol2 · Pose 1",
        "ligand_a.mol2 · Pose 2",
        "ligand_b.mol2",
    )


def test_index_labels_use_global_pose_order_to_avoid_cross_file_collisions(
    multi_source_result,
):
    labels = observation_labels(
        multi_source_result,
        mode="docking",
        label_mode="index",
    )

    assert tuple(labels.values()) == ("Pose 1", "Pose 2", "Pose 3")


def test_md_labels_preserve_explicit_frame_and_replica_identity(
    multi_source_result,
):
    series = ObservationSeries(
        mode="md",
        points=(
            ObservationPoint(
                "S000001:P0001:R001",
                0,
                frame_index=10,
                time_ns=1.0,
                replica_id="replica-A",
            ),
            ObservationPoint(
                "S000001:P0002:R001",
                1,
                frame_index=12,
                time_ns=1.2,
                replica_id="replica-A",
            ),
            ObservationPoint(
                "S000002:P0001:R001",
                2,
                frame_index=10,
                time_ns=1.0,
                replica_id="replica-B",
            ),
        ),
        time_step_ns=0.1,
    )

    labels = observation_labels(
        multi_source_result,
        mode="md",
        label_mode="index",
        series=series,
    )

    assert tuple(labels.values()) == (
        "Frame 10 · replica-A",
        "Frame 12 · replica-A",
        "Frame 10 · replica-B",
    )


def test_unknown_label_mode_is_rejected(multi_source_result):
    import pytest

    with pytest.raises(ValueError, match="label_mode"):
        observation_labels(
            multi_source_result,
            mode="docking",
            label_mode="unsafe",
        )
