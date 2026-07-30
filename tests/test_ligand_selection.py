"""Ligand/uploaded-file chart-scope contracts."""

from __future__ import annotations

import pytest

from docklens.ligand_selection import (
    default_md_series_for_result,
    ligand_groups,
    subset_observation_series,
    subset_run_result,
)
from docklens.observation_series import ObservationPoint, ObservationSeries


def test_catalog_groups_multipose_files_without_conflating_ligand_ids(
    multi_source_result,
):
    groups = ligand_groups(multi_source_result)

    assert tuple(group.key for group in groups) == ("S000001", "S000002")
    assert groups[0].observation_count == 2
    assert groups[0].ligand_ids == ("LIG-A",)
    assert groups[1].observation_count == 1
    assert "ligand_a.mol2" in groups[0].label


def test_subset_preserves_empty_poses_denominator_and_original_result(
    multi_source_result,
):
    selected = subset_run_result(multi_source_result, "S000001")

    assert tuple(item.pose_id for item in selected.summaries) == (
        "S000001:P0001:R001",
        "S000001:P0002:R001",
    )
    assert tuple(item.pose_id for item in selected.details) == (
        "S000001:P0001:R001",
    )
    assert len(selected.input_qc) == 1
    assert len(multi_source_result.summaries) == 3


def test_subset_rejects_unknown_group_instead_of_selecting_silently(
    multi_source_result,
):
    with pytest.raises(ValueError, match="unknown"):
        subset_run_result(multi_source_result, "missing")


def test_explicit_series_subset_preserves_replica_frames_time_and_gaps():
    series = ObservationSeries(
        mode="md",
        points=(
            ObservationPoint("a1", 0, 1, 0.1, "A"),
            ObservationPoint("a3", 1, 3, 0.3, "A"),
            ObservationPoint("b1", 2, 1, 0.1, "B"),
        ),
        time_step_ns=0.1,
    )

    selected = subset_observation_series(series, ("a1", "a3"))

    assert selected.observation_ids == ("a1", "a3")
    assert tuple(point.frame_index for point in selected.points) == (1, 3)
    assert tuple(point.time_ns for point in selected.points) == (0.1, 0.3)
    assert selected.transition_pairs() == ()


def test_automatic_md_axis_never_crosses_uploaded_file_boundaries(
    multi_source_result,
):
    series = default_md_series_for_result(
        multi_source_result,
        time_step_ns=0.25,
    )

    assert series.transition_pairs() == (
        ("S000001:P0001:R001", "S000001:P0002:R001"),
    )
    assert len({point.replica_id for point in series.points}) == 2
