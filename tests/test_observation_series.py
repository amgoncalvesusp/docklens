"""Explicit trajectory-axis contracts for temporal analyses."""

from __future__ import annotations

import pytest

from docklens.observation_series import (
    ObservationPoint,
    ObservationSeries,
    default_observation_series,
    observation_series_from_dataframe,
)
import pandas as pd


def test_series_preserves_explicit_order_instead_of_sorting_identifiers():
    series = ObservationSeries(
        mode="md",
        points=(
            ObservationPoint("f1", 0, frame_index=1, time_ns=0.0),
            ObservationPoint("f10", 1, frame_index=2, time_ns=0.1),
            ObservationPoint("f2", 2, frame_index=3, time_ns=0.2),
        ),
        time_step_ns=0.1,
    )

    assert series.observation_ids == ("f1", "f10", "f2")


def test_series_rejects_duplicate_ids_and_regressive_frame_or_time():
    with pytest.raises(ValueError, match="unique"):
        ObservationSeries(
            mode="md",
            points=(
                ObservationPoint("f1", 0),
                ObservationPoint("f1", 1),
            ),
        )
    with pytest.raises(ValueError, match="frame_index"):
        ObservationSeries(
            mode="md",
            points=(
                ObservationPoint("f1", 0, frame_index=2),
                ObservationPoint("f2", 1, frame_index=1),
            ),
        )
    with pytest.raises(ValueError, match="time_ns"):
        ObservationSeries(
            mode="md",
            points=(
                ObservationPoint("f1", 0, time_ns=0.2),
                ObservationPoint("f2", 1, time_ns=0.1),
            ),
        )


def test_replica_boundaries_allow_restarted_indices_and_are_never_transitions():
    series = ObservationSeries(
        mode="md",
        points=(
            ObservationPoint("a1", 0, frame_index=1, replica_id="A"),
            ObservationPoint("a2", 1, frame_index=2, replica_id="A"),
            ObservationPoint("b1", 2, frame_index=1, replica_id="B"),
            ObservationPoint("b2", 3, frame_index=2, replica_id="B"),
        ),
    )

    assert series.transition_pairs() == (("a1", "a2"), ("b1", "b2"))


def test_frame_gaps_are_not_silently_counted_as_transitions():
    series = ObservationSeries(
        mode="md",
        points=(
            ObservationPoint("f1", 0, frame_index=1),
            ObservationPoint("f3", 1, frame_index=3),
            ObservationPoint("f4", 2, frame_index=4),
        ),
    )

    assert series.transition_pairs(require_consecutive_frames=True) == (
        ("f3", "f4"),
    )
    assert series.transition_pairs(require_consecutive_frames=False) == (
        ("f1", "f3"),
        ("f3", "f4"),
    )


def test_docking_series_rejects_temporal_transition_pairs():
    series = default_observation_series(("p1", "p2"), mode="docking")

    with pytest.raises(ValueError, match="MD"):
        series.transition_pairs()


def test_dataframe_map_preserves_replicas_gaps_and_optional_time():
    frame = pd.DataFrame(
        {
            "observation_id": ("a1", "a3", "b4", "b5"),
            "replica_id": ("A", "A", "B", "B"),
            "frame_index": (1, 3, 4, 5),
            "time_ns": (0.0, 0.5, 1.0, 1.25),
        }
    )

    series = observation_series_from_dataframe(
        frame,
        expected_ids=("a1", "a3", "b4", "b5"),
        default_time_step_ns=0.25,
    )

    assert series.observation_ids == ("a1", "a3", "b4", "b5")
    assert series.transition_pairs() == (("b4", "b5"),)
    assert series.points[1].time_ns == pytest.approx(0.5)


def test_dataframe_map_derives_time_and_requires_exact_observation_ids():
    frame = pd.DataFrame(
        {
            "observation_id": ("f2", "f4"),
            "replica_id": ("run-1", "run-1"),
            "frame_index": (2, 4),
        }
    )

    series = observation_series_from_dataframe(
        frame,
        expected_ids=("f2", "f4"),
        default_time_step_ns=0.1,
    )

    assert tuple(point.time_ns for point in series.points) == (0.2, 0.4)
    with pytest.raises(ValueError, match="exactly"):
        observation_series_from_dataframe(
            frame,
            expected_ids=("f2", "missing"),
            default_time_step_ns=0.1,
        )
