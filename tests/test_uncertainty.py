"""Statistical confidence contracts for saved-frame interaction occupancy."""

from __future__ import annotations

import pandas as pd
import pytest

from docklens.analytics import AnalysisContext
from docklens.observation_series import ObservationPoint, ObservationSeries
from docklens.uncertainty import (
    block_bootstrap_difference,
    block_bootstrap_occupancy,
)


def _matrix(values, feature=("SER70", "hbond")):
    return pd.DataFrame(
        [[value] for value in values],
        index=[f"f{index}" for index in range(len(values))],
        columns=pd.MultiIndex.from_tuples(
            [feature],
            names=("receptor_residue", "interaction_type"),
        ),
        dtype=bool,
    )


def test_block_bootstrap_is_deterministic_and_reports_method_parameters():
    matrix = _matrix([1, 1, 1, 1, 0, 0, 0, 0] * 2)

    first = block_bootstrap_occupancy(
        matrix,
        AnalysisContext(mode="md", time_step_ns=0.1),
        iterations=200,
        block_size=4,
        seed=77,
    )
    second = block_bootstrap_occupancy(
        matrix,
        AnalysisContext(mode="md", time_step_ns=0.1),
        iterations=200,
        block_size=4,
        seed=77,
    )
    row = first.iloc[0]

    pd.testing.assert_frame_equal(first, second)
    assert row.occupancy_pct == 50.0
    assert row.ci_low_pct <= row.occupancy_pct <= row.ci_high_pct
    assert row.block_size == 4
    assert row.iterations == 200
    assert row.method == "circular moving-block bootstrap"
    assert not bool(row.insufficient_data)


def test_bootstrap_rejects_docking_and_invalid_parameters():
    matrix = _matrix([1, 0, 1, 0, 1, 0, 1, 0])

    with pytest.raises(ValueError, match="MD"):
        block_bootstrap_occupancy(
            matrix,
            AnalysisContext(mode="docking"),
        )
    with pytest.raises(ValueError, match="iterations"):
        block_bootstrap_occupancy(
            matrix,
            AnalysisContext(mode="md"),
            iterations=20,
        )
    with pytest.raises(ValueError, match="block_size"):
        block_bootstrap_occupancy(
            matrix,
            AnalysisContext(mode="md"),
            block_size=0,
        )


def test_short_trajectory_is_explicitly_flagged_instead_of_overstated():
    frame = block_bootstrap_occupancy(
        _matrix([1, 0, 1, 0]),
        AnalysisContext(mode="md"),
        iterations=100,
        seed=5,
    )

    assert bool(frame.iloc[0].insufficient_data)
    assert "fewer than 8" in frame.iloc[0].warning


def test_independent_block_bootstrap_difference_reports_b_minus_a_interval():
    system_a = _matrix([1, 0] * 8)
    system_b = _matrix([1, 1, 1, 0] * 4)

    frame = block_bootstrap_difference(
        system_a,
        system_b,
        AnalysisContext(mode="md"),
        iterations=200,
        block_size_a=2,
        block_size_b=4,
        seed=11,
    )
    row = frame.iloc[0]

    assert row.occupancy_a_pct == 50.0
    assert row.occupancy_b_pct == 75.0
    assert row.delta_pct_points == 25.0
    assert row.ci_low_pct_points <= row.delta_pct_points
    assert row.ci_high_pct_points >= row.delta_pct_points
    assert row.delta_definition == "B - A"


def test_block_bootstrap_resamples_each_replica_without_crossing_boundaries():
    matrix = _matrix([0, 0, 0, 0, 1, 1, 1, 1])
    series = ObservationSeries(
        mode="md",
        points=tuple(
            ObservationPoint(
                f"f{index}",
                index,
                frame_index=index % 4,
                replica_id="A" if index < 4 else "B",
            )
            for index in range(8)
        ),
    )

    frame = block_bootstrap_occupancy(
        matrix,
        AnalysisContext(mode="md"),
        iterations=100,
        block_size=2,
        seed=3,
        series=series,
    )
    row = frame.iloc[0]

    assert row.occupancy_pct == 50.0
    assert row.ci_low_pct == 50.0
    assert row.ci_high_pct == 50.0
    assert row.replica_count == 2
