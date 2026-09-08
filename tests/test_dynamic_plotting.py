"""Publication artifact contracts for interaction states and uncertainty."""

from __future__ import annotations

import matplotlib
import pandas as pd

matplotlib.use("Agg")

from docklens.analytics import AnalysisContext
from docklens.dynamic_plotting import (
    build_difference_uncertainty_chart,
    build_state_population_chart,
    build_state_timeline_chart,
    build_transition_chart,
    build_uncertainty_chart,
)
from docklens.dynamic_states import interaction_state_analysis
from docklens.uncertainty import block_bootstrap_occupancy
from docklens.uncertainty import block_bootstrap_difference


def _matrix():
    return pd.DataFrame(
        [
            (1, 0),
            (1, 0),
            (0, 1),
            (0, 1),
            (1, 0),
            (1, 0),
            (0, 1),
            (0, 1),
        ],
        index=[f"f{index}" for index in range(8)],
        columns=pd.MultiIndex.from_tuples(
            [("SER70", "hbond"), ("TRP105", "pi_stacking")],
            names=("receptor_residue", "interaction_type"),
        ),
        dtype=bool,
    )


def test_state_timeline_and_transition_artifacts_expose_exact_source_rows():
    context = AnalysisContext(mode="md", time_step_ns=0.25)
    analysis = interaction_state_analysis(_matrix(), context, threshold=1.0)

    timeline = build_state_timeline_chart(analysis)
    transitions = build_transition_chart(analysis)

    assert timeline.kind == "interaction-state-timeline"
    assert len(timeline.data) == 8
    assert timeline.metadata["mode"] == "md"
    assert timeline.metadata["time_step_ns"] == 0.25
    assert timeline.figure.axes[0].get_xlabel() == "Simulation time (ns)"
    assert transitions.kind == "observed-state-transitions"
    assert transitions.metadata["interpretation"] == (
        "descriptive observed transitions; not a validated Markov model"
    )
    assert transitions.data["transition_count"].sum() == 7


def test_docking_state_population_uses_pose_family_terminology():
    analysis = interaction_state_analysis(
        _matrix(),
        AnalysisContext(mode="docking"),
        threshold=1.0,
    )

    artifact = build_state_population_chart(analysis)

    assert artifact.kind == "pose-family-population"
    assert artifact.figure.axes[0].get_title(loc="left") == (
        "Pose-family population"
    )
    assert "frequency" in artifact.figure.axes[0].get_xlabel().lower()


def test_uncertainty_chart_exports_bootstrap_parameters():
    context = AnalysisContext(mode="md")
    intervals = block_bootstrap_occupancy(
        _matrix(),
        context,
        iterations=100,
        block_size=2,
        seed=9,
    )

    artifact = build_uncertainty_chart(intervals)

    assert artifact.kind == "md-occupancy-confidence-intervals"
    assert artifact.metadata["bootstrap_method"] == (
        "circular moving-block bootstrap"
    )
    assert artifact.metadata["iterations"] == 100
    assert artifact.metadata["block_size"] == 2
    assert set(artifact.data.columns) >= {
        "occupancy_pct",
        "ci_low_pct",
        "ci_high_pct",
    }


def test_difference_uncertainty_chart_labels_independent_b_minus_a_design():
    context = AnalysisContext(mode="md")
    intervals = block_bootstrap_difference(
        _matrix(),
        _matrix().iloc[::-1],
        context,
        iterations=100,
        block_size_a=2,
        block_size_b=2,
        seed=13,
    )

    artifact = build_difference_uncertainty_chart(intervals)

    assert artifact.kind == "md-differential-occupancy-confidence-intervals"
    assert artifact.metadata["design"] == "independent B - A"
    assert artifact.metadata["multiplicity_adjustment"] == "none"
    assert "B - A" in artifact.figure.axes[0].get_xlabel()


def test_timeline_display_limit_preserves_full_md_population_and_original_time():
    matrix = pd.DataFrame(True, index=[f"f{i}" for i in range(150)],
                          columns=pd.MultiIndex.from_tuples([("SER70", "hbond")]))
    analysis = interaction_state_analysis(matrix, AnalysisContext(mode="md", time_step_ns=.5),
                                          max_training_observations=100)
    display_ids = [f"f{i}" for i in range(51, 150)]
    timeline = build_state_timeline_chart(analysis, display_observation_ids=display_ids)
    population = build_state_population_chart(analysis)
    assert timeline.metadata["total_observations"] == 150
    assert timeline.metadata["observations_displayed"] == 99
    assert timeline.data["ordinal"].tolist() == list(range(51, 150))
    assert timeline.data["time_ns"].iloc[0] == 25.5
    assert population.metadata["population_denominator"] == 150
    assert population.data["population_count"].sum() == 150
