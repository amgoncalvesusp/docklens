"""Contracts for interaction states, pose families and MD transitions."""

from __future__ import annotations

import pandas as pd
import pytest

from docklens.analytics import AnalysisContext
from docklens.dynamic_states import (
    interaction_state_analysis,
    state_assignment_frame,
    state_summary_frame,
    state_transition_frame,
)


def _matrix(rows):
    columns = pd.MultiIndex.from_tuples(
        [
            ("SER70", "hbond"),
            ("GLU166", "saltbridge"),
            ("TRP105", "pi_stacking"),
        ],
        names=("receptor_residue", "interaction_type"),
    )
    return pd.DataFrame(
        [values for _name, values in rows],
        index=[name for name, _values in rows],
        columns=columns,
        dtype=bool,
    )


def test_states_are_deterministic_and_keep_representative_and_consensus_features():
    matrix = _matrix(
        [
            ("f1", (1, 1, 0)),
            ("f2", (1, 1, 0)),
            ("f3", (0, 0, 1)),
            ("f4", (0, 0, 1)),
        ]
    )

    analysis = interaction_state_analysis(
        matrix,
        AnalysisContext(mode="md", time_step_ns=0.25),
        threshold=0.8,
    )

    assert tuple(state.state_id for state in analysis.states) == ("S1", "S2")
    assert analysis.states[0].members == ("f1", "f2")
    assert analysis.states[0].representative == "f1"
    assert analysis.states[0].characteristic_features == (
        ("GLU166", "saltbridge"),
        ("SER70", "hbond"),
    )
    assert analysis.states[0].population_pct == 50.0


def test_md_state_assignments_preserve_frame_order_time_and_dwell_statistics():
    matrix = _matrix(
        [
            ("f1", (1, 0, 0)),
            ("f2", (1, 0, 0)),
            ("f3", (0, 0, 1)),
            ("f4", (0, 0, 1)),
            ("f5", (1, 0, 0)),
        ]
    )

    analysis = interaction_state_analysis(
        matrix,
        AnalysisContext(mode="md", time_step_ns=0.5),
        threshold=1.0,
    )
    assignments = state_assignment_frame(analysis)
    summaries = state_summary_frame(analysis).set_index("state_id")

    assert assignments["observation_id"].tolist() == [
        "f1",
        "f2",
        "f3",
        "f4",
        "f5",
    ]
    assert assignments["time_ns"].tolist() == [0.0, 0.5, 1.0, 1.5, 2.0]
    assert summaries.loc["S1", "episode_count"] == 2
    assert summaries.loc["S1", "mean_dwell_observations"] == 1.5
    assert summaries.loc["S1", "longest_dwell_ns"] == 1.0


def test_md_transition_probabilities_include_self_transitions_and_normalize_rows():
    matrix = _matrix(
        [
            ("f1", (1, 0, 0)),
            ("f2", (1, 0, 0)),
            ("f3", (0, 0, 1)),
            ("f4", (0, 0, 1)),
            ("f5", (1, 0, 0)),
        ]
    )
    analysis = interaction_state_analysis(
        matrix,
        AnalysisContext(mode="md"),
        threshold=1.0,
    )

    transitions = state_transition_frame(analysis)
    s1 = transitions[transitions["from_state"] == "S1"]

    assert dict(zip(s1["to_state"], s1["transition_count"])) == {
        "S1": 1,
        "S2": 1,
    }
    assert s1["transition_probability_pct"].sum() == pytest.approx(100.0)


def test_docking_families_never_claim_temporal_transition_statistics():
    matrix = _matrix([("p1", (1, 0, 0)), ("p2", (0, 0, 1))])
    analysis = interaction_state_analysis(
        matrix,
        AnalysisContext(mode="docking"),
        threshold=1.0,
    )

    assert all(state.mean_dwell_observations is None for state in analysis.states)
    with pytest.raises(ValueError, match="MD"):
        state_transition_frame(analysis)


def test_large_analysis_discloses_training_sample_and_marks_unrepresented_outliers():
    rows = [(f"f{index}", (1, 0, 0)) for index in range(1, 11)]
    rows[5] = ("f6", (0, 0, 1))
    analysis = interaction_state_analysis(
        _matrix(rows),
        AnalysisContext(mode="md"),
        threshold=1.0,
        max_training_observations=4,
    )

    assignments = state_assignment_frame(analysis)

    assert analysis.sampled
    assert analysis.training_observations == 4
    assert analysis.total_observations == 10
    assert "OUTLIER" in set(assignments["state_id"])
    summaries = state_summary_frame(analysis)
    assert summaries["population_count"].sum() == 10
    assert "OUTLIER" in set(summaries["state_id"])


def test_family_clustering_does_not_merge_dissimilar_endpoints_by_chaining():
    matrix = _matrix(
        [
            ("a", (1, 0, 0)),
            ("b", (1, 1, 0)),
            ("c", (0, 1, 0)),
        ]
    )

    analysis = interaction_state_analysis(
        matrix,
        AnalysisContext(mode="docking"),
        threshold=0.5,
    )

    assert len(analysis.states) == 2
    assert analysis.states[0].members == ("a", "b")
    assert analysis.method == "complete-link threshold clustering"


def test_state_ids_and_medoids_do_not_depend_on_dataframe_row_order():
    matrix = _matrix(
        [
            ("f1", (1, 1, 0)),
            ("f2", (1, 1, 0)),
            ("f3", (0, 0, 1)),
        ]
    )

    forward = interaction_state_analysis(
        matrix, AnalysisContext(mode="docking"), threshold=1.0
    )
    reverse = interaction_state_analysis(
        matrix.iloc[::-1], AnalysisContext(mode="docking"), threshold=1.0
    )

    forward_map = dict(
        zip(
            state_assignment_frame(forward)["observation_id"],
            state_assignment_frame(forward)["state_id"],
        )
    )
    reverse_map = dict(
        zip(
            state_assignment_frame(reverse)["observation_id"],
            state_assignment_frame(reverse)["state_id"],
        )
    )
    assert forward_map == reverse_map
    assert tuple(state.representative for state in forward.states) == (
        "f1",
        "f3",
    )
