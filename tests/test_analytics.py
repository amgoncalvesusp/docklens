"""Scientific contracts for DockLens 1.0 analytical views."""

from __future__ import annotations

import pandas as pd
import pytest

from docklens.analytics import (
    AnalysisContext,
    differential_prevalence,
    docking_md_retention,
    episode_statistics,
    fingerprint_clusters,
    fingerprint_matrix,
    fingerprint_similarity,
    residue_type_prevalence,
)
from docklens.observation_series import ObservationPoint, ObservationSeries
from docklens.results import Detail, Endpoint, Summary, make_result


def _endpoint(side, residue, atom):
    letters = "".join(char for char in residue if char.isalpha())
    digits = "".join(char for char in residue if char.isdigit())
    return Endpoint(
        side=side,
        kind="atom",
        atom_name=atom,
        atom_serials=(1,),
        resname="LIG" if side == "ligand" else letters,
        resseq="1" if side == "ligand" else digits,
        chain="",
    )


def _result(observations, events, *, scores=None, key_residues=()):
    scores = scores or {}
    summaries = tuple(
        Summary(
            ligand_id="LIG",
            source_file=f"{pose_id}.pdb",
            sol=None,
            pose=index,
            docking_score=scores.get(pose_id),
            n_total_interactions=0,
            n_key_residue_interactions=0,
            counts={},
            source_id=f"source-{index}",
            pose_id=pose_id,
        )
        for index, pose_id in enumerate(observations, 1)
    )
    details = tuple(
        Detail(
            ligand_id="LIG",
            source_file=f"{pose_id}.pdb",
            interaction_type=kind,
            subtype="",
            ligand=_endpoint("ligand", "LIG1", ligand_atom),
            receptor=_endpoint("receptor", residue, receptor_atom),
            distance_A=distance,
            source_id=f"source-{observations.index(pose_id) + 1}",
            pose_id=pose_id,
            interaction_id=f"event-{index}",
            pose=observations.index(pose_id) + 1,
            is_key_residue=residue in key_residues,
        )
        for index, (
            pose_id,
            residue,
            kind,
            ligand_atom,
            receptor_atom,
            distance,
        ) in enumerate(events, 1)
    )
    return make_result(
        summaries=summaries,
        details=details,
        key_residues=key_residues,
        receptor_residues={event[1] for event in events},
    )


def test_residue_prevalence_consolidates_atom_pairs_per_observation():
    result = _result(
        ["p1", "p2", "p3"],
        [
            ("p1", "GLU166", "hbond", "H1", "OE1", 2.8),
            ("p1", "GLU166", "hbond", "H2", "OE2", 2.9),
            ("p2", "GLU166", "hbond", "H3", "OE1", 3.0),
            ("p2", "SER70", "saltbridge", "N1", "OG", 3.8),
        ],
        key_residues=("GLU166",),
    )

    frame = residue_type_prevalence(result)
    glu = frame.query(
        "receptor_residue == 'GLU166' and interaction_type == 'hbond'"
    ).iloc[0]

    assert glu.observation_count == 2
    assert glu.total_observations == 3
    assert glu.prevalence_pct == pytest.approx(66.6667, rel=1e-4)
    assert bool(glu.is_key_residue)


def test_fingerprint_matrix_similarity_and_threshold_clusters_are_deterministic():
    result = _result(
        ["p1", "p2", "p3"],
        [
            ("p1", "SER70", "hbond", "H1", "OG", 2.8),
            ("p1", "GLU166", "saltbridge", "N1", "OE1", 3.5),
            ("p2", "SER70", "hbond", "H2", "OG", 2.9),
            ("p2", "GLU166", "saltbridge", "N2", "OE2", 3.6),
            ("p3", "TRP105", "pi_stacking", "C1", "CZ2", 4.8),
        ],
    )

    matrix = fingerprint_matrix(result)
    similarity = fingerprint_similarity(matrix)
    clusters = fingerprint_clusters(matrix, threshold=0.75)

    assert matrix.shape == (3, 3)
    assert matrix.loc["p1", ("SER70", "hbond")]
    assert similarity.loc["p1", "p2"] == 1.0
    assert similarity.loc["p1", "p3"] == 0.0
    assert tuple(cluster.members for cluster in clusters) == (("p1", "p2"), ("p3",))
    assert clusters[0].medoid == "p1"


def test_empty_fingerprints_are_identical_only_to_each_other():
    matrix = pd.DataFrame(
        [[False], [False], [True]],
        index=["empty-a", "empty-b", "hit"],
        columns=pd.MultiIndex.from_tuples([("SER70", "hbond")]),
    )

    similarity = fingerprint_similarity(matrix)

    assert similarity.loc["empty-a", "empty-b"] == 1.0
    assert similarity.loc["empty-a", "hit"] == 0.0


def test_md_episode_statistics_report_occupancy_runs_and_time():
    result = _result(
        ["f1", "f2", "f3", "f4", "f5", "f6"],
        [
            ("f1", "GLU166", "hbond", "H1", "OE1", 2.8),
            ("f2", "GLU166", "hbond", "H2", "OE1", 2.9),
            ("f4", "GLU166", "hbond", "H3", "OE2", 3.1),
            ("f5", "GLU166", "hbond", "H4", "OE2", 3.0),
            ("f6", "GLU166", "hbond", "H5", "OE2", 2.7),
        ],
    )

    metric = episode_statistics(
        result, AnalysisContext(mode="md", time_step_ns=0.25)
    )[0]

    assert metric.observation_count == 5
    assert metric.total_observations == 6
    assert metric.occupancy_pct == pytest.approx(83.3333, rel=1e-4)
    assert metric.episode_count == 2
    assert metric.longest_episode_observations == 3
    assert metric.mean_episode_observations == 2.5
    assert metric.longest_episode_ns == 0.75
    assert metric.mean_distance_A == pytest.approx(2.9)


def test_md_episodes_do_not_cross_replica_boundaries_or_frame_gaps():
    result = _result(
        ["a1", "a2", "b1", "b3"],
        [
            ("a2", "GLU166", "hbond", "H1", "OE1", 2.8),
            ("b1", "GLU166", "hbond", "H2", "OE1", 2.9),
            ("b3", "GLU166", "hbond", "H3", "OE2", 3.0),
        ],
    )
    series = ObservationSeries(
        mode="md",
        points=(
            ObservationPoint("a1", 0, frame_index=1, replica_id="A"),
            ObservationPoint("a2", 1, frame_index=2, replica_id="A"),
            ObservationPoint("b1", 2, frame_index=1, replica_id="B"),
            ObservationPoint("b3", 3, frame_index=3, replica_id="B"),
        ),
        time_step_ns=0.25,
    )

    metric = episode_statistics(
        result,
        AnalysisContext(mode="md", time_step_ns=0.25),
        series=series,
    )[0]

    assert metric.episode_count == 3
    assert metric.longest_episode_observations == 1


def test_differential_prevalence_uses_independent_denominators():
    system_a = _result(
        ["a1", "a2"],
        [("a1", "SER70", "hbond", "H1", "OG", 2.8)],
    )
    system_b = _result(
        ["b1", "b2", "b3", "b4"],
        [
            ("b1", "SER70", "hbond", "H1", "OG", 2.8),
            ("b2", "SER70", "hbond", "H2", "OG", 2.9),
            ("b3", "SER70", "hbond", "H3", "OG", 3.0),
        ],
    )

    row = differential_prevalence(system_a, system_b).iloc[0]

    assert row.prevalence_a_pct == 50.0
    assert row.prevalence_b_pct == 75.0
    assert row.delta_pct_points == 25.0


def test_docking_to_md_retention_distinguishes_all_categories():
    docking = _result(
        ["p1"],
        [
            ("p1", "SER70", "hbond", "H1", "OG", 2.8),
            ("p1", "GLU166", "saltbridge", "N1", "OE1", 3.5),
            ("p1", "LYS73", "hbond", "H2", "NZ", 2.9),
        ],
    )
    md = _result(
        ["f1", "f2", "f3", "f4"],
        [
            ("f1", "SER70", "hbond", "H1", "OG", 2.8),
            ("f2", "SER70", "hbond", "H2", "OG", 2.9),
            ("f3", "SER70", "hbond", "H3", "OG", 3.0),
            ("f1", "GLU166", "saltbridge", "N1", "OE1", 3.5),
            ("f2", "ASN132", "hbond", "H4", "OD1", 2.9),
        ],
    )

    frame = docking_md_retention(docking, md, retained_threshold_pct=50)
    categories = {
        (row.receptor_residue, row.interaction_type): row.category
        for row in frame.itertuples()
    }

    assert categories[("SER70", "hbond")] == "retained"
    assert categories[("GLU166", "saltbridge")] == "intermittent"
    assert categories[("LYS73", "hbond")] == "lost"
    assert categories[("ASN132", "hbond")] == "gained"


def test_analysis_context_validates_mode_and_time_step():
    with pytest.raises(ValueError, match="mode"):
        AnalysisContext(mode="trajectory")
    with pytest.raises(ValueError, match="time_step_ns"):
        AnalysisContext(mode="md", time_step_ns=0)
