"""Deterministic interaction states for MD and pose families for docking.

The clustering model is trained on a disclosed, evenly spaced subset when the
number of observations exceeds ``max_training_observations``.  Every remaining
observation is assigned to its most similar training medoid, provided that it
meets the same similarity threshold; otherwise it is explicitly reported as an
``OUTLIER`` instead of being forced into an unsupported state.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TypeAlias

import numpy as np
import pandas as pd

from .analytics import (
    AnalysisContext,
    fingerprint_similarity,
)
from .observation_series import (
    ObservationSeries,
    default_observation_series,
)


Feature: TypeAlias = tuple[str, str]


@dataclass(frozen=True)
class StateAssignment:
    """Assignment of one pose/frame to a learned interaction state."""

    observation_id: str
    state_id: str
    similarity_to_representative: float
    is_training_observation: bool
    is_outlier: bool
    ordinal: int
    frame_index: int | None
    time_ns: float | None
    replica_id: str


@dataclass(frozen=True)
class InteractionState:
    """One MD interaction state or docking pose family."""

    state_id: str
    members: tuple[str, ...]
    representative: str
    characteristic_features: tuple[Feature, ...]
    population_count: int
    population_pct: float
    mean_similarity: float
    episode_count: int | None
    longest_dwell_observations: int | None
    mean_dwell_observations: float | None
    longest_dwell_ns: float | None
    mean_dwell_ns: float | None


@dataclass(frozen=True)
class InteractionStateAnalysis:
    """Auditable result of state/family discovery and full-data assignment."""

    context: AnalysisContext
    series: ObservationSeries
    states: tuple[InteractionState, ...]
    assignments: tuple[StateAssignment, ...]
    threshold: float
    method: str
    total_observations: int
    training_observations: int
    training_observation_ids: tuple[str, ...]
    sampled: bool
    sampling_method: str
    outlier_observations: tuple[str, ...]

    @property
    def mode(self) -> str:
        return self.context.mode

    @property
    def group_term(self) -> str:
        """Use scientifically correct terminology for the data modality."""
        return "pose family" if self.mode == "docking" else "interaction state"


def _normalise_matrix(matrix: pd.DataFrame) -> pd.DataFrame:
    if not isinstance(matrix, pd.DataFrame):
        raise TypeError("matrix must be a pandas DataFrame")
    normalised = matrix.copy(deep=True)
    normalised.index = tuple(str(value) for value in normalised.index)
    if not normalised.index.is_unique:
        raise ValueError("matrix observation IDs must be unique")
    if not isinstance(normalised.columns, pd.MultiIndex):
        raise ValueError("matrix columns must identify residue/type features")
    if normalised.columns.nlevels != 2:
        raise ValueError("matrix features must contain residue and interaction type")
    normalised.columns = pd.MultiIndex.from_tuples(
        tuple((str(left), str(right)) for left, right in normalised.columns),
        names=("receptor_residue", "interaction_type"),
    )
    return normalised.astype(bool)


def _ordered_matrix(
    matrix: pd.DataFrame,
    context: AnalysisContext,
    series: ObservationSeries | None,
) -> tuple[pd.DataFrame, ObservationSeries]:
    resolved = series or default_observation_series(
        matrix.index,
        mode=context.mode,
        time_step_ns=context.time_step_ns,
    )
    if resolved.mode != context.mode:
        raise ValueError("observation series mode must match the analysis context")
    matrix_ids = tuple(str(value) for value in matrix.index)
    if set(resolved.observation_ids) != set(matrix_ids):
        raise ValueError("observation series IDs must match the matrix")
    return matrix.loc[list(resolved.observation_ids)].copy(deep=True), resolved


def _training_indices(total: int, maximum: int) -> tuple[int, ...]:
    if maximum < 1:
        raise ValueError("max_training_observations must be at least one")
    if total <= maximum:
        return tuple(range(total))
    if maximum == 1:
        return (0,)
    return tuple(
        int(round(index * (total - 1) / (maximum - 1)))
        for index in range(maximum)
    )


def _jaccard(left: np.ndarray, right: np.ndarray) -> float:
    union = np.logical_or(left, right).sum()
    if not union:
        return 1.0
    return float(np.logical_and(left, right).sum() / union)


def _assign_states(
    matrix: pd.DataFrame,
    training_ids: tuple[str, ...],
    *,
    threshold: float,
) -> tuple[
    tuple[tuple[str, str, float], ...],
    tuple[tuple[str, tuple[str, ...], str], ...],
]:
    training = matrix.loc[list(training_ids)]
    similarity = fingerprint_similarity(training)
    clusters = _complete_link_clusters(similarity, threshold=threshold)
    cluster_records = tuple(
        (f"S{cluster_id}", members, medoid)
        for cluster_id, (members, medoid) in enumerate(clusters, 1)
    )
    medoid_vectors = tuple(
        (state_id, medoid, matrix.loc[medoid].to_numpy(dtype=bool, copy=True))
        for state_id, _members, medoid in cluster_records
    )
    assignments = []
    for observation_id, row in matrix.iterrows():
        vector = row.to_numpy(dtype=bool, copy=True)
        scores = tuple(
            (state_id, _jaccard(vector, medoid_vector))
            for state_id, _medoid, medoid_vector in medoid_vectors
        )
        best_state, best_score = max(
            scores,
            key=lambda item: (item[1], -int(item[0][1:])),
        )
        assigned_state = best_state if best_score >= threshold else "OUTLIER"
        assignments.append((str(observation_id), assigned_state, best_score))
    return tuple(assignments), cluster_records


def _complete_link_clusters(
    similarity: pd.DataFrame,
    *,
    threshold: float,
) -> tuple[tuple[tuple[str, ...], str], ...]:
    """Agglomerate only clusters whose every cross-pair meets the threshold."""
    clusters = [tuple([str(name)]) for name in similarity.index]
    while True:
        candidates = []
        for left_index, left in enumerate(clusters):
            for right_index in range(left_index + 1, len(clusters)):
                right = clusters[right_index]
                complete_link = min(
                    float(similarity.loc[left_name, right_name])
                    for left_name in left
                    for right_name in right
                )
                if complete_link >= threshold:
                    merged = tuple(sorted(left + right))
                    candidates.append(
                        (
                            -complete_link,
                            merged,
                            left_index,
                            right_index,
                        )
                    )
        if not candidates:
            break
        _score, merged, left_index, right_index = min(candidates)
        clusters = [
            cluster
            for index, cluster in enumerate(clusters)
            if index not in {left_index, right_index}
        ] + [merged]
    canonical = tuple(
        sorted(
            (tuple(sorted(cluster)) for cluster in clusters),
            key=lambda members: (-len(members), members),
        )
    )
    results = []
    for members in canonical:
        within = similarity.loc[list(members), list(members)]
        averages = within.mean(axis=1)
        best = float(averages.max())
        medoid = min(
            name for name in members if np.isclose(float(averages[name]), best)
        )
        results.append((members, medoid))
    return tuple(results)


def _episode_lengths(
    state_id: str,
    assignment_by_id: dict[str, str],
    series: ObservationSeries,
) -> tuple[int, ...]:
    grouped: dict[str, list] = {}
    for point in series.points:
        grouped.setdefault(point.replica_id, []).append(point)
    runs = []
    for points in grouped.values():
        current = 0
        previous = None
        for point in points:
            has_gap = (
                previous is not None
                and previous.frame_index is not None
                and point.frame_index is not None
                and point.frame_index - previous.frame_index != 1
            )
            if has_gap or assignment_by_id[point.observation_id] != state_id:
                if current:
                    runs.append(current)
                current = 0
            if assignment_by_id[point.observation_id] == state_id:
                current += 1
            previous = point
        if current:
            runs.append(current)
    return tuple(runs)


def _characteristic_features(
    matrix: pd.DataFrame,
    members: tuple[str, ...],
) -> tuple[Feature, ...]:
    if not members:
        return ()
    consensus = matrix.loc[list(members)].all(axis=0)
    return tuple(sorted(feature for feature, present in consensus.items() if present))


def interaction_state_analysis(
    matrix: pd.DataFrame,
    context: AnalysisContext | None = None,
    *,
    threshold: float = 0.65,
    max_training_observations: int = 300,
    series: ObservationSeries | None = None,
) -> InteractionStateAnalysis:
    """Discover MD states or docking pose families from binary fingerprints."""
    context = context or AnalysisContext(mode="md")
    if not 0 <= threshold <= 1:
        raise ValueError("threshold must be between zero and one")
    normalised = _normalise_matrix(matrix)
    ordered, resolved_series = _ordered_matrix(normalised, context, series)
    total = len(ordered)
    if total == 0:
        return InteractionStateAnalysis(
            context=context,
            series=resolved_series,
            states=(),
            assignments=(),
            threshold=threshold,
            method="complete-link threshold clustering",
            total_observations=0,
            training_observations=0,
            training_observation_ids=(),
            sampled=False,
            sampling_method="all observations",
            outlier_observations=(),
        )

    training_indices = _training_indices(total, max_training_observations)
    training_ids = tuple(str(ordered.index[index]) for index in training_indices)
    assigned, cluster_records = _assign_states(
        ordered,
        training_ids,
        threshold=threshold,
    )
    point_by_id = resolved_series.point_map()
    training_set = frozenset(training_ids)
    assignment_objects = tuple(
        StateAssignment(
            observation_id=observation_id,
            state_id=state_id,
            similarity_to_representative=similarity,
            is_training_observation=observation_id in training_set,
            is_outlier=state_id == "OUTLIER",
            ordinal=point_by_id[observation_id].ordinal,
            frame_index=point_by_id[observation_id].frame_index,
            time_ns=point_by_id[observation_id].time_ns,
            replica_id=point_by_id[observation_id].replica_id,
        )
        for observation_id, state_id, similarity in assigned
    )
    assignment_by_id = {
        assignment.observation_id: assignment.state_id
        for assignment in assignment_objects
    }
    similarity_by_id = {
        assignment.observation_id: assignment.similarity_to_representative
        for assignment in assignment_objects
    }
    effective_step = resolved_series.time_step_ns or context.time_step_ns
    outliers = tuple(
        sorted(
            assignment.observation_id
            for assignment in assignment_objects
            if assignment.is_outlier
        )
    )
    states = []
    for state_id, _training_members, medoid in cluster_records:
        members = tuple(sorted(
            assignment.observation_id
            for assignment in assignment_objects
            if assignment.state_id == state_id
        ))
        runs = (
            _episode_lengths(state_id, assignment_by_id, resolved_series)
            if context.mode == "md"
            else ()
        )
        longest = max(runs, default=0) if context.mode == "md" else None
        mean_run = float(np.mean(runs)) if runs else (
            0.0 if context.mode == "md" else None
        )
        states.append(
            InteractionState(
                state_id=state_id,
                members=members,
                representative=medoid,
                characteristic_features=_characteristic_features(
                    ordered,
                    members,
                ),
                population_count=len(members),
                population_pct=100.0 * len(members) / total,
                mean_similarity=float(
                    np.mean([similarity_by_id[member] for member in members])
                ),
                episode_count=(len(runs) if context.mode == "md" else None),
                longest_dwell_observations=longest,
                mean_dwell_observations=mean_run,
                longest_dwell_ns=(
                    longest * effective_step if longest is not None else None
                ),
                mean_dwell_ns=(
                    mean_run * effective_step if mean_run is not None else None
                ),
            )
        )
    if outliers:
        representative = outliers[0]
        representative_vector = ordered.loc[representative].to_numpy(
            dtype=bool,
            copy=True,
        )
        similarities = tuple(
            _jaccard(
                ordered.loc[member].to_numpy(dtype=bool, copy=True),
                representative_vector,
            )
            for member in outliers
        )
        runs = (
            _episode_lengths("OUTLIER", assignment_by_id, resolved_series)
            if context.mode == "md"
            else ()
        )
        longest = max(runs, default=0) if context.mode == "md" else None
        mean_run = float(np.mean(runs)) if runs else (
            0.0 if context.mode == "md" else None
        )
        states.append(
            InteractionState(
                state_id="OUTLIER",
                members=outliers,
                representative=representative,
                characteristic_features=_characteristic_features(
                    ordered,
                    outliers,
                ),
                population_count=len(outliers),
                population_pct=100.0 * len(outliers) / total,
                mean_similarity=float(np.mean(similarities)),
                episode_count=(len(runs) if context.mode == "md" else None),
                longest_dwell_observations=longest,
                mean_dwell_observations=mean_run,
                longest_dwell_ns=(
                    longest * effective_step if longest is not None else None
                ),
                mean_dwell_ns=(
                    mean_run * effective_step if mean_run is not None else None
                ),
            )
        )
    return InteractionStateAnalysis(
        context=context,
        series=resolved_series,
        states=tuple(states),
        assignments=assignment_objects,
        threshold=threshold,
        method="complete-link threshold clustering",
        total_observations=total,
        training_observations=len(training_ids),
        training_observation_ids=training_ids,
        sampled=len(training_ids) < total,
        sampling_method=(
            "evenly spaced observations"
            if len(training_ids) < total
            else "all observations"
        ),
        outlier_observations=outliers,
    )


def state_assignment_frame(analysis: InteractionStateAnalysis) -> pd.DataFrame:
    """Return ordered, fully disclosed state/family assignments."""
    columns = [
        "observation_id",
        "state_id",
        "similarity_to_representative",
        "is_training_observation",
        "is_outlier",
        "ordinal",
        "frame_index",
        "time_ns",
        "replica_id",
    ]
    return pd.DataFrame(
        [
            {
                field: getattr(assignment, field)
                for field in columns
            }
            for assignment in analysis.assignments
        ],
        columns=columns,
    )


def state_summary_frame(analysis: InteractionStateAnalysis) -> pd.DataFrame:
    """Return one row per state/family, using modality-safe terminology."""
    columns = [
        "state_id",
        "group_term",
        "representative",
        "population_count",
        "population_pct",
        "mean_similarity",
        "characteristic_features",
        "episode_count",
        "longest_dwell_observations",
        "mean_dwell_observations",
        "longest_dwell_ns",
        "mean_dwell_ns",
    ]
    return pd.DataFrame(
        [
            {
                "state_id": state.state_id,
                "group_term": analysis.group_term,
                "representative": state.representative,
                "population_count": state.population_count,
                "population_pct": state.population_pct,
                "mean_similarity": state.mean_similarity,
                "characteristic_features": state.characteristic_features,
                "episode_count": state.episode_count,
                "longest_dwell_observations": state.longest_dwell_observations,
                "mean_dwell_observations": state.mean_dwell_observations,
                "longest_dwell_ns": state.longest_dwell_ns,
                "mean_dwell_ns": state.mean_dwell_ns,
            }
            for state in analysis.states
        ],
        columns=columns,
    )


def state_transition_frame(
    analysis: InteractionStateAnalysis,
    *,
    lag: int = 1,
) -> pd.DataFrame:
    """Return observed MD transitions without crossing replicas or frame gaps."""
    if analysis.mode != "md":
        raise ValueError("temporal state transitions are available only for MD")
    assignment_by_id = {
        assignment.observation_id: assignment.state_id
        for assignment in analysis.assignments
    }
    counts: dict[tuple[str, str], int] = {}
    for left, right in analysis.series.transition_pairs(lag=lag):
        edge = (assignment_by_id[left], assignment_by_id[right])
        counts[edge] = counts.get(edge, 0) + 1
    totals: dict[str, int] = {}
    for (source, _target), count in counts.items():
        totals[source] = totals.get(source, 0) + count
    order = {
        state.state_id: index for index, state in enumerate(analysis.states)
    }
    order["OUTLIER"] = len(order)
    rows = [
        {
            "from_state": source,
            "to_state": target,
            "transition_count": count,
            "transition_probability_pct": 100.0 * count / totals[source],
            "lag_observations": lag,
        }
        for (source, target), count in sorted(
            counts.items(),
            key=lambda item: (
                order.get(item[0][0], len(order) + 1),
                order.get(item[0][1], len(order) + 1),
            ),
        )
    ]
    return pd.DataFrame(
        rows,
        columns=[
            "from_state",
            "to_state",
            "transition_count",
            "transition_probability_pct",
            "lag_observations",
        ],
    )


__all__ = [
    "InteractionState",
    "InteractionStateAnalysis",
    "StateAssignment",
    "interaction_state_analysis",
    "state_assignment_frame",
    "state_summary_frame",
    "state_transition_frame",
]
