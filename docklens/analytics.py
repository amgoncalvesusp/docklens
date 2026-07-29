"""Auditable analytical transforms for docking poses and saved MD frames.

Aggregate charts use one binary event per observation, receptor residue and
interaction type. Raw atom-pair rows remain available in ``RunResult.details``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Iterable, Tuple

import numpy as np
import pandas as pd

from .results import RunResult

if TYPE_CHECKING:
    from .observation_series import ObservationSeries


_MODES = frozenset({"docking", "md"})
_RETENTION_CATEGORIES = ("retained", "intermittent", "lost", "gained")


@dataclass(frozen=True)
class AnalysisContext:
    """Describe how ordered observations should be interpreted."""

    mode: str = "docking"
    time_step_ns: float = 1.0

    def __post_init__(self):
        if self.mode not in _MODES:
            raise ValueError("mode must be 'docking' or 'md'")
        if self.time_step_ns <= 0:
            raise ValueError("time_step_ns must be greater than zero")

    @property
    def observation_label(self) -> str:
        return "pose" if self.mode == "docking" else "frame"

    @property
    def aggregate_label(self) -> str:
        return "frequency" if self.mode == "docking" else "occupancy"


@dataclass(frozen=True)
class FingerprintCluster:
    """A deterministic similarity-threshold component and its representative."""

    cluster_id: int
    members: Tuple[str, ...]
    medoid: str
    mean_similarity: float


@dataclass(frozen=True)
class EpisodeStatistic:
    """Persistence statistics for one residue/interaction channel."""

    receptor_residue: str
    interaction_type: str
    observation_count: int
    total_observations: int
    occupancy_pct: float
    episode_count: int
    longest_episode_observations: int
    mean_episode_observations: float
    longest_episode_ns: float
    mean_episode_ns: float
    mean_distance_A: float | None
    is_key_residue: bool


def observation_ids(result: RunResult) -> Tuple[str, ...]:
    """Return stable observation identifiers, including empty poses/frames."""
    ordered = []
    seen = set()
    for pose_id in (
        [summary.pose_id for summary in result.summaries]
        + [detail.pose_id for detail in result.details]
    ):
        if pose_id not in seen:
            seen.add(pose_id)
            ordered.append(pose_id)
    return tuple(ordered)


def _presence_records(result: RunResult) -> Tuple[tuple, ...]:
    """Consolidate duplicate atom pairs into semantic chart observations."""
    unique = {
        (
            detail.pose_id,
            detail.receptor_residue,
            detail.interaction_type,
            bool(detail.is_key_residue),
        )
        for detail in result.details
    }
    return tuple(sorted(unique, key=lambda value: (value[0], value[1], value[2])))


def residue_type_prevalence(result: RunResult) -> pd.DataFrame:
    """Return prevalence per receptor residue and interaction type."""
    columns = [
        "receptor_residue",
        "interaction_type",
        "observation_count",
        "total_observations",
        "prevalence_pct",
        "is_key_residue",
    ]
    total = len(observation_ids(result))
    grouped = {}
    for pose_id, residue, kind, is_key in _presence_records(result):
        key = (residue, kind)
        entry = grouped.setdefault(key, {"observations": set(), "is_key": False})
        entry["observations"].add(pose_id)
        entry["is_key"] = entry["is_key"] or is_key
    rows = [
        {
            "receptor_residue": residue,
            "interaction_type": kind,
            "observation_count": len(entry["observations"]),
            "total_observations": total,
            "prevalence_pct": (
                100.0 * len(entry["observations"]) / total if total else 0.0
            ),
            "is_key_residue": entry["is_key"],
        }
        for (residue, kind), entry in sorted(grouped.items())
    ]
    return pd.DataFrame(rows, columns=columns)


def fingerprint_matrix(result: RunResult) -> pd.DataFrame:
    """Build a boolean observation × (residue, interaction type) matrix."""
    observations = observation_ids(result)
    features = tuple(
        sorted(
            {
                (detail.receptor_residue, detail.interaction_type)
                for detail in result.details
            }
        )
    )
    columns = pd.MultiIndex.from_tuples(
        features, names=("receptor_residue", "interaction_type")
    )
    matrix = pd.DataFrame(False, index=observations, columns=columns, dtype=bool)
    for pose_id, residue, kind, _is_key in _presence_records(result):
        matrix.loc[pose_id, (residue, kind)] = True
    matrix.index.name = "observation_id"
    return matrix


def fingerprint_similarity(matrix: pd.DataFrame) -> pd.DataFrame:
    """Calculate pairwise Jaccard/Tanimoto similarity for binary fingerprints."""
    names = tuple(str(value) for value in matrix.index)
    values = matrix.to_numpy(dtype=bool, copy=True)
    similarities = np.eye(len(values), dtype=float)
    for left in range(len(values)):
        for right in range(left + 1, len(values)):
            union = np.logical_or(values[left], values[right]).sum()
            intersection = np.logical_and(values[left], values[right]).sum()
            score = float(intersection / union) if union else 1.0
            similarities[left, right] = score
            similarities[right, left] = score
    return pd.DataFrame(similarities, index=names, columns=names)


def fingerprint_clusters(
    matrix: pd.DataFrame,
    *,
    threshold: float = 0.65,
    similarity: pd.DataFrame | None = None,
) -> Tuple[FingerprintCluster, ...]:
    """Group threshold-connected fingerprints and select deterministic medoids."""
    if not 0 <= threshold <= 1:
        raise ValueError("threshold must be between zero and one")
    similarity = (
        similarity.copy(deep=True)
        if similarity is not None
        else fingerprint_similarity(matrix)
    )
    if tuple(similarity.index) != tuple(str(value) for value in matrix.index):
        raise ValueError("similarity observations must match the matrix")
    names = tuple(similarity.index)
    order = {name: index for index, name in enumerate(names)}
    remaining = set(names)
    components = []
    while remaining:
        seed = min(remaining, key=order.get)
        queue = [seed]
        component = []
        remaining.remove(seed)
        while queue:
            current = queue.pop(0)
            component.append(current)
            neighbours = [
                candidate
                for candidate in remaining
                if similarity.loc[current, candidate] >= threshold
            ]
            for candidate in sorted(neighbours, key=order.get):
                remaining.remove(candidate)
                queue.append(candidate)
        components.append(tuple(sorted(component, key=order.get)))
    components.sort(key=lambda values: (-len(values), order[values[0]]))
    clusters = []
    for cluster_id, members in enumerate(components, 1):
        within = similarity.loc[list(members), list(members)]
        averages = within.mean(axis=1)
        best_score = float(averages.max())
        medoid = min(
            (name for name in members if float(averages[name]) == best_score),
            key=order.get,
        )
        clusters.append(
            FingerprintCluster(
                cluster_id=cluster_id,
                members=members,
                medoid=medoid,
                mean_similarity=best_score,
            )
        )
    return tuple(clusters)


def _consecutive_run_lengths(indices: Iterable[int]) -> Tuple[int, ...]:
    ordered = tuple(sorted(set(indices)))
    if not ordered:
        return ()
    runs = []
    length = 1
    for previous, current in zip(ordered, ordered[1:]):
        if current == previous + 1:
            length += 1
        else:
            runs.append(length)
            length = 1
    runs.append(length)
    return tuple(runs)


def episode_statistics(
    result: RunResult,
    context: AnalysisContext | None = None,
    *,
    series: ObservationSeries | None = None,
) -> Tuple[EpisodeStatistic, ...]:
    """Summarize saved-frame occupancy and continuous interaction episodes."""
    context = context or AnalysisContext(mode="md")
    observations = observation_ids(result)
    if series is not None:
        if series.mode != "md":
            raise ValueError("episode statistics require an MD observation series")
        if set(series.observation_ids) != set(observations):
            raise ValueError(
                "observation series IDs must match the analytical result"
            )
        observations = series.observation_ids
    index = {pose_id: position for position, pose_id in enumerate(observations)}
    grouped = {}
    for pose_id, residue, kind, is_key in _presence_records(result):
        entry = grouped.setdefault(
            (residue, kind), {"indices": [], "is_key": False}
        )
        entry["indices"].append(index[pose_id])
        entry["is_key"] = entry["is_key"] or is_key
    distances = {}
    for detail in result.details:
        if detail.distance_A is not None:
            distances.setdefault(
                (detail.receptor_residue, detail.interaction_type), []
            ).append(detail.distance_A)
    total = len(observations)
    metrics = []
    for (residue, kind), entry in sorted(grouped.items()):
        if series is None:
            runs = _consecutive_run_lengths(entry["indices"])
        else:
            active_ids = {
                observations[position] for position in set(entry["indices"])
            }
            runs = _series_run_lengths(active_ids, series)
        count = len(set(entry["indices"]))
        longest = max(runs, default=0)
        mean_run = float(np.mean(runs)) if runs else 0.0
        geometry = distances.get((residue, kind), ())
        metrics.append(
            EpisodeStatistic(
                receptor_residue=residue,
                interaction_type=kind,
                observation_count=count,
                total_observations=total,
                occupancy_pct=(100.0 * count / total if total else 0.0),
                episode_count=len(runs),
                longest_episode_observations=longest,
                mean_episode_observations=mean_run,
                longest_episode_ns=longest * context.time_step_ns,
                mean_episode_ns=mean_run * context.time_step_ns,
                mean_distance_A=(
                    float(np.mean(geometry)) if geometry else None
                ),
                is_key_residue=entry["is_key"],
            )
        )
    return tuple(metrics)


def _series_run_lengths(
    active_ids: set[str], series: ObservationSeries
) -> Tuple[int, ...]:
    """Measure active runs without joining replicas or missing frame indices."""
    runs = []
    current_length = 0
    previous = None
    for point in series.points:
        boundary = (
            previous is None
            or point.replica_id != previous.replica_id
            or (
                point.frame_index is not None
                and previous.frame_index is not None
                and point.frame_index != previous.frame_index + 1
            )
        )
        if boundary and current_length:
            runs.append(current_length)
            current_length = 0
        if point.observation_id in active_ids:
            current_length += 1
        elif current_length:
            runs.append(current_length)
            current_length = 0
        previous = point
    if current_length:
        runs.append(current_length)
    return tuple(runs)


def differential_prevalence(
    system_a: RunResult, system_b: RunResult
) -> pd.DataFrame:
    """Compare systems with independent observation denominators."""
    columns = [
        "receptor_residue",
        "interaction_type",
        "prevalence_a_pct",
        "prevalence_b_pct",
        "delta_pct_points",
        "is_key_residue",
    ]

    def indexed(result):
        frame = residue_type_prevalence(result)
        return {
            (row.receptor_residue, row.interaction_type): row
            for row in frame.itertuples()
        }

    left = indexed(system_a)
    right = indexed(system_b)
    rows = []
    for residue, kind in sorted(set(left) | set(right)):
        a_row = left.get((residue, kind))
        b_row = right.get((residue, kind))
        a_value = float(a_row.prevalence_pct) if a_row else 0.0
        b_value = float(b_row.prevalence_pct) if b_row else 0.0
        rows.append(
            {
                "receptor_residue": residue,
                "interaction_type": kind,
                "prevalence_a_pct": a_value,
                "prevalence_b_pct": b_value,
                "delta_pct_points": b_value - a_value,
                "is_key_residue": bool(
                    (a_row and a_row.is_key_residue)
                    or (b_row and b_row.is_key_residue)
                ),
            }
        )
    return pd.DataFrame(rows, columns=columns)


def docking_md_retention(
    docking: RunResult,
    md: RunResult,
    *,
    retained_threshold_pct: float = 50.0,
) -> pd.DataFrame:
    """Classify docking contacts against their saved-frame MD occupancy."""
    if not 0 <= retained_threshold_pct <= 100:
        raise ValueError("retained_threshold_pct must be between zero and 100")
    columns = [
        "receptor_residue",
        "interaction_type",
        "docking_prevalence_pct",
        "md_occupancy_pct",
        "category",
        "is_key_residue",
    ]
    docking_rows = {
        (row.receptor_residue, row.interaction_type): row
        for row in residue_type_prevalence(docking).itertuples()
    }
    md_rows = {
        (row.receptor_residue, row.interaction_type): row
        for row in residue_type_prevalence(md).itertuples()
    }
    rows = []
    for residue, kind in sorted(set(docking_rows) | set(md_rows)):
        docking_row = docking_rows.get((residue, kind))
        md_row = md_rows.get((residue, kind))
        docking_value = (
            float(docking_row.prevalence_pct) if docking_row else 0.0
        )
        md_value = float(md_row.prevalence_pct) if md_row else 0.0
        if docking_value == 0 and md_value > 0:
            category = "gained"
        elif docking_value > 0 and md_value == 0:
            category = "lost"
        elif docking_value > 0 and md_value >= retained_threshold_pct:
            category = "retained"
        else:
            category = "intermittent"
        rows.append(
            {
                "receptor_residue": residue,
                "interaction_type": kind,
                "docking_prevalence_pct": docking_value,
                "md_occupancy_pct": md_value,
                "category": category,
                "is_key_residue": bool(
                    (docking_row and docking_row.is_key_residue)
                    or (md_row and md_row.is_key_residue)
                ),
            }
        )
    frame = pd.DataFrame(rows, columns=columns)
    if not frame.empty:
        frame["category"] = pd.Categorical(
            frame["category"], categories=_RETENTION_CATEGORIES, ordered=True
        )
    return frame


def docking_score_metrics(result: RunResult) -> pd.DataFrame:
    """Pair docking scores with auditable per-pose interaction metrics."""
    presence = fingerprint_matrix(result)
    key_features = {
        (detail.receptor_residue, detail.interaction_type)
        for detail in result.details
        if detail.is_key_residue
    }
    rows = []
    for summary in result.summaries:
        vector = (
            presence.loc[summary.pose_id]
            if summary.pose_id in presence.index
            else pd.Series(dtype=bool)
        )
        active = {
            feature for feature, value in vector.items() if bool(value)
        }
        rows.append(
            {
                "pose_id": summary.pose_id,
                "ligand_id": summary.ligand_id,
                "pose": summary.pose,
                "docking_score": summary.docking_score,
                "distinct_residue_type_contacts": len(active),
                "distinct_residues": len({feature[0] for feature in active}),
                "key_residue_type_contacts": len(active & key_features),
                "key_residues": len(
                    {feature[0] for feature in active & key_features}
                ),
            }
        )
    return pd.DataFrame(rows)


__all__ = [
    "AnalysisContext",
    "EpisodeStatistic",
    "FingerprintCluster",
    "differential_prevalence",
    "docking_md_retention",
    "docking_score_metrics",
    "episode_statistics",
    "fingerprint_clusters",
    "fingerprint_matrix",
    "fingerprint_similarity",
    "observation_ids",
    "residue_type_prevalence",
]
